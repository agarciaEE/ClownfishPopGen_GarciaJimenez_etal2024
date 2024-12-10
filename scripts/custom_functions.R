##########################
### REQUIRED LIBRARIES ###
##########################

library(ggh4x)
library(ggtext)
library(ggplot2)
library("ggsci")
library(tidyverse)
library(ggExtra)
library(patchwork)
library(ggpattern)
library(ggpubr)
library(dplyr)
library(tidyr)
library(sf)
library(ggnewscale)
library(CausalImpact)
library(fdANOVA)
library(fda)
library(ggrepel) 
library(sjPlot)
library(gtsummary)
library(gt)
library(mgcv)         
library(lme4)         
library(nlme)
library(broom.mixed)
library(phytools)


########################
### CUSTOM FUNCTIONS ###
########################

##### functions for msmc2 post analysis
#######################################

##### Harmonic mean
harmonic_mean <- function(x, w = NULL, na.rm = TRUE) {
  
  n <- length(x)
  
  if (n == 0) {
    return(NA)
  }
  
  # If there are any zeroes in x, the harmonic mean is undefined
  if (any(x == 0, na.rm = TRUE)) {
    stop("Harmonic mean is undefined if any element is zero.")
  }
  
  if (!is.null(w)) {
    # Ensure the length of weights matches the length of x
    if (length(w) != n) {
      stop("Length of weights vector must be equal to the length of the data vector.")
    }
    
    # Compute weighted harmonic mean
    whm <- sum(w, na.rm = na.rm) / sum(w / x, na.rm = na.rm)
    return(whm)
    
  } else {
    # Compute harmonic mean without weights
    hm <- n / sum(1 / x, na.rm = na.rm)
    return(hm)
  }
}

##### msmc summary in intervals

summary_Ne_in_interval <- function(left, right, df) {
  
  # Filter species that overlap with the current interval
  overlapping_reconstructions <- df %>%
    filter(left_time_boundary < right & right_time_boundary > left)
  
  if (nrow(overlapping_reconstructions) == 0) {
    return(data.frame(min_Ne = NA, median_Ne = NA, mean_Ne = NA, sd_Ne = NA, max_Ne = NA))
  }
  
  # Calculate mean, sd, min, and max for overlapping species
  return(data.frame(
    min_Ne = min(overlapping_reconstructions$Ne, na.rm = TRUE),
    median_Ne = median(overlapping_reconstructions$Ne, na.rm = TRUE),
    mean_Ne = mean(overlapping_reconstructions$Ne, na.rm = TRUE),
    sd_Ne = sd(overlapping_reconstructions$Ne, na.rm = TRUE),
    max_Ne = max(overlapping_reconstructions$Ne, na.rm = TRUE)
  ))
}

##### Average msmc reconstruction

summary_by_intervals <- function(df, time_points = NULL, logarithmic = TRUE, num_time_points = 25) {
  
  # get time points
  if (is.null(time_points)){
    max_time <- max(df$time, na.rm = TRUE)
    min_time <- min(df$time, na.rm = TRUE)
    if (logarithmic){
      time_points <- exp(seq(log(min_time), log(max_time), length.out = num_time_points))
    } else {
      time_points <- seq(round(min_time/1E+6)*1E+6, round(max_time/1E+6)*1E+6, length.out = num_time_points)
    }
  }

  # Create non-overlapping time intervals based on unique time points
  intervals <- data.frame(
    left_time_boundary = head(time_points, -1),
    right_time_boundary = tail(time_points, -1)
  )
  
  # Apply the function to each interval
  stats <- mapply(summary_Ne_in_interval, 
                  intervals$left_time_boundary, 
                  intervals$right_time_boundary, 
                  MoreArgs = list(df = df), 
                  SIMPLIFY = FALSE)
  
  # Combine the results into the intervals data frame
  stats_df <- do.call(rbind, stats)
  intervals <- cbind(intervals, stats_df)
  
  return(intervals)
}

eval_model <- function(model) {
  
  # Extract the residuals and fitted values from the model
  residuals <- residuals(model)
  fitted_values <- fitted(model)
  
  # Perform Shapiro-Wilk test
  sw_residuals <- if(length(residuals) > 5000) { sample(residuals, 4999) } else { residuals }
  sw_test <- shapiro.test(sw_residuals)
  
  # Set up the 2x2 plotting area
  par(mfrow = c(2, 2), mar = c(5,5,5,5))
  
  # Plot 1: Residuals vs. Fitted Values with Scale-Location inset
  plot(fitted_values, residuals, pch = 19,
       xlab = "Fitted Values",
       ylab = "Residuals",
       main = "Residuals vs. Fitted Values")
  abline(h = 0, col = "red")
  
  # Add Scale-Location plot as an inset
  points(fitted_values, sqrt(abs(residuals)), pch = 19, col = alpha("red", 0.2), cex = 0.5, 
         xlab = "", ylab = "",
         main = "Scale-Location",
         cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)
  axis(4)
  mtext(expression(sqrt("|Residuals|")), side = 4, line = 2, cex = 0.7)
  
  # Plot 2: Normal Q-Q Plot
  qqnorm(residuals, main = "Normal Q-Q Plot of Residuals")
  qqline(residuals, col = "red")
  
  # Add Shapiro-Wilk test results to Q-Q plot
  sw_text <- sprintf("Shapiro-Wilk test:\nW = %.4f, p = %.4f", 
                     sw_test$statistic, sw_test$p.value)
  mtext(sw_text, side = 3, line = -2, adj = 0.05, cex = 0.7)
  
  # Plot 3: Histogram of Residuals
  hist(residuals, breaks = 30, main = "Histogram of Residuals",
       xlab = "Residuals", col = "lightblue")
  
  # Plot 4: Autocorrelation of Residuals
  acf(residuals, main = "Autocorrelation of Residuals")
}

# Function to extract model info: logLik, AIC, df
extract_model_info <- function(model) {
  logLik_val <- logLik(model) # Log-likelihood
  aic_val <- AIC(model)       # AIC
  df_val <- attr(logLik_val, "df") # Degrees of freedom
  aicc <- aic_val + (2 * df_val * (df_val + 1)) / (nobs(model) - df_val - 1)
  
  return(c(logLik = as.numeric(logLik_val), AIC = aic_val, df = df_val, AICc = aicc))
}


## EEMS

plot_log_posterior <- function (mcmcpaths) {
  
  message("Generate posterior probability trace. ", "See plots$pilog01.")
  rleid <- function(x) {
    r <- rle(x)
    rep(seq_along(r$lengths), r$lengths)
  }
  pl_df <- NULL
  for (path in mcmcpaths) {
    demes <- as.integer(gsub("ndemes", "", strsplit(path, "-")[[1]][3]))
    chain <- as.integer(gsub("chain", "", strsplit(path, "-")[[1]][4]))
    pl <- reemsplots2:::read_matrix(file.path(path, "mcmcpilogl.txt"))
    pl <- cbind(pl, chain, demes)
    pl_df <- dplyr::bind_rows(pl_df, as.data.frame(pl) %>% dplyr::mutate(path))
  }
  pl_df <- pl_df %>% setNames(c("pi", "logl", "chain", "demes", "path")) %>% 
    dplyr::mutate(mcmcpath = factor(rleid(path))) %>% group_by(mcmcpath) %>% 
    dplyr::mutate(iter = row_number(), pilogl = pi + logl)
  pl_df$demes <- factor(pl_df$demes, levels = sort(as.integer(unique(pl_df$demes))))
  pl_df$chain <- factor(pl_df$chain, levels = sort(as.integer(unique(pl_df$chain))))
  ggplot(pl_df, aes(x = iter, y = pilogl, group = path, color = chain)) + 
    geom_path() + labs(x = "MCMC iteration  (after burn-in and thinning)", 
                       y = "log posterior") + 
    theme_bw() + theme(panel.grid.minor = element_blank(), 
                       panel.grid.major.x = element_blank(),
                       text = element_text(size = 16),
                       strip.background = element_rect(fill = "black"),
                       strip.text = element_text(color = "white", size = 16)) +
    facet_wrap(.~demes)
}

# get max geodesic distance of a dataframe
maxDist <- function(df){
  df <- expand.grid(range(df[,1]), range(df[,2]))
  comb <- combn(1:4, 2)
  maxd <- max(apply(comb, 2, function(i) oce::geodDist(df[i[1],1], df[i[1],2],
                                                       df[i[2],1], df[i[2],2])), na.rm = T) 
  return(maxd)
}

# modified make_eems_plots function
make_eems_plots2 <- function (mcmcpath, longlat = TRUE, dpi = 250, add_grid = FALSE, 
                              col_grid = "#BBBBBB", add_demes = FALSE, col_demes = "#000000", 
                              add_outline = FALSE, col_outline = "#FFFFFF", eems_colors = NULL, 
                              prob_level = 0.9, m_colscale = NULL, q_colscale = NULL, add_abline = FALSE) 
{
  reemsplots2:::check_mcmcpath_contents(mcmcpath)
  func_params <- list(add_grid = add_grid, add_demes = add_demes, 
                      add_outline = add_outline, eems_colors = eems_colors, 
                      prob_level = prob_level, col_grid = col_grid, col_demes = col_demes, 
                      col_outline = col_outline, m_colscale = m_colscale, q_colscale = q_colscale, 
                      add_abline = add_abline)
  plot_params <- reemsplots2:::check_plot_params(func_params)
  dimns <- reemsplots2:::read_dimns(mcmcpath[1], longlat, dpi)
  plots <- list()
  p <- reemsplots2:::eems_contours(mcmcpath, dimns, longlat, plot_params, 
                                   is_mrates = TRUE)
  plots$mrates01 <- p[[1]]
  plots$mrates02 <- p[[2]]
  p <- reemsplots2:::eems_contours(mcmcpath, dimns, longlat, plot_params, 
                                   is_mrates = FALSE)
  plots$qrates01 <- p[[1]]
  plots$qrates02 <- p[[2]]
  plots
}

##### modified functions from plotmaps package
##############################################
# get Year based on IBD segment legnths (from Al-Alsadi 2018)
ibdL2nG <- function(l1, l2, gT = NULL) {
  
  nG <- 300/4 * (1/l1 + 1/l2)
  
  if (!is.null(gT)) nG * gT else nG
}

generate_rounded_breaks <- function(limits) {
  
  range <- diff(limits)
  magnitude <- 10^floor(log10(range))
  
  # Determine the number of decimal places needed
  max_decimals <- ifelse(magnitude > 1, 0, log10(magnitude) - 1)
  # Try different numbers of breaks
  for (n_breaks in c(4, 5, 6)) {
    step <- range / (n_breaks - 1)
    rounded_step <- round(step / magnitude) * magnitude
    
    while (rounded_step == 0) {
      # If rounded_step is 0, use a smaller step size
      magnitude <- magnitude / 10
      rounded_step <- round(step / magnitude) * magnitude
    }
    
    # Determine the number of decimal places for rounding
    digits <- max(0, -floor(log10(rounded_step)))
    
    # Adjust start to center the breaks
    start <- floor(limits[1] / rounded_step) * rounded_step
    
    breaks <- seq(start, by = rounded_step, length.out = n_breaks)
    
    # Check if breaks cover the range and are within limits
    if (max(breaks) >= round(limits[2], digits) && min(breaks) <= round(limits[1], digits)) {
      breaks <- unique(round(breaks, digits = max(digits, max_decimals)))
      
      return(breaks)
    }
  }
  # If no suitable breaks found, return 5 equally spaced breaks
  return(round(seq(limits[1], limits[2], length.out = 5), digits = max_decimals))
}

add_pts <- function (g, col_pts = "black", fill_pts = "black", const_size = F) 
{
  tbl <- table(g$ipmap)
  ind <- as.numeric(names(tbl))
  sizes <- as.vector(tbl)
  df <- data.frame(x = g$demes[ind, 1], y = g$demes[ind, 2], 
                   sizes = sizes)
  if (const_size) {
    pts <- geom_point(aes(x = x, y = y), data = df, color = col_pts, 
                      fill = fill_pts, size = 1.5)
  }
  else {
    pts <- geom_point(aes(x = x, y = y, size = sizes), data = df, 
                      colour = col_pts, fill = fill_pts, pch = 21, show.legend = FALSE)
  }
}

exp_color_gradient <- function(colors, n = 100) {
  # Generate exponential sequence
  exp_seq <- exp(seq(0, log(n), length.out = n)) / exp(log(n))
  
  # Interpolate colors
  colorRampPalette(colors)(n)[round(exp_seq * (n-1) + 1)]
}

get_dist_metric <- function (mcmcpath) 
{
  dist_metric <- "euclidean"
  lines <- tolower(readLines(file.path(mcmcpath, "mapsrun.txt")))
  for (line in lines) {
    if (grepl("\\s*distance\\s*=\\s*", line)) 
      dist_metric <- gsub("\\s*distance\\s*=\\s*(\\w+)", 
                          "\\1", line)
  }
  if (dist_metric != "euclidean" && dist_metric != "greatcirc") 
    stop("mapsrun.txt should specify `euclidean` or `greatcirc` distance.")
  dist_metric
}

compute_summary_statistic <- function (params, dimns, lower_quantile = 0.025, upper_quantile = 0.975) 
{
  rslts <- plotmaps:::compute_rates_each_pixel(params, dimns)
  
  mean.rate = NA
  med.rate = NA
  Pgt = Plw = NA
  
  if (params$plot.mean) {
    # Calculate mean rates across iterations
    mean.rate <- apply(rslts, c(2, 3), mean)
    var <- apply(rslts, c(2, 3), function(x) sd(x))
  }
  if (params$plot.median | params$plot.sign) {
    # Calculate median rates across iterations
    med.rate <- apply(rslts, c(2, 3), median)
    var <- apply(rslts, c(2, 3), function(x) IQR(x))
  }
  # Compute Confidence intervals
  # Compute the lower and upper bounds of the CI for the differences
  lower_ci <- apply(rslts, c(2, 3), function(x) quantile(x, probs = lower_quantile))
  upper_ci <- apply(rslts, c(2, 3), function(x) quantile(x, probs = upper_quantile))
  
  # Compute probability of rate being higher than the overall mean rate
  # Calculate the base rate
  base.rate <- apply(rslts, 1, median)
  Pgt <- array(dim = dim(rslts))
  Plw <- array(dim = dim(rslts))
  for (i in seq_along(base.rate)) {
    # Calculate the proportion of values greater than the corresponding base.rate value
    Pgt[i, , ] <- apply(rslts[i,,], 2, function(x) (x > base.rate[i]))
    Plw[i, , ] <- apply(rslts[i,,], 2, function(x) (x < base.rate[i]))
  }
  Pgt <- apply(Pgt, c(2, 3), mean)
  Plw <- apply(Plw, c(2, 3), mean)
  # same scaling from original plotmaps function
  if (!params$plot.difference) {
    mean.rate <- 10^mean.rate
    med.rate <- 10^med.rate
    lower_ci <- 10^lower_ci
    upper_ci <- 10^upper_ci
    var <- 10^var
    if (!params$is.mrates) {
      mean.rate <- 1/(2 * mean.rate)
      med.rate <- 1/(2 * med.rate)
      var <- 1/(2 * var)
      lower.ci <- lower_ci
      P_gt <- Pgt
      lower_ci <- 1/(2 * upper_ci) # values are inverted for N
      upper_ci <- 1/(2 * lower.ci) # values are inverted for N
      Pgt <- Plw # values are inverted for N
      Plw <- P_gt # values are inverted for N
    }
  }
  
  if (params$correct.by.area){
    l <- plotmaps:::compute_scaling(params$mcmcpath[1])
    
    if (params$is.mrates & !params$plot.difference) {
      mean.rate <- sqrt(mean.rate * l$m.scalingfactor)
      med.rate <- sqrt(med.rate * l$m.scalingfactor)
      lower_ci <- sqrt(lower_ci * l$m.scalingfactor)
      upper_ci <- sqrt(upper_ci * l$m.scalingfactor)
    }
    else if (!params$is.mrates & !params$plot.difference) {
      mean.rate <- mean.rate * l$N.scalingfactor
      med.rate <- med.rate * l$N.scalingfactor
      lower_ci <- sqrt(lower_ci * l$m.scalingfactor)
      upper_ci <- sqrt(upper_ci * l$m.scalingfactor)
    }
  }
  
  return(list(avg = mean.rate, med = med.rate, var = var,
              upper_ci = upper_ci, lower_ci = lower_ci, 
              Pgt = Pgt, Plw = Plw))
}

get_title <- function (params) 
{
  m <- expression(m)
  N <- expression(N)
  diff.m <- bquote(paste(log10, bgroup("(", frac(m^"'", m), 
                                       ")")))
  diff.N <- bquote(paste(log10, bgroup("(", frac(N^"'", N), 
                                       ")")))
  sigma <- expression(sigma)
  D <- expression("D"[e])
  diff.sigma <- bquote(paste(log10, bgroup("(", frac(sigma^"'", 
                                                     sigma), ")")))
  diff.D <- bquote(paste(log10, bgroup("(", frac("D"[e]^"'", 
                                                 "D"[e]), ")")))
  if (!params$plot.difference) {
    if (params$correct.by.area) {
      if (params$is.mrates) {
        return(sigma)
      }
      else {
        return(D)
      }
    }
    else {
      if (params$is.mrates) {
        return(m)
      }
      else {
        return(N)
      }
    }
  }
  else {
    if (params$correct.by.area) {
      if (params$is.mrates) {
        return(diff.sigma)
      }
      else {
        return(diff.D)
      }
    }
    else {
      if (params$is.mrates) {
        return(diff.m)
      }
      else {
        return(diff.N)
      }
    }
  }
}

add_contours <- function (params, dimns, g, summary_stats = NULL) 
{
  x <- dimns$marks[, 1]
  y <- dimns$marks[, 2]
  
  n_runs <- length(params$mcmcpath)
  
  if (is.null(summary_stats)){
    summary_stats <- compute_summary_statistic(params, dimns)
  }
  
  if (params$plot.mean) {
    summary_stat <- summary_stats$avg
  } else {
    summary_stat <- summary_stats$med
  }
  
  base.rate <- mean(summary_stat)
  
  df <- data.frame(x = x, y = y, ss = c(summary_stat), 
                   var = c(summary_stats$var),
                   upper.ci = c(summary_stats$upper_ci), 
                   lower.ci = c(summary_stats$lower_ci),
                   filter = dimns$filter)
  print(summary(df$ss))
  df <- df[df$filter, ]
  
  legend.title <- get_title(params)
  trans <- plotmaps:::get_trans(params)
  
  if (params$set.range) {
    if (params$is.mrates) {
      limits <- params$m.limits
    }
    else {
      limits <- params$N.limits
    }
    if (min(df$ss, na.rm = TRUE) < limits[1]){
      warning("Lower limit is higher than the minimum value: ", min(df$ss, na.rm = TRUE))
    }
    if (max(df$ss, na.rm = TRUE) > limits[2]){
      warning("Upper limit is lower than the maximum value: ", max(df$ss, na.rm = TRUE))
    }
  } else {
    limits <- range(df$ss)
  }
  # if values vary in folds use log10 trans for plotting
  magnitude <- diff(log10(limits))
  if (magnitude > 2) { trans = "log10" }
  
  # get breaks and labels
  if (trans == "log10" ) { 
    if (min(df$ss) == 0) {
      # Find the smallest non-zero value
      min_nonzero <- min(df$ss[df$ss > 0], na.rm = TRUE)
      # Set min.val to one order of magnitude below the smallest non-zero value
      min.val <- 10^(floor(log10(min_nonzero)) - 1)
    } else {
      # Set min.val to one order of magnitude below the minimum value
      min.val <- 10^(floor(log10(min(df$ss))) - 1)
    }
    limits = ifelse(limits == 0, min.val, limits) 
    breaks = 10^seq(floor(log10(limits[1])), ceiling(log10(limits[2])))
    labels = format(breaks, digits = 2, scientific = TRUE)
    labels[1] = 0
  } else {
    breaks <- generate_rounded_breaks(limits)
    if (any(breaks > 1e+3) |  any(breaks < 1e-2)){
      labels = format(breaks, digits = 2, scientific = TRUE)
    } else {
      labels = breaks
    }
  }
  
  if (is.null(params$col.palette)) {
    if (params$plot.sign) {
      col.palette <- plotmaps:::default_eems_colors() 
    } else {
      if (params$is.mrates) {
        col.palette <- colorRampPalette(c("#F7FBFF", "#08306B"))(n=100)
      }
      else {
        col.palette <- colorRampPalette(c("#F7FBFF", "#00441B"))(n=100)
      }
    }
  } else {  col.palette <- params$col.palette }
  
  if (params$plot.sign) {
    
    probs <- (summary_stats$Pgt - summary_stats$Plw + 1)/2
    probs[probs < 0] <- 0
    probs[probs > 1] <- 1
    
    sign <- data.frame(x = x, y = y, 
                       P = c(probs), 
                       filter = dimns$filter)
    
    sign <- sign[sign$filter, ]
    
    breaks <- sort(c(params$alpha, 1 - params$alpha))
    if (params$is.mrates) {
      r <- "m"
    } else { r <- "N" }
    
    labels <- c(bquote("P(" ~ .(r) ~ "<" ~ bar(.(r)) ~ ") =" ~ .(max(breaks))),
                bquote("P(" ~ .(r) ~ ">" ~ bar(.(r)) ~ ") =" ~ .(max(breaks))), expression())
    
    g <-  g + geom_raster(data = sign, aes(x = x, y = y, fill = P), 
                          alpha = 1) + 
      scale_fill_gradientn(colours = col.palette, 
                           name = "", na.value = "lightgray", 
                           limits = c(0, 1), breaks = breaks, labels = labels) +
      geom_contour(data = sign, aes(x = x, y = y, z = P), breaks =  breaks, color = "white") + 
      theme(legend.key.width = unit(0.75, "cm"),
            legend.text = element_text(size = 15),
            legend.key.height = unit(2.25, "cm"),
            legend.title = element_text(size = 15)) + coord_fixed() 
    
  } else {
    if (params$plot.ci == "upper") {
      # plot upper CI
      g <- g + geom_raster(data = df, aes(x = x, y = y, fill = upper.ci), 
                           alpha = 1) + 
        scale_fill_gradientn(colours = col.palette, 
                             name = legend.title, trans = trans, na.value = "lightgray") +
        theme(legend.key.width = unit(0.75, "cm"),
              legend.text = element_text(size = 15),
              legend.key.height = unit(2.25, "cm"),
              legend.title = element_text(size = 15)) + coord_fixed() 
    } else if (params$plot.ci == "lower") {
      # plot lower CI
      g <- g + geom_raster(data = df, aes(x = x, y = y, fill = lower.ci), 
                           alpha = 1) + 
        scale_fill_gradientn(colours = col.palette, 
                             name = legend.title, trans = trans, na.value = "lightgray") +
        theme(legend.key.width = unit(0.75, "cm"),
              legend.text = element_text(size = 15),
              legend.key.height = unit(2.25, "cm"),
              legend.title = element_text(size = 15)) + coord_fixed() 
      
    } else if (params$plot.ci == "var") {
      if (params$is.mrates) { r <- "m" } else { r <- "N" }
      # plot range of uncertainty (CI range per pixel)
      var.legend.title <- ifelse(params$plot.mean, c(bquote(sigma[.(r)]), expression()), c(bquote("IQR"[.(r)]), expression()))
      g <- g + geom_raster(data = df, aes(x = x, y = y, fill = var), 
                           alpha = 1) + 
        scale_fill_gradientn(colours =  rev(plotmaps:::default_eems_colors()),
                             name = var.legend.title, trans = "identity", na.value = "lightgray") +
        theme(legend.key.width = unit(0.75, "cm"),
              legend.text = element_text(size = 15),
              legend.key.height = unit(2.25, "cm"),
              legend.title = element_text(size = 15)) + coord_fixed() 
    }
    else {
      g <- g + geom_raster(data = df, aes(x = x, y = y, fill = ss), 
                           alpha = 1) + 
        geom_contour(data = df, aes(x = x, y = y, z = ss), breaks = breaks, color = "white") + 
        metR::geom_text_contour(data = df, aes(x = x, y = y, z = ss), 
                                breaks = breaks, stroke = 0.2, skip = 0, check_overlap = TRUE) +
        scale_fill_gradientn(colours = col.palette, limits = limits, breaks = breaks, labels = labels,
                             name = legend.title, trans = trans, na.value = "white") +
        theme(legend.key.width = unit(0.75, "cm"),
              legend.text = element_text(size = 15),
              legend.key.height = unit(2.25, "cm"),
              legend.title = element_text(size = 15)) + coord_fixed() 
    }
  } 
  
  graph <- read_graph(params$mcmcpath[1], params$longlat)
  if (params$add.graph) {
    g <- g + plotmaps:::add_graph(graph, color = params$col.graph)
  }
  if (params$add.pts) {
    g <- g + add_pts(graph, col_pts = params$col.pts, fill_pts = params$fill.pts, const_size = params$const.size)
  }
  return(g)
}

plot_contour <- function (params, summary_stats = NULL) 
{
  dimns <- plotmaps:::read_dimns(params$mcmcpath[1], params$longlat)
  g <- make_base(dimns, params)
  g <- add_contours(params, dimns, g, summary_stats = summary_stats)
  if (params$add.countries) {
    g <- add_map(dimns, params, g)
  }
  filename <- params$outpath
  if (params$is.mrates) {
    filename <- paste0(filename, "/mrates-")
  }
  else {
    filename <- paste0(filename, "/Nsizes-")
  }
  if (params$plot.mean) {
    if (params$plot.sign) {
      filename <- paste0(filename, "mean-sign.pdf")
    }
    else {
      filename <- paste0(filename, "mean", params$plot.ci, ".pdf")
    }
  }
  else {
    if (params$plot.sign) {
      filename <- paste0(filename, "median-sign.pdf")
    }
    else {
      filename <- paste0(filename, "median", params$plot.ci, ".pdf")
    }
  }
  if (params$save.plot) { ggsave(filename, width = params$width, height = params$height) }
  return(g)
}

plot_trace <- function (mcmcpaths, outpath) 
{
  
  message("Generate posterior probability trace.")
  rleid <- function(x) {
    r <- rle(x)
    rep(seq_along(r$lengths), r$lengths)
  }
  pl_df <- NULL
  for (path in mcmcpaths) {
    demes <- as.integer(gsub("ndemes", "", grep("ndemes", strsplit(path, "-")[[1]], value = TRUE)))
    if (length(demes) == 0 || is.na(demes)){ demes = ""}
    chain <-  as.integer(gsub("chain", "", grep("chain", strsplit(path, "-")[[1]], value = TRUE)))
    if (length(chain) == 0 || is.na(chain)){ chain = which(path %in% mcmcpaths)}
    stopifnot(file.exists(path))
    pl <- matrix(scan(file.path(path, "mcmcpilogl.txt"), what = numeric(), quiet = TRUE), ncol = 2, byrow = TRUE)
    pl <- cbind(pl, chain, demes)
    pl_df <- dplyr::bind_rows(pl_df, as.data.frame(pl) %>% dplyr::mutate(path))
  }
  pl_df <- pl_df %>% setNames(c("pi", "logl", "chain", "demes", "path")) %>% 
    dplyr::mutate(mcmcpath = factor(rleid(path))) %>% group_by(mcmcpath) %>% 
    dplyr::mutate(iter = row_number(), pilogl = pi + logl)
  pl_df$demes <- factor(pl_df$demes, levels = sort(as.integer(unique(pl_df$demes))))
  pl_df$chain <- factor(pl_df$chain, levels = sort(as.integer(unique(pl_df$chain))))
  g <- ggplot(pl_df, aes(x = iter, y = pilogl, group = path, color = chain)) + 
    geom_path() + labs(x = "MCMC iteration  (after burn-in and thinning)", 
                       y = "log posterior") + 
    scale_color_brewer() +
    theme_bw() + theme(panel.grid.minor = element_blank(), 
                       panel.grid.major.x = element_blank(),
                       plot.background = element_rect(color = "white"),
                       text = element_text(size = 16),
                       strip.background = element_rect(fill = "black"),
                       strip.text = element_text(color = "white", size = 16)) +
    facet_wrap(.~demes)
  if (!missing(outpath)) ggsave(paste0(outpath, "/logll_trace.pdf"), g)
  return(g)
}

plot_fit_data <- function (mcmcpath, outpath, longlat) 
{
  
  oDemes <- scan(paste0(mcmcpath[1], "/rdistoDemes.txt"), quiet = TRUE)
  oDemes <- matrix(oDemes, ncol = 3, byrow = TRUE)
  sizes <- as.matrix(oDemes[, 3])
  nPops <- nrow(oDemes)
  Demes <- seq(nPops)
  Sobs <- as.matrix(read.table(paste0(mcmcpath[1], "/rdistJtDobsJ.txt"), 
                               header = FALSE))
  Shat = matrix(nrow = nrow(Sobs), ncol = ncol(Sobs), 0)
  for (path in mcmcpath) {
    Shat <- Shat + as.matrix(read.table(paste0(path, "/rdistJtDhatJ.txt"), 
                                        header = FALSE))
  }
  Shat <- Shat/length(mcmcpath)
  colnames(Sobs) <- Demes
  rownames(Sobs) <- Demes
  colnames(Shat) <- Demes
  rownames(Shat) <- Demes
  
  if (!longlat) {
    long <- oDemes[, 2]
    lat <- oDemes[, 1]
  }
  else {
    long <- oDemes[, 1]
    lat <- oDemes[, 2]
  }
  x <- cbind(long, lat)
  Dist <- sp::spDists(x, x, longlat = TRUE)
  Dist <- Dist[upper.tri(Dist, diag = FALSE)]
  
  Sizes <- sizes %*% t(sizes)
  diag(Sizes) <- (sizes * (sizes - 1))/2
  
  df.between <- data.frame(Dist = Dist, Sobs = Sobs[upper.tri(Sobs)], 
                           Shat = Shat[upper.tri(Shat)], Sizes = Sizes[upper.tri(Sizes)], 
                           row.names = NULL)
  df.within <- data.frame(Sobs = diag(Sobs), Shat = diag(Shat), 
                          Sizes = diag(Sizes), row.names = NULL)
  
  # plot obs vs fitted genetic similarities between demes
  of_btw <- ggplot(df.between) + geom_point(aes(y = Sobs, x = Shat, size = Sizes), alpha = 0.6) +
    scale_size_continuous(range = c(1, 10), name = "# pairs") + 
    theme_classic() + geom_abline(intercept = 0) + ylab("genetic (observed) similarity between demes") + 
    xlab("fitted similarity between demes")
  
  ## plot semivariogram
  vg <- ggplot(df.between, aes(y = Sobs, x = Dist)) + geom_point(alpha = 0.6) +
    geom_smooth() +
    theme_classic() + geom_abline(intercept = 0) + theme(legend.position = 0) + 
    ylab("genetic similarity") + xlab("geographic distance")
  
  # plot obs vs fitted genetic similarities withing demes
  of_wth <- ggplot(df.within) + 
    geom_point(aes(y = Sobs, x = Shat, size = Sizes), alpha = 0.6) +
    scale_size_continuous(range = c(1, 10), name = "# pairs") + 
    theme_classic() + geom_abline(intercept = 0) + ylab("genetic (observed) similarity within demes") + 
    xlab("fitted similarity within demes")
  
  # save plots
  ggsave(filename = paste0(outpath, "/observed_vs_fitted-between.pdf"), plot = of_btw, width = 4, height = 4)
  ggsave(paste0(outpath, "/semivariogam.pdf"), plot = vg, width = 4, height = 4)
  ggsave(paste0(outpath, "/observed_vs_fitted-within.pdf"), plot = of_wth, width = 4, height = 4)
  
  return(list(OF_between = of_btw, OF_within = of_wth, semivg = vg))
}

plot_MAPS <- function (mcmcpath, outpath = "out_MAPS",  longlat = TRUE, width = 10, 
                       height = 6, add.pts = TRUE, add.graph = FALSE, add.countries = FALSE, 
                       plot.mean = TRUE, plot.difference = FALSE, correct.by.area = FALSE, col.palette = NULL,
                       set.range = FALSE, m.limits = NA, N.limits = NA, alpha = 0.1,
                       col.graph = "grey90",  col.pts = "black", fill.pts = "black",
                       const.size = FALSE, trans = NA, save.plot = TRUE) 
{
  out <- list()
  dir.create(file.path(outpath), showWarnings = FALSE)
  files <- c("/mcmcmtiles.txt", "/mcmcmrates.txt", "/mcmcxcoord.txt", 
             "/mcmcycoord.txt", "/mcmcqtiles.txt", "/mcmcqrates.txt", 
             "/mcmcpilogl.txt", "/rdistoDemes.txt", "/rdistJtDobsJ.txt", 
             "/rdistJtDhatJ.txt")
  mcmcpath <- check_files_at_path(files, mcmcpath)
  if (plot.mean) {
    plot.median = FALSE
  }
  else {
    plot.median = TRUE
  }
  params <- list(mcmcpath = mcmcpath, outpath = outpath, longlat = longlat, 
                 is.mrates = TRUE, plot.mean = plot.mean, plot.median = plot.median, correct.by.area = correct.by.area,
                 plot.sign = FALSE, plot.ci = '', width = width, height = height, add.countries = add.countries, 
                 add.graph = add.graph, palette = col.palette, 
                 col.graph = col.graph, add.pts = add.pts, col.pts = col.pts, fill.pts = fill.pts,
                 const.size = const.size, plot.difference = plot.difference, set.range = set.range, alpha = alpha, m.limits = m.limits, 
                 N.limits = N.limits, trans = trans, save.plot = save.plot)
  if (!is.na(params$trans)) {
    if (!(params$trans %in% c("identity", "log10"))) {
      stop("trans incorrectly specified, trans = identity or log10")
    }
  }
  if (params$set.range) {
    if (length(params$m.limits) < 2 | length(params$N.limits) < 
        2) {
      stop("set.range = TRUE but m.limits or N.limits not specified")
    }
  }
  dimns <- plotmaps:::read_dimns(params$mcmcpath[1], params$longlat)
  
  message("Computing migration statistics")
  m_summary_stats <- compute_summary_statistic(params, dimns)
  message("plotting migration surface")
  out$m <- plot_contour(params, m_summary_stats)
  
  message("plotting upper boundary migration surface")
  params$plot.ci = "upper"
  out$m_upper <- plot_contour(params, m_summary_stats)
  message("plotting lower boundary migration surface")
  params$plot.ci = "lower"
  out$m_lower <- plot_contour(params, m_summary_stats)
  message("plotting migration surface variation")
  params$plot.ci = "var"
  out$m_var <- plot_contour(params, m_summary_stats)
  
  if (!params$plot.difference) {
    message("plotting migration surface sign plot")
    params$plot.sign = TRUE
    out$m_sign <- plot_contour(params, m_summary_stats)
  }
  
  # reset params for population size plots
  params$plot.ci = ''
  params$is.mrates = FALSE
  params$plot.sign = FALSE
  
  message("Computing population size statistics")
  N_summary_stats <- compute_summary_statistic(params, dimns)
  message("plotting population-size surface")
  out$N <- plot_contour(params, N_summary_stats)
  
  message("plotting upper boundary migration surface")
  params$plot.ci = "upper"
  out$N_upper <- plot_contour(params, N_summary_stats)
  message("plotting lower boundary migration surface")
  params$plot.ci = "lower"
  out$N_lower <- plot_contour(params, N_summary_stats)
  message("plotting migration surface variation")
  params$plot.ci = "var"
  out$N_var <- plot_contour(params, N_summary_stats)
  
  if (!params$plot.difference) {
    message("plotting population-size surface sign plot")
    params$plot.sign = TRUE
    out$N_sign <- plot_contour(params, N_summary_stats)
  }
  
  message("plotting diagonostics of model fit and MCMC convergence")
  out$trace <- plot_trace(mcmcpath, outpath)
  out$fit <- plot_fit_data(mcmcpath, outpath, params$longlat)
  
  out$m_summary <- m_summary_stats
  out$N_summary <- N_summary_stats
  
  return(out)
}

select_best_loglik <- function(loglik_values, threshold = 0.05, min_group_size = 3) {
  # Sort log-likelihood values in descending order (best to worst)
  sorted_loglik <- sort(loglik_values, decreasing = TRUE)
  
  # Calculate differences between consecutive values
  diffs <- diff(sorted_loglik)
  
  # Find the largest gap that separates at least min_group_size values
  for (i in 1:(length(sorted_loglik) - min_group_size + 1)) {
    group <- sorted_loglik[i:(i + min_group_size - 1)]
    next_value <- sorted_loglik[i + min_group_size]
    
    if (is.na(next_value)) {
      # If we've reached the end, select all remaining values
      return(sorted_loglik[1:i + min_group_size - 1])
    }
    
    gap <- group[min_group_size] - next_value
    
    if (gap > threshold * abs(group[1])) {
      # If we find a significant gap, return the group above it
      return(sorted_loglik[1:(i + min_group_size - 1)])
    }
  }
  
  # If no significant gap is found, return all values
  return(sorted_loglik)
}

