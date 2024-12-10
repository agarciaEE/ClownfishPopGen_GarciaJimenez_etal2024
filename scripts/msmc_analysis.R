library(ggh4x)
library(ggtext)
library(ggplot2)
library(tidyverse)
library(ggExtra)
library(patchwork)
library(ggpattern)
library(ggpubr)
library(dplyr)
library(tidyr)
library(sf)
library(ggnewscale)
library(colorspace)

# working directory
setwd("~/Unil/Research/3_EcoPopGen/project_repository/")

### Custom Functions
source("./scripts/custom_functions.R")

#####
## Load samples info
#####
data <- read.csv("./raw_data/samples_info.csv")

colnames(data) <- c("sample_ID", "species", "abbrev", 
                    "country", "country_ID",
                    "location", "loc_ID", 
                    "lon", "lat", "category")

# host categories
guilds <- c("Generalist", "RM specialist", "EQ specialist", "SD specialist")
# species vector
species <- c("AKA", "AKY", "CLK", "LAT", "CRP", "MEL", "EPH", "PRD", "POL", "SAN")
# population vector
popLevels <- c("KNY", "MYT", "MDG", "MDV", "THA", "IDN", "PHI", "TAI",
               "PAU", "PNG", "SOI", "AUS", "NCD", "FIJ", "POL")
# population names
popnames <- setNames(c("Kenya", "Mayotte", "Madagascar", "Maldives", "Thailand", "Indonesia", "Philippines", "Taiwan",
                       "Palau", "Papua New Guinea", "Solomon Islands", "Australia", "New Caledonia", "Fiji", "French Polynesia"), 
                     popLevels)
# set factors
data$category <- factor(data$category, levels = guilds)
data$abbrev <- factor(data$abbrev, levels = species)
data$country_ID <- factor(data$country_ID, levels = popLevels)

# colors
popcols <- setNames(colorRampPalette(viridis::viridis(9))(n=length(popLevels)), popLevels)
guildcols <- setNames(c("black", "orange", "red", "skyblue"), guilds)

# generation time
gen = 5
# mutation rate
mu = 4E-08

#####
## sea level data
#####
SLdf <- read.table("./raw_data/sea_level_fluctuations/Sea-level_Miller_et_al.txt", header= TRUE)
SLdf <- SLdf[order(SLdf$Age, decreasing = TRUE),] # order by Age
SLdf$Age <- SLdf$Age * 1000000 # express in MY

# function to interpolate sea level relative to time
interp_func <- approxfun(SLdf$Age, SLdf$Sea.level, rule=1)

# subset Sea Level dataset to the time frame of msmc reconstructions
SLdf <- SLdf[SLdf$Age <= 2e6,]

# average change over periods of 25K years for each time point
time_period <- 5e4
SL_change25Ky <- sapply(SLdf$Age, function(age) {
  slc <- SLdf$Sea.level[SLdf$Age <= age & SLdf$Age >= age - time_period]
  mean(slc)
})

SL_change <- SLdf$Age[SL_change25Ky < -50]
SL_change <- c(SL_change, last(SL_change)-5e3)

# Initialize an empty data frame for storing periods
periods <- data.frame(start = numeric(0), end = numeric(0))
start_drop <- SL_change[3]
for (start_drop in SL_change){
  next_drops <- SL_change[SL_change < start_drop]
  obs_diff <- (next_drops - start_drop)
  exp_diff <- (max(diff(SLdf$Age)) * (1:length(SL_change[SL_change < start_drop])))
  d <- round(obs_diff / exp_diff, 2)
  if (any(d == 1)){
    end_drop = min(next_drops[d == 1])
    periods <- rbind(periods, data.frame(start = start_drop, end = end_drop))
    start_drop <-  SL_change[SL_change < end_drop][1]
  } 
}
periods <- periods[!duplicated(periods$end),]

# Plot the drop periods
ggplot() +
  geom_rect(data = periods, aes(xmin = start, xmax = end, ymin = min(SLdf$Sea.level), ymax = max(SLdf$Sea.level)), 
            size = 4, color = NA, fill = "darkred", alpha = 0.25) +
  geom_line(data = SLdf, aes(x = Age, y = Sea.level)) +
  scale_x_log10() +
  labs(title = "Drop Periods", x = "Time", y = "Sea Level") +
  theme_minimal()

# get the last sea-level drop
SL_lastdrop <- data.frame(left_time_boundary = last(periods$end),
                          right_time_boundary = last(periods$start))

## Read demographic reconstructions
###################################
# individual bootstraps Ne per population 
msmc10_data_filtered <- read.csv(file = "./analyses/msmc2/ind_bootstrap_msmc_data_filtered.csv", row.names = 1)
# multihet bootstraps Ne per population 
msmc50_data_filtered <- read.csv(file = "./analyses/msmc2/multihet_bootstrap_msmc_data_filtered.csv", row.names = 1)

# convert to factors
msmc10_data_filtered$species <- factor(msmc10_data_filtered$species, levels = species)
msmc10_data_filtered$pop <- factor(msmc10_data_filtered$pop, levels = popLevels)
msmc10_data_filtered$guild <- factor(msmc10_data_filtered$guild, levels = guilds)

# replace Inf in right time boundary time by left_time_boudnary 
msmc10_data_filtered$right_time_boundary[!is.finite(msmc10_data_filtered$right_time_boundary)] = msmc10_data_filtered$left_time_boundary[!is.finite(msmc10_data_filtered$right_time_boundary)] 

# get max time
maxTime <- max(c(msmc10_data_filtered$left_time_boundary, msmc10_data_filtered$right_time_boundary))

# get time period to weight Ne
msmc10_data_filtered$period <- msmc10_data_filtered$right_time_boundary - msmc10_data_filtered$left_time_boundary
msmc10_data_filtered$w <- msmc10_data_filtered$period/maxTime # relative to max time

# add sea level based on time points
msmc10_data_filtered$seaLevel <- interp_func(msmc10_data_filtered$time)

# axis labels and breaks for ploting
####################################
# x-axis 
exps <- (nchar(round(range(msmc10_data_filtered$time))))
xlims <- 10 ^ exps
xlabs <- sapply(exps[1]:exps[2], function(i) 10^i)
xticks <- unlist(lapply(2:length(xlabs), function(i) seq(xlabs[i-1], xlabs[i], length.out = 10)[1:9]))
xtickslabs <- sapply(xticks, function(i) ifelse(any(i == xlabs), format(i, scientific = TRUE, digits = 10), ""))
xlab <- bquote(paste("Years (g=", .(gen), ";", ~ mu, "=", .(format(mu, scientific = TRUE)), ")"))

# y-axis
ylims <- ceiling(range(msmc10_data_filtered$Ne) + c(-1, 0))
ylab <- expression("Effective population size" ~(x10^4))
labels <- sapply(levels(msmc10_data_filtered$species), function(i) paste0("*",i, "*"))

# visualize data
ggplot(msmc10_data_filtered, aes(x = time, y = Ne)) +
  geom_step(aes(group = interaction(species, pop, bootstrap), col = guild), linewidth = 1) +
  scale_x_log10(xlab,  breaks = xticks, labels = xtickslabs, expand = c(0,0)) +
  scale_y_continuous(ylab, expand = c(0,0), minor_breaks = seq(ylims[1], ylims[2], 0.25), guide = "axis_minor") + # this is added to the original code) 
  scale_color_manual("Ecological guild", values = guildcols) + 
  scale_fill_manual("Ecological guild", values = colorspace::lighten(guildcols, 0.8)) + 
  theme_linedraw() + 
  guides(colour = guide_legend(override.aes = list(linewidth = 2)), fill = "none") +
  theme(text = element_text(size = 18), 
        ggh4x.axis.ticks.length.minor = rel(-3/2),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        strip.text.x = element_text(size = 30),
        panel.spacing = unit(2, "lines"),
        legend.position = "none") +
  facet_grid(guild~., scales = "free")

###########################################################################################
#-----------------------------------------------------------------------------------------#
#                                  Analyses & Figures                                     #
#-----------------------------------------------------------------------------------------#
###########################################################################################

dir.create("Figures", showWarnings = FALSE)

##################################
#--------------------------------#
# All reconstructions per guild  #
#--------------------------------#
##################################

custom_theme <- theme_linedraw() + 
  theme(text = element_text(size = 20), 
        plot.margin = ggplot2::margin(0.25, 0.25, 0.25, 0.5, unit = "in"),
        plot.title = element_text(hjust = 0.5),
        ggh4x.axis.ticks.length.minor = rel(-3/2),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.text.x = element_text(size = 24, face = "italic"),
        strip.background = element_rect(fill = "black"),
        legend.text = element_text(size = 22), 
        legend.title = element_text(size = 28, face = "bold"),
        legend.key = element_rect(fill = NA, color = NA))

# Average the 10 individual-based bootstraps
summ_msmc10_bypop <- msmc10_data_filtered %>%
  group_by(guild, sci_name, species, pop) %>%
  group_modify(~ summary_by_intervals(.x))

summ_msmc10_bypop$time <- (summ_msmc10_bypop$left_time_boundary + summ_msmc10_bypop$right_time_boundary)/2

msmc_bootstrap_plotbyguild_List <- lapply(guilds, function(guild) {
  
  ggplot(msmc50_data_filtered[msmc50_data_filtered$guild == guild,], aes(time, Ne, col = pop)) +
    geom_step(aes(group = interaction(pop, bootstrap)), alpha = 0.25, linewidth = 0.25) +
    scale_color_manual(values = colorspace::lighten(popcols, 0.25), drop = FALSE, guide = "none") +
    new_scale_color() +
    geom_step(data = summ_msmc10_bypop[summ_msmc10_bypop$guild == guild,], aes(y = mean_Ne, col = pop), linewidth = 1, show.legend = TRUE) +
    scale_color_manual("", values = popcols, drop = FALSE, labels = popnames) +
    scale_x_log10(xlab, breaks = xticks, labels = xtickslabs) +
    scale_y_continuous(ylab, limits = ylims, minor_breaks = seq(ylims[1], ylims[2], 0.25), guide = "axis_minor") + 
    coord_cartesian(xlim = xlims, ylim = ylims) +
    guides(colour = guide_legend(override.aes = list(linewidth = 5, 
                                                     values = popcols, 
                                                     drop = FALSE, 
                                                     labels = popnames))) +
    facet_wrap(.~sci_name) + custom_theme + theme(strip.background = element_rect(fill = guildcols[guild]), 
                                                  strip.text = element_text(color = ifelse(guild == "Generalist", "white", "black")))
})
names(msmc_bootstrap_plotbyguild_List) <- guilds

p <- ggarrange(msmc_bootstrap_plotbyguild_List$Generalist,
               msmc_bootstrap_plotbyguild_List$`RM specialist`,
               msmc_bootstrap_plotbyguild_List$`EQ specialist`,
               msmc_bootstrap_plotbyguild_List$`SD specialist`, ncol = 1, heights = c(1.66, 1, 1, 1), 
               common.legend = TRUE, legend = "bottom", labels = letters[1:4], hjust = -0.5,
               font.label = list(size = 20, face = "bold", color ="black")) + bgcolor("white") + border(color = NA)
ggsave("./Figures/all_species_pop_msmc_bootstrap_mean.png", plot = p, width = 15, height = 20)

##################################
#--------------------------------#
# Ne vs Sea Level across guilds  #
#--------------------------------#
##################################

## (a) sea level fluctuations
l_ylims <- round(range(SLdf$Sea.level) + 10 * c(-1,0))

l <- ggplot() +
  geom_rect_pattern(data = cbind(start = periods$start[1],
                                 end = last(periods$end)), aes(xmin=start, xmax=end, ymin=l_ylims[1], ymax=l_ylims[2]), 
                    pattern = "stripe", fill = "grey90", pattern_fill ="white",
                    colour = "white", pattern_colour = "white", 
                    pattern_spacing = 0.10) + 
  geom_smooth(data = msmc10_data_filtered, 
              method = "loess",
              aes(x = time, y = seaLevel), 
              fullrange = TRUE, span = 0.075, na.rm = TRUE, se = F, col = "black") +
  #geom_line(data = seaLevel1, aes(x = Age, y = Sea.level )) + 
  scale_y_continuous("Sea level (m)", expand = c(0,0), limits = l_ylims) + 
  scale_x_log10("",  breaks = xticks, labels = xtickslabs, expand = c(0,0)) +
  coord_cartesian(xlim = c(1e+2, 1e+6), ylim = l_ylims) +
  theme_light() +
  theme(plot.margin =  margin(t = 0.5, r = 0.5, b = 0, l = 0.2, unit = "in"),
        text = element_text(size = 18, color = 'black'), 
        axis.text.x = element_text(size = 16, color = 'black'),
        axis.text.y = element_text(size = 16, color = 'black'),
        axis.ticks = element_line(color = "black"),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey85"), 
        panel.grid.minor.y = element_line(color = "grey85"),
        panel.border = element_rect(colour = "black", linewidth = 1))

# Add geological times
geological_periods <- data.frame(epoch = c("Holocene", "Pleistocene"),
                                 min_year = c(0, 12e3),
                                 max_year = c(12e3, 2.6e6),
                                 y = l_ylims[2]*0.35)

l <- l +
  geom_rect(data = geological_periods, 
            aes(xmin = min_year, xmax = max_year, ymin = 5, ymax = l_ylims[2], fill = epoch), col = "black",
            inherit.aes = FALSE) + 
  geom_text(data = geological_periods, 
            aes(x = (min_year + max_year) / 2, y = y, label = epoch), 
            inherit.aes = FALSE, 
            size = 5, vjust = 0, hjust = 2.25) +
  scale_fill_manual(values = c("#FFDDC1", "#E2FFC1")) + theme(legend.position = "none")


## (b) Ne average reconstructions by species per guild
# Average Ne across species per guild
summ_res_guild <- na.exclude(msmc10_data_filtered %>%
                               group_by(guild) %>%
                               group_modify(~ summary_by_intervals(.x)))

# Average Ne across populations per species
summ_res_sp <- na.exclude(msmc10_data_filtered %>%
                            group_by(guild, species) %>%
                            group_modify(~ summary_by_intervals(.x)))

# add time
summ_res_guild$time <- (summ_res_guild$left_time_boundary + summ_res_guild$right_time_boundary)/2
summ_res_sp$time <- (summ_res_sp$left_time_boundary + summ_res_sp$right_time_boundary)/2

s <- ggplot() +
  geom_rect_pattern(data = cbind(start = periods$start[1],
                                 end = last(periods$end)), aes(xmin=start, xmax=end, ymin=ylims[1], ymax=ylims[2]),
                    fill = "grey90", pattern_fill ="white",
                    colour = "white", pattern_colour = "white") +
  geom_rect(data = summ_res_guild, aes(xmin = left_time_boundary, xmax = right_time_boundary, 
                                       ymin = min_Ne, ymax = max_Ne, fill = guild)) +
  geom_step(data = summ_res_sp, aes(x = time, y = mean_Ne, group = species, col = guild), linewidth = 1) +
  scale_x_log10(xlab,  breaks = xticks, labels = xtickslabs, expand = c(0,0)) +
  coord_cartesian(xlim = c(1e+2, 1e+6), ylim = c(0, 8)) +
  scale_y_continuous(ylab, limits = c(ylims[1], ylims[2]), expand = c(0,0), minor_breaks = seq(ylims[1], ylims[2], 0.25), guide = "axis_minor") + # this is added to the original code) 
  scale_color_manual("Ecological guild", values = c("black", "orange", "red", "skyblue")) + 
  scale_fill_manual("Ecological guild", values = colorspace::lighten(c("black", "orange", "red", "skyblue"), 0.8)) + 
  theme_linedraw() + 
  guides(colour = guide_legend(override.aes = list(linewidth = 2)), fill = "none") +
  theme(text = element_text(size = 18), 
        ggh4x.axis.ticks.length.minor = rel(-3/2),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        strip.text.x = element_text(size = 30),
        panel.spacing = unit(2, "lines"),
        legend.position = "none") +
  facet_grid(guild~., scales = "free")

## (c) Ne differences between raise and drop periods per guild

# differentiate sea-level rise (Holocene) from sea-level fluctuation (Pleistocene)
msmc10_data_filtered$type[msmc10_data_filtered$time < as.numeric(SL_lastdrop$left_time_boundary)] = "raise"
msmc10_data_filtered$type[msmc10_data_filtered$time >= as.numeric(SL_lastdrop$left_time_boundary)] = "drop"
msmc10_data_filtered$type <- factor(msmc10_data_filtered$type, levels = c("raise", "drop"))

# compute harmonic mean grouped by type, guild, sp_pop, and bootstrap
WHM_Ne_drop <- msmc10_data_filtered %>%
  dplyr::group_by(type, guild, species, pop, bootstrap) %>%
  dplyr::summarize(WHM_Ne = harmonic_mean(Ne, w), .groups = 'drop')

WHM_Ne_drop$behvaior = factor(ifelse(WHM_Ne_drop$guild == "Generalist", "Generalist", "Specialist"),
                              levels = c("Generalist", "Specialist"))

# Test raise vs drop differences
# Function to choose between Wilcoxon and t-test
perform_test <- function(raise, drop) {
  # Check for normality using Shapiro-Wilk test
  raise_normal <- shapiro.test(raise)$p.value > 0.05
  drop_normal <- shapiro.test(drop)$p.value > 0.05
  
  # Choose the appropriate test based on normality
  if (raise_normal & drop_normal) {
    # Perform paired t-test if both distributions are normal
    test_result <- t.test(raise, drop, paired = TRUE)
    test_type <- "t-test"
  } else {
    # Perform Wilcoxon test if one or both are not normal
    test_result <- wilcox.test(raise, drop, paired = TRUE)
    test_type <- "Wilcoxon"
  }
  
  # Return a data frame with the statistic, p-value, and test type
  data.frame(
    Ne_raise = mean(raise, na.rm = TRUE),
    Ne_drop = mean(drop, na.rm = TRUE),
    statistic = test_result$statistic,
    p_value = test_result$p.value,
    test_type = test_type
  )
}

# Apply the test for each guild and combine results
raise_drop_test_results <- lapply(guilds, function(g) {
  # Reshape the data from long to wide format
  df <- WHM_Ne_drop %>%
    spread(type, WHM_Ne)
  raise <- df$raise[df$guild == g]
  drop <- df$drop[df$guild == g]
  # Perform the appropriate test
  result <- perform_test(raise, drop)
  # Add guild information to the result
  result$guild <- g
  result
})

raise_drop_test_results <- do.call(rbind, raise_drop_test_results)

# % Difference between recent and past Ne per group
(raise_drop_test_results$Ne_raise - raise_drop_test_results$Ne_drop)/raise_drop_test_results$Ne_drop * 100 

# custom statistic label to add to the plot
raise_drop_test_results$label <- paste0(ifelse(raise_drop_test_results$test_type == "Wilcoxon", "W", "t"), " = ", format(raise_drop_test_results$statistic, digits = 3),
                                        ifelse(raise_drop_test_results$p_value < 0.001, paste0("; p < 0.001"),
                                               ifelse(raise_drop_test_results$p_value < 0.01, paste0("; p < 0.01"),
                                                      ifelse(raise_drop_test_results$p_value < 0.05, paste0("; p < 0.05"),
                                                             paste0("; p = ", round(raise_drop_test_results$p_value, 3))))))

# ploting parameters
dodge <- position_dodge(width = 0.8)
symnum.args <- list(cutpoints = c(0,0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "n.s."))

N <- ggplot(WHM_Ne_drop, aes(x = guild, y = WHM_Ne,  ill=interaction(guild,type))) 
for (i in 1:4) {
  N <- N + geom_rect_pattern(xmin=i, xmax=i + 0.4, ymin=0, ymax=4.5, 
                             pattern = "stripe", fill = "grey90", pattern_fill ="white",
                             colour = "white", pattern_colour = "white", 
                             pattern_spacing = 0.10)
}
N <- N + geom_boxplot(lwd = 2, position = dodge, aes(col = guild, fill = guild), outlier.size = 3) +
  geom_boxplot(fatten = 5, position = dodge, col = "white", aes(fill = guild)) +
  scale_x_discrete(name = "", labels = gsub("_", " ", levels(WHM_Ne_drop$guild))) +
  scale_y_continuous(name = expression("Historical" ~ "N"[e] ~(x10^4)), expand = c(0,0)) + 
  scale_fill_manual(values = guildcols) +
  scale_color_manual(values = guildcols) +
  coord_cartesian(ylim = c(0, 5)) +
  theme_light() +
  theme(axis.text.y = element_text(size = 16, color = 'black'),
        axis.ticks.x = element_line(color = "black", linewidth = 2),
        axis.text.x = element_text(size = 16, color = 'black'),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey85"), 
        panel.grid.minor.y = element_line(color = "grey85"),
        panel.border = element_rect(colour = "black", linewidth = 1),
        legend.position = "none", 
        axis.title.y = element_text(vjust = -0.25, size = 18)) +
  geom_text(data = raise_drop_test_results, aes(x = guild, y = 4.75, 
                                                label = label), 
            size = 4, inherit.aes = FALSE)

fig <- l/s/N + 
  plot_layout(nrow = 3, heights = c(0.4,2,0.5)) + 
  plot_annotation(tag_levels = "a", tag_suffix = "  ")  + 
  theme(plot.tag = element_text(size = 14))

ggsave(filename = "./Figures/msmc_seaLevel_guilds_with_bootstraps.png", 
       plot = fig, width = 10, height = 17)

##################################
#--------------------------------#
# Plot map with msmc per country #
#--------------------------------#
##################################

# get populations lon lat coordinates
pop_coordinates <- data.frame(
  pop = popLevels,
  x = sapply(popLevels, function(i) mean(data$lon[data$country_ID == i], na.rm = TRUE)),
  y = sapply(popLevels, function(i) mean(data$lat[data$country_ID == i], na.rm = TRUE))
)
pop_coordinates <- na.exclude(pop_coordinates)
pop_coordinates$x[pop_coordinates$x < 0] <- pop_coordinates$x[pop_coordinates$x < 0] + 360 # convert longitudes from -180,180 to 0,360

# add country names
pop_coordinates$country <- popnames[pop_coordinates$pop]

# get world map marine regions
regfile <- "./raw_data/MEOW_ECOS/"
regions <- sf::st_read(regfile)
regions <- sf::st_transform(regions, "+init=epsg:4326") # transform to WGS84
regions <- sf::st_shift_longitude(regions) # convert longitudes from -180,180 to 0,360
# subset regions of interest
regions <- regions[regions$REALM %in% c("Central Indo-Pacific", "Western Indo-Pacific", "Eastern Indo-Pacific",
                                        "Temperate Australasia") | regions$ECOREGION %in% c("East China Sea", "Central Kuroshio Current"),]
regions <- regions[!regions$PROVINCE %in% "Easter Island",]
regions <- regions[!regions$PROVINCE %in% unique(regions$PROVINCE[regions$REALM == "Temperate Australasia" & !regions$PROVINCE %in% c("West Central Australian Shelf", "East Central Australian Shelf")]),]

# Get world map data
world_map <- sf::st_as_sf(maps::map("world", wrap = c(0,360), plot = FALSE, fill = TRUE)) 

# Average Ne bootstraps per populations per species
summ_res_pop <- na.exclude(msmc10_data_filtered %>%
                             group_by(guild, species, pop) %>%
                             group_modify(~ summary_by_intervals(.x)))

# add country name
summ_res_pop$country <- sapply(summ_res_pop$pop, function(pop) pop_coordinates$country[pop_coordinates$pop == pop])
# add time
summ_res_pop$time <- (summ_res_pop$left_time_boundary + summ_res_pop$right_time_boundary)/2

# create a grid to divide ploting areas
grid <- expand.grid(seq(10, 210, 10), seq(-50, 50, 10))

buffer = 20 # separation between each location-specific msmsc plot
for (i in 1:nrow(pop_coordinates)) {
  x <- pop_coordinates$x[i]
  y <- pop_coordinates$y[i]
  # remove grid cells overlaying lownfish distributions
  d <- apply(grid, 1, function(i) sqrt((i[1]-x)^2 + (i[2]-y)^2))
  idx <- which(d < buffer)
  if (length(idx)>0){
    grid <- grid[-idx,]
  }
}

# create map baseplot
wmap <- ggplot(world_map) +
  geom_sf(data = regions, fill = "lightskyblue", col = NA, alpha = 0.75) +
  geom_sf(fill = "grey90", size = 1) +
  coord_sf(xlim = c(1, 250), ylim = c(-60,60), expand = FALSE) 

# order popLevels for visualization purposes
for (pop in popLevels[c(1:6, 9,7,8,10:15)]) {
  
  x <- pop_coordinates$x[pop_coordinates$pop == pop]
  y <- pop_coordinates$y[pop_coordinates$pop == pop]
  
  # get closest grid cell to location
  d <- apply(grid, 1, function(i) sqrt((i[1]-x)^2 + (i[2]-y)^2))
  idx <- which.min(d)
  
  x0 <- grid$Var1[idx]
  y0 <- grid$Var2[idx]
  
  # remove grid cells within a distance from selected cell
  d <- apply(grid, 1, function(i) sqrt((i[1]-x0)^2 + (i[2]-y0)^2))
  idx <- which(d < (25))
  grid <- grid[-idx,]
  
  # subset reconstruction for selected location
  df <- summ_res_pop[summ_res_pop$pop == pop,]
  
  # plot population msmc 
  p <- ggplot(df, aes(x = time, y = mean_Ne)) +
    geom_step(aes(group = species, col = guild), linewidth = 1) +
    scale_x_log10("",  breaks = xticks, limits = xlims,
                  labels = xtickslabs, expand = c(0,0.5)) +
    scale_y_continuous("", expand = c(0.2,0), limits = ylims, breaks = seq(ylims[1], ylims[2], by = 2)) +
    coord_cartesian(xlim = c(1e+2, 1e+6), ylim = ylims) +
    scale_color_manual("Host category", values = guildcols) +
    theme_linedraw() +
    theme(
      plot.background = element_rect(fill='transparent', color=NA),
      plot.title = element_text(hjust = 0.5),
      axis.ticks.length=unit(-0.1, "cm"),
      axis.text.x = element_text(size = 8, vjust = 8),
      axis.text.y = element_text(size = 8, margin = ggplot2::margin(r = -10)),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      strip.text.x = element_text(size = 12),
      legend.position = "none")  +
    facet_wrap(.~country, scales = "free")
  
  # get distance from location to the grid cell
  xy0dist <- sqrt((x-x0)^2 + (y-y0)^2)
  
  # if distance is higher than 15, draw a line connector
  if (xy0dist > 15) {
    wmap <- wmap + geom_segment(x = x, y = y, xend = x0, yend = y0, colour = "grey50") +
      geom_point(x = x, y = y, size = 3, col = "grey20", shape = 21, fill = "yellow")
  }
  # add plot to the map at the selected grid cell
  wmap <- wmap +  annotation_custom(grob = ggplotGrob(p),
                                    xmin = x0 - 15.5, xmax = x0 + 15.5,
                                    ymin = y0 - 13.5, ymax = y0 + 13.5) + 
    theme(legend.position = "none")
}

ggsave("./Figures/map_msms_by_pop.png", plot = wmap, height = 10, width = 20)

##################################
#--------------------------------#
#    Regression models of Ne     #
#--------------------------------#
##################################

##########
# Fit GLMM
##########

glmm_model_guild_SL <- glmer(Ne ~ scale(seaLevel) * guild * scale(time) + pop + (1 | species/pop/bootstrap), 
                             family = Gamma(link = "log"),
                             data = msmc10_data_filtered,
                             control = glmerControl(optimizer = "bobyqa", 
                                                    optCtrl = list(maxfun = 1e5)))
summary(glmm_model_guild_SL)

# plot
plot_model(glmm_model_guild_SL, type = "pred", terms = c("time", "seaLevel", "guild"), line.size = 2, pred.type = "fe", colors = c("#EF8A62", "grey30", "#67A9CF")) + 
  coord_cartesian(ylim = c(0, 8))
ggsave("./Figures/msmc_res_glmm_guilds_model.png")

# table
tidy(glmm_model_guild_SL, conf.int = TRUE) %>%
  filter(effect == "fixed") %>%
  mutate(
    estimate = sprintf("%.3f", estimate),
    conf.low = sprintf("%.3f", conf.low),
    conf.high = sprintf("%.3f", conf.high),
    significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""),
    p.value = format.pval(p.value, digits = 3),
  ) %>%
  dplyr::select(term, estimate, conf.low, conf.high, p.value, significance) %>%
  gt() %>% 
  gt::gtsave("./Figures/msmc_res_glmm_guilds_regression_table.html")
pdf("./Figures/msmc_res_glmm_guild_residuals.pdf", width = 10, height = 7)
eval_model(glmm_model_guild_SL)
dev.off()

################################################
#----------------------------------------------#
# Functional Data Analysis (FDA) across guild  #
#----------------------------------------------#
################################################

# Interpolation of time points to a common grid
################################################
# create a common time grid
min_time <- 0
max_time <- round(max(msmc10_data_filtered$time)/1E+6)*1E+6 # max time rounded to the million
common_time_grid <- sort(unique(as.numeric(sapply(c(1:6), function(i) seq(0, 10, by = 1) * 10 ^ i)))) # even time grid in log scale
common_time_grid <- common_time_grid[common_time_grid <= max_time]

# Interpolate population sizes for each species/population/bootstrap combination
interpolated_data <- msmc10_data_filtered %>%
  group_by(guild, species, pop, bootstrap) %>%
  do({
    msmc_res_Ne <- approx(.$time, .$Ne, xout = common_time_grid, rule = 2)$y 
    data.frame(time = common_time_grid, Ne = msmc_res_Ne)
  }) %>%
  ungroup()

interpolated_data$id = paste0(interpolated_data$species, "_", interpolated_data$pop, ".", interpolated_data$bootstrap)
# create binary category differentiating generalists and specialists
interpolated_data$behavior <- factor(ifelse(interpolated_data$guild == "Generalist", "generalist", "specialist"), levels = c("generalist", "specialist"))
# interpolate sea level relative to time
interpolated_data$seaLevel <- interp_func(interpolated_data$time)

# visualize interpolated data
ggplot(na.exclude(interpolated_data), aes(x = time, y = Ne)) +
  geom_step(aes(group = interaction(species, pop, bootstrap), col = guild), linewidth = 1) +
  scale_x_log10(xlab,  breaks = xticks, labels = xtickslabs, limits = xlims) +
  scale_y_continuous(ylab, limits = c(0, ceiling(max(interpolated_data$Ne, na.rm = TRUE))), expand = c(0,0), 
                     minor_breaks = seq(0, ceiling(max(interpolated_data$Ne, na.rm = TRUE)), 0.25), guide = "axis_minor") + # this is added to the original code) 
  scale_color_manual("Host category", values = guildcols) + 
  scale_fill_manual("Host category", values = colorspace::lighten(guildcols, 0.8)) + 
  coord_cartesian(xlim = c(1e+2, 1e+6), ylim = ylims) +
  theme_linedraw() + 
  guides(colour = guide_legend(override.aes = list(linewidth = 2)), fill = "none") +
  theme(text = element_text(size = 18), 
        ggh4x.axis.ticks.length.minor = rel(-3/2),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        strip.text.x = element_text(size = 30),
        panel.spacing = unit(2, "lines"),
        legend.position = "none") +
  facet_grid(guild~., scales = "free")

###########################
# Functional ANOVA (FANOVA)
###########################
# create a matrix with time points as rows and species-population-bootstrap reconstructions as columns
interpolated_data_sprd <- as.data.frame(spread(interpolated_data[,c("time", "Ne", "id")], id, Ne))
rownames(interpolated_data_sprd) <- interpolated_data_sprd$time
interpolated_data_sprd <- interpolated_data_sprd[,-1] # remove time column

# get host category per column
guild_labels <- sapply(names(interpolated_data_sprd), function(i) interpolated_data$guild[interpolated_data$species == substr(i, 1, 3)][1])

plotFANOVA(interpolated_data_sprd, group.label = guild_labels, means = TRUE) + 
  scale_color_manual("Host category", values = guildcols) +
  scale_x_continuous(xlab,  breaks = which(log10(common_time_grid) %% 1 == 0), 
                     labels = c(format(common_time_grid[which(log10(common_time_grid) %% 1 == 0)], scientific = TRUE)))  +
  scale_y_continuous(expression("Effective population size" ~(x10^4))) +
  theme_linedraw() + 
  guides(colour = guide_legend(override.aes = list(linewidth = 2)), fill = "none") +
  theme(text = element_text(size = 24), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.spacing = unit(2, "lines"))
ggsave("./Figures/interpolated_FANOVA_guilds.png")

# run FANOVA tests
fanova_result_guild <- fanova.tests(na.exclude(interpolated_data_sprd), group.label = guild_labels)
# Print FANOVA results
cat(paste0(capture.output(print(fanova_result_guild)), "\n"), file = "./Figures/interpolated_FANOVA_behavior_tests_guilds.txt")
print(fanova_result_guild)

# get specialist vs generalist vector
behavior_labels <- sapply(names(interpolated_data_sprd), function(i) interpolated_data$behavior[interpolated_data$species == substr(i, 1, 3)][1])

plotFANOVA(interpolated_data_sprd, group.label = behavior_labels, means = TRUE) +
  scale_color_manual("Mutualistic\n behavior", values = c("#E8C241", "#3171BD")) +
  scale_x_continuous(xlab,  breaks = which(log10(common_time_grid) %% 1 == 0), 
                     labels = c(format(common_time_grid[which(log10(common_time_grid) %% 1 == 0)], scientific = TRUE)))  +
  scale_y_continuous(expression("Effective population size" ~(x10^4))) +
  theme_linedraw() + 
  guides(colour = guide_legend(override.aes = list(linewidth = 2)), fill = "none") +
  theme(text = element_text(size = 24), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.spacing = unit(2, "lines"))
ggsave("./Figures/interpolated_FANOVA_behavior.png")

# run FANOVA tests
fanova_result_behavior <- fanova.tests(interpolated_data_sprd, group.label = behavior_labels)
# Print FANOVA results
cat(paste0(capture.output(print(fanova_result_behavior)), "\n"), file = "./Figures/interpolated_FANOVA_behavior_tests_behavior.txt")
print(fanova_result_behavior)

