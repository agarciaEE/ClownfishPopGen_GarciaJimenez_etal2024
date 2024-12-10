library(tidyverse)
library(ggplot2)
library(grid)
library(gridExtra)
library(data.table)
library(plyr)
library(dplyr)
library(ggpubr)
library(ggrepel)
library(adegenet)
library(reshape2) 
library(geosphere)
library(ggnewscale)
library(ggthemes)
library(colorspace)
library(broom)
library(sf)
library(nlme)
library(sjPlot)

# working directory
setwd("~/Unil/Research/3_EcoPopGen/project_repository/")

### Custom Functions
source("./scripts/custom_functions.R")

#####
## Load samples info
#####
data <- read.csv("./raw_data/samples_info.csv")

# shift longitudes to 0-360
data$lon[!is.na(data$lon) & data$lon < 0] = data$lon[!is.na(data$lon) & data$lon < 0]  + 360

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


# rename host categories (guilds) according to new nomenclature
data$category <- factor(data$category, levels = guilds)

# colors
popcols <- setNames(colorRampPalette(viridis::viridis(9))(n=length(popLevels)), popLevels)
guildcols <- setNames(c("black", "orange", "red", "skyblue"), guilds)

pop_count <- data %>%
  group_by(abbrev, country_ID) %>%
  dplyr::summarise(lon = mean(lon, na.rm = T),
                   lat = mean(lat, na.rm = T),
                   n = length(country))

# host category comparisons
mycomp <-  combn(levels(data$category), 2)
mycomp <- lapply(1:ncol(mycomp), function(i) mycomp[,i])
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "n.s."))

###########################################
## ------------------------------------- ##
##              Fst, pi, Dxy             ##
## ------------------------------------- ##
###########################################
Avg_popFst_List <- list()
Avg_popDxy_List <- list()
Avg_popPi_List <- list()

for (s in species){
  
  path <- paste0("./analyses/popgen/", s, "/")

  pops <- unique(data$country_ID[data$abbrev == s])
  pops <- popLevels[popLevels %in% pops]
  
  # check all files are present
  if (s %in% c("EPH", "LAT")) {
    missing.FST = FALSE
  } else {
    # check missing files
    comb_pops <- apply(combn(pops, 2), 2, function(i) paste(i, collapse = "-"))
    
    # read files
    pop_fst_files <- list.files(path, pattern = "[A-Z].weir.fst")
    pop_fst_comps <- gsub("MVD", "MDV", gsub(paste0(s, "_"), "", 
                                             sapply(pop_fst_files, function(i) str_split(i, "\\.")[[1]][1])))
    
    missing.FST = FALSE
    if (!all(comb_pops %in% pop_fst_comps)) {
      missing_pop_comps <- sapply(comb_pops[!comb_pops %in% pop_fst_comps], function(i) paste0(substr(i, 5,7), "-", substr(i, 1,3)))
      if (!all(missing_pop_comps %in% pop_fst_comps)){
        cat("\nMissing ", s, " fst pop comparisons:\n", 
            comb_pops[!comb_pops %in% pop_fst_comps])
        missing.FST = TRUE
      }
    }
  }
  
  pop_pi_files <- list.files(pattern = "[A-Z].sites.pi")
  pops_pi <- gsub("MVD", "MDV", gsub(paste0(s, "_"), "", 
                                     sapply(pop_pi_files, function(i) str_split(i, "\\.")[[1]][1])))
  
  missing.PI = FALSE
  if (!all(pops %in% pops_pi)) {
    cat("\nMissing ", s, " pi pop files:\n", 
        pops[!pops %in% pops_pi])
    missing.PI = TRUE
  }
  
  if (!missing.PI | missing.FST){
    
    if (!s %in% c("EPH", "LAT")) {
      
      # Fst
      ######
      # read .weir.fst files
      popFst <- lapply(paste0(path, pop_fst_files), fread, header=T)
      names(popFst) <- pop_fst_comps
      
      Avg_popFst_df <- data.frame(t(combn(pops, 2)))
      Avg_popFst_df$value <- apply(Avg_popFst_df, 1, function(row) {
        comp <- c(paste(row, collapse = "-"), paste(row[2:1], collapse = "-"))
        comp <- comp[comp %in% pop_fst_comps]
        mean(popFst[[comp]][[3]], na.rm = TRUE)
      })
      colnames(Avg_popFst_df) <- c("pop1", "pop2", "value")
      
      Avg_popFst_df$pop1 <- factor(Avg_popFst_df$pop1, levels = pops)
      Avg_popFst_df$pop2 <- factor(Avg_popFst_df$pop2, levels = pops)
      
      Avg_popFst_df$value[Avg_popFst_df$value < 0 ] = 0
      Avg_popFst_List[[s]] <- Avg_popFst_df

      # Dxy
      ######
      # read files
      popDxy <- as.data.frame(fread(paste0(path, s, "_dxy.csv"), header=T))
      colnames(popDxy) <- gsub("MVD", "MDV", colnames(popDxy))      # correct Maldives population ID
      popDxy <- popDxy[, c(1:5, grep("dxy", colnames(popDxy)))]
      
      # get comparisons
      pop_Dxy_comps <- gsub("dxy_", "", grep("dxy", colnames(popDxy), value = TRUE))
      colnames(popDxy) <- gsub("dxy_", "", colnames(popDxy))
      
      Avg_popDxy_df <- data.frame(t(combn(pops, 2)))
      Avg_popDxy_df$value <- apply(Avg_popDxy_df, 1, function(row) {
        comp <- c(paste(row, collapse = "_"), paste(row[2:1], collapse = "_"))
        comp <- comp[comp %in% pop_Dxy_comps]
        mean(popDxy[,comp], na.rm = TRUE)
      })
      colnames(Avg_popDxy_df) <- c("pop1", "pop2", "value")
      
      Avg_popDxy_df$pop1 <- factor(Avg_popDxy_df$pop1, levels = pops)
      Avg_popDxy_df$pop2 <- factor(Avg_popDxy_df$pop2, levels = pops)
      
      Avg_popDxy_df$value[Avg_popDxy_df$value < 0 ] = 0
      Avg_popDxy_List[[s]] <- Avg_popDxy_df
      
    }
    
    # Pi
    ######
    pop_pi_files <- list.files(path, pattern = "[A-Z].sites.pi")
    popPi <- lapply(paste0(path, pop_pi_files), fread, header=T)
    names(popPi) <- pops_pi
    
    # remove species file
    popPi <- popPi[names(popPi) != s]
    
    # compute pi per 10 Kb windows
    site.pi2window <- function(pi_data, window_size = 1e4){
      
      # Summarize to get the sum of PI for each window and divide by window size
      windowed_pi <- pi_data %>%
        mutate(WINDOW = floor(POS / window_size) * window_size) %>%
        dplyr::group_by(CHROM, WINDOW) %>%
        dplyr::summarise(
          n_poly_sites = sum(PI > 0),
          n_sites = n(),
          avg_pi_site = sum(PI, na.rm = TRUE) / window_size,
          avg_pi_window = n_poly_sites / window_size # S / L
        )
      
      return(windowed_pi)
    }
    
    popPi_window <- sapply(popPi, site.pi2window, simplify = FALSE)
    
    Avg_popPi_df <- data.frame(pop = pops,
                               value = sapply(pops, function(pop) mean(popPi_window[[pop]]$avg_pi_site, na.rm = TRUE))
    )

    Avg_popPi_List[[s]] <- Avg_popPi_df
  
  }
}

#########
#### plot
#########
##############
# Heatmaps (fst, dxy, pi)
##############
pi_lim <- ceiling(max(unlist(lapply(Avg_popPi_List, function(i) i[,2]))) * 1e2) / 1e2

p_List <- list()
for (s in  species[!species %in% c("EPH", "LAT")]){
  
  pops <- popLevels[popLevels %in% data$country_ID[data$abbrev == s]]
  spname <- sub("Amphiprion_", "A. ", data$species[data$abbrev == s][1])
  
  p <- ggplot() +
    # Upper diagonal for Dxy
    geom_tile(data = Avg_popDxy_List[[s]],
              aes(pop1, pop2, fill = value), color = "white", linewidth = 0.5) +
    geom_text(data = Avg_popDxy_List[[s]],  size = 8,
              aes(pop1, pop2, label = ifelse(is.na(value), "NA", round(value, 2)), 
                  color = ifelse(is.na(value), "white", "black"))) +  
    scale_color_manual(values = c("white" = "white", "black" = "black")) +
    scale_fill_gradientn(colors = c("white", "#ED8F5F"), 
                         limits = c(0, 0.5),
                         name = expression(italic(d[xy])), 
                         guide = "none", 
                         na.value = "black") + 
    # Lower diagonal for FST
    new_scale_fill() +
    geom_tile(data = Avg_popFst_List[[s]],
              aes(pop2, pop1, fill = value), color = "white", linewidth = 0.5) +
    geom_text(data = Avg_popFst_List[[s]], size = 8,
              aes(pop2, pop1, label = ifelse(is.na(value), "NA", round(value, 2)), 
                  color = ifelse(is.na(value), "white", "black"))) +  
    scale_color_manual(values = c("white" = "white", "black" = "black")) +
    scale_fill_gradientn(colors = c("white", "violetred4"), 
                         limits = c(0, 0.5),
                         name = expression(F[ST]), 
                         guide = "none", 
                         na.value = "black") +
    # Diagonal for PI
    new_scale_fill() +
    geom_tile(data = Avg_popPi_List[[s]] %>%
                filter(pop %in% pops),
              aes(pop, pop, fill = value), color = "black", linewidth = 1.5) +
    geom_text(data = Avg_popPi_List[[s]] %>%
                filter(pop %in% pops),  size = 8,
              aes(pop, pop, label = round(value, 3)), color = "black") +  # Add FST values inside tiles
    scale_fill_gradientn(colors = c("white", "steelblue4"),
                         limits = c(0, pi_lim),
                         name = expression(pi),
                         guide = "none", 
                         na.value = "black") +
    scale_y_discrete(position = "right", labels = popnames[pops]) +
    scale_x_discrete(labels = popnames[pops]) +
    labs(title = spname, x = "", y = "") +
    theme_bw() +
    guides(color = "none") + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 14),
          axis.text.y = element_text(hjust = 1, size = 14),
          plot.title = element_text(size = 20, face = "italic", vjust = 1, hjust = 0.5),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_blank(),
          panel.background = element_rect(fill = "transparent"),  
          plot.margin = ggplot2::margin(0.5,0.5,0.5,0.5, unit = "cm"),
          plot.background = element_rect(fill = "white"),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16, face = "italic"))

  p_List[[s]] <- p
}

# Extract legends as grobs
extract_legend <- function(color, limits = c(0, 0.5), 
                           name = expression(italic(d[xy])), size = 6, 
                           width = 15, height = 15) {
  
  plot <- ggplot() +
    geom_tile(data = Avg_popDxy_List[[s]],
              aes(pop1, pop2, fill = value), color = "white") +
    scale_fill_gradientn(colors = c("white", color), 
                         limits = limits,
                         name = name, 
                         guide = guide_colorbar(position = "top", title.position = "right", 
                                                theme = theme(legend.key.width  = unit(width, "cm"),
                                                              legend.key.height  = unit(height, "cm"),
                                                              legend.title = element_text(size = size),
                                                              legend.text = element_text(size = size * 0.6))), 
                         na.value = "black")
  ggplotGrob(plot)$grobs[[which(sapply(ggplotGrob(plot)$grobs, function(x) x$name) == "guide-box")]]
}

legend_dxy <- extract_legend("#ED8F5F", limits = c(0, 0.5), name = expression(italic(d[xy])), 
                             size = 20,  width = 15, height = 2)
legend_fst <- extract_legend("violetred4", limits = c(0, 0.5), name = expression(italic(F[ST])), 
                             size = 20, width = 15, height = 2)
legend_pi <- extract_legend("steelblue4", limits = c(0, pi_lim), name = expression(pi), 
                            size = 20, width = 15, height = 2)

# Arrange the legends side-by-side
grid_legends <- ggarrange(legend_dxy, legend_fst, legend_pi, ncol = 1) + bgcolor("white") + border("white")

pL <- ggarrange(p_List[["CLK"]], p_List[["CRP"]],
          p_List[["AKY"]], p_List[["MEL"]],
          p_List[["AKA"]], p_List[["PRD"]],
          p_List[["POL"]], p_List[["SAN"]],
          grid_legends, ncol = 3, nrow = 3, font.label = list(size = 26, face = "bold", color ="black"),
          labels = letters[1:8])

ggsave(paste0("./Figures/pi_dxy_Fst.png"), plot = pL, width = 30, height = 30, limitsize = FALSE)

##############
# Adjust by geo distance and population split time
##############
geodist_splittime_df <- read.csv("./analyses/msmc2/pop_geodist_divtime.csv")

##
# Fst
##
Avg_popFst_data <- na.exclude(plyr::ldply(Avg_popFst_List, .id = "species"))
names(Avg_popFst_data)[4] = "fst"

for (i in 1:nrow(Avg_popFst_data)) {
  sp <- Avg_popFst_data$species[i]
  pop1 <- Avg_popFst_data$pop1[i]
  pop2 <- Avg_popFst_data$pop2[i]
  
  km <- c(geodist_splittime_df$geoDist[geodist_splittime_df$species == sp & geodist_splittime_df$pop1 == pop1 & geodist_splittime_df$pop2 == pop2],
          geodist_splittime_df$geoDist[geodist_splittime_df$species == sp & geodist_splittime_df$pop1 == pop2 & geodist_splittime_df$pop2 == pop1])
  div <- c(geodist_splittime_df$splitTime[geodist_splittime_df$species == sp & geodist_splittime_df$pop1 == pop1 & geodist_splittime_df$pop2 == pop2],
           geodist_splittime_df$splitTime[geodist_splittime_df$species == sp & geodist_splittime_df$pop1 == pop2 & geodist_splittime_df$pop2 == pop1])
  
  Avg_popFst_data$dist[i] = ifelse(length(km) == 0, NA, km)
  Avg_popFst_data$div[i] = ifelse(length(div) == 0, NA, div)
  
}
Avg_popFst_data$guild <- factor(sapply(Avg_popFst_data$species, function(i) data$category[data$abbrev == i][1]), 
                                levels = guilds)
Avg_popFst_data$behavior <- factor(ifelse(Avg_popFst_data$guild == "Generalist", "Generalist", "Specialist"),
                                   levels = c("Generalist", "Specialist"))
Avg_popDxy_data$species <- factor(Avg_popDxy_data$species, levels=species)
Avg_popDxy_data$pop1 <- factor(Avg_popDxy_data$pop1, levels=popLevels)
Avg_popDxy_data$pop2 <- factor(Avg_popDxy_data$pop2, levels=popLevels)

lm_fst_geo_time <- lmerTest::lmer(fst ~ scale(dist) * scale(div) * behavior + (1 | species), data = Avg_popFst_data)
summary(lm_fst_geo_time)
anova(lm_fst_geo_time)

# Extract fixed-effect estimates for dist and div
fixed_effects <- nlme::fixef(lm_fst_geo_time)

pdf("./Figures/lm_fst_vs_geodistANDdivtime_residuals.pdf", width = 10, height = 7)
eval_model(lm_fst_geo_time)
dev.off()

sjPlot::tab_model(lm_fst_geo_time, auto.label = FALSE, file = "./Figures/lm_fst_vs_geodistANDdivtime_regression_table.html")

sjPlot::plot_model(lm_fst_geo_time, type = "pred", terms = c("dist", "div", "behavior"),  line.size = 2, colors = c("#EF8A62", "grey30", "#67A9CF")) +
  theme_light() +
  labs(x = "Geographical distance (km)", y = expression("Predicted" ~ F[ST]), title = "", 
       color = "Divergence Time") +
  coord_cartesian(ylim = c(0,1)) +
  theme(text = element_text(size = 18, color = 'black'), 
        axis.text.y = element_text(size = 16, color = 'black'),
        axis.text.x = element_text(size = 16, color = 'black'),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey85"), 
        panel.grid.minor.y = element_line(color = "grey85"),
        panel.border = element_rect(colour = "black", size = 1),
        legend.position = "top", 
        axis.title.y = element_text(vjust = -0.25, size = 18)) 
ggsave(paste0("./Figures/lm_Fst_vs_geodistANDdivtime.png"), width = 10, height = 9)

# Adjust Fst by removing the effects of dist and div
Avg_popFst_data <- Avg_popFst_data %>%
  dplyr::mutate(
    dist_effect = fixed_effects["scale(dist)"] * scale(dist),
    div_effect = fixed_effects["scale(div)"] * scale(div),
    adjusted_fst = fst - dist_effect - div_effect,
  )

##
# dxy
##
Avg_popDxy_data <- na.exclude(plyr::ldply(Avg_popDxy_List, .id = "species"))[,1:4]
names(Avg_popDxy_data)[4] = "Dxy"

for (i in 1:nrow(Avg_popDxy_data)) {
  sp <- Avg_popDxy_data$species[i]
  pop1 <- Avg_popDxy_data$pop1[i]
  pop2 <- Avg_popDxy_data$pop2[i]
  
  km <- c(geodist_splittime_df$geoDist[geodist_splittime_df$species == sp & geodist_splittime_df$pop1 == pop1 & geodist_splittime_df$pop2 == pop2],
          geodist_splittime_df$geoDist[geodist_splittime_df$species == sp & geodist_splittime_df$pop1 == pop2 & geodist_splittime_df$pop2 == pop1])
  div <- c(geodist_splittime_df$splitTime[geodist_splittime_df$species == sp & geodist_splittime_df$pop1 == pop1 & geodist_splittime_df$pop2 == pop2],
           geodist_splittime_df$splitTime[geodist_splittime_df$species == sp & geodist_splittime_df$pop1 == pop2 & geodist_splittime_df$pop2 == pop1])
  
  Avg_popDxy_data$dist[i] = ifelse(length(km) == 0, NA, km)
  Avg_popDxy_data$div[i] = ifelse(length(div) == 0, NA, div)
  
}
Avg_popDxy_data$guild <- factor(sapply(Avg_popDxy_data$species, function(i) data$category[data$abbrev == i][1]), 
                                levels = guilds)
Avg_popDxy_data$behavior <- factor(ifelse(Avg_popDxy_data$guild == "Generalist", "Generalist", "Specialist"),
                                   levels = c("Generalist", "Specialist"))

Avg_popDxy_data$species <- factor(Avg_popDxy_data$species, levels=species)
Avg_popDxy_data$pop1 <- factor(Avg_popDxy_data$pop1, levels=popLevels)
Avg_popDxy_data$pop2 <- factor(Avg_popDxy_data$pop2, levels=popLevels)

lm_Dxy_geo_time <- lmerTest::lmer(Dxy ~ scale(dist) * scale(div) * behavior + (1 | species), data = Avg_popDxy_data)
summary(lm_Dxy_geo_time)
anova(lm_Dxy_geo_time)

# Extract fixed-effect estimates for dist and div
fixed_effects <- nlme::fixef(lm_Dxy_geo_time)

pdf("./Figures/lm_Dxy_vs_geodistANDdivtime__residuals.pdf", width = 10, height = 7)
eval_model(lm_Dxy_geo_time)
dev.off()

sjPlot::tab_model(lm_Dxy_geo_time, auto.label = FALSE, file = "./Figures/lm_Dxy_vs_geodistANDdivtime_regression_table.html")

sjPlot::plot_model(lm_Dxy_geo_time, type = "pred", terms = c("dist", "div", "behavior"),  line.size = 2, colors = c("#EF8A62", "grey30", "#67A9CF")) +
  theme_light() +
  labs(x = "Geographical distance (km)", y = expression("Predicted" ~ d[xy]), title = "", 
       color = "Divergence Time") +
  coord_cartesian(ylim = c(0,1)) +
  theme(text = element_text(size = 18, color = 'black'), 
        axis.text.y = element_text(size = 16, color = 'black'),
        axis.text.x = element_text(size = 16, color = 'black'),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey85"), 
        panel.grid.minor.y = element_line(color = "grey85"),
        panel.border = element_rect(colour = "black", size = 1),
        legend.position = "top", 
        axis.title.y = element_text(vjust = -0.25, size = 18)) 
ggsave(paste0("./Figures/lm_Dxy_vs_geodistANDdivtime.png"), width = 10, height = 9)

# Adjust Dxy by removing the effects of dist and div
Avg_popDxy_data <- Avg_popDxy_data %>%
  dplyr::mutate(
    dist_effect = fixed_effects["scale(dist)"] * scale(dist),
    div_effect = fixed_effects["scale(div)"] * scale(div),
    adjusted_Dxy = Dxy - dist_effect - div_effect,
  )

##
# pi summary
##
Avg_popPi_data <- plyr::ldply(Avg_popPi_List, .id = "species")
Avg_popPi_data <- na.exclude(Avg_popPi_data)
colnames(Avg_popPi_data)[3] <- "pi"

Avg_popPi_data$guild <- factor(sapply(Avg_popPi_data$species, function(i) data$category[data$abbrev == i][1]), 
                                levels = guilds)
Avg_popPi_data$behavior <- factor(ifelse(Avg_popPi_data$guild == "Generalist", "Generalist", "Specialist"),
                                   levels = c("Generalist", "Specialist"))
Avg_popPi_data$species <- factor(Avg_popPi_data$species, levels=species)
Avg_popPi_data$pop <- factor(Avg_popPi_data$pop, levels=popLevels)

# Filter out populations with 1 individual
Avg_popPi_data$nsamples = sapply(1:nrow(Avg_popPi_data), function(i) 
  length(data$sample_ID[data$abbrev == Avg_popPi_data$species[i] & data$country_ID == Avg_popPi_data$pop[i]]))
Avg_popPi_data <- Avg_popPi_data[Avg_popPi_data$nsamples > 1,]

##############
# Boxplots
##############
##############
# Fst
##############
NFst_sp <- ggplot(Avg_popFst_data, aes(x = species, y = adjusted_fst, fill = guild)) +
  geom_point(data = Avg_popFst_data[Avg_popFst_data$species == "AKY",], shape = 21, size = 3) +
  geom_boxplot(data = Avg_popFst_data[Avg_popFst_data$species %in% c("CRP", "CLK", "AKA", "PRD", "MEL"),], lwd = 2, outlier.size = 3, outlier.color = NA) +
  geom_boxplot(data = Avg_popFst_data[Avg_popFst_data$species %in% c("CRP", "CLK", "AKA", "PRD", "MEL"),], fatten = 3, col = "white", outlier.color = NA) +
  geom_point(data = Avg_popFst_data[Avg_popFst_data$species == "POL",], shape = 21, size = 3) +
  geom_boxplot(data = Avg_popFst_data[Avg_popFst_data$species == "SAN",], lwd = 2, outlier.size = 3, outlier.color = NA) +
  geom_boxplot(data = Avg_popFst_data[Avg_popFst_data$species == "SAN",], fatten = 3, col = "white", outlier.color = NA) +
  geom_jitter(data = Avg_popFst_data[Avg_popFst_data$n_points > 1,], fill = "white", col = "black", shape = 21, size = 3) +
  scale_fill_manual("Host category", values = guildcols) +
  scale_x_discrete(
    breaks = levels(Avg_popDxy_data$species),
    labels = sapply(levels(Avg_popDxy_data$species), function(i) sub("Amphiprion_", "A. ", data$species[data$abbrev ==i][1]))
  ) +
  theme_light() + 
  coord_cartesian(ylim = c(0, 1)) +
  theme(text = element_text(size = 18, color = 'black'), 
        axis.text.y = element_text(size = 20, color = 'black'),
        axis.ticks.x = element_line(color = "black", linewidth = 2),
        axis.text.x = element_text(size = 20, color = 'black', face = "italic", angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey85"), 
        panel.grid.minor.y = element_line(color = "grey85"),
        plot.margin = ggplot2::margin(0.75,0.2,0.2,1, unit = "cm"),
        panel.border = element_rect(colour = "black", size = 1),
        legend.position = "none", 
        axis.title.y = element_text(vjust = 1, size = 22)) +
  labs(x = "", 
       y = expression("Adjusted" ~ italic(F[ST])))

Fst_class_inset <- ggplot(Avg_popFst_data, aes(x = behavior, y = adjusted_fst, fill = behavior)) +
  geom_boxplot(lwd = 2, outlier.size = 3, outlier.color = NA) +
  geom_boxplot(fatten = 3, col = "white", outlier.color = NA) +
  geom_jitter(fill = "white", col = "black", shape = 21, size = 3) +
  scale_fill_manual(values = c( "black", "grey")) +
  coord_cartesian(ylim = c(0, 0.65)) +
  theme_light() + 
  theme(text = element_text(size = 18, color = 'black'), 
        axis.text.y = element_text(size = 20, color = 'black'),
        axis.ticks.x = element_line(color = "black", linewidth = 2),
        axis.text.x = element_text(size = 20, color = 'black'),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey85"), 
        panel.grid.minor.y = element_line(color = "grey85"),
        plot.background = element_rect(fill = "transparent", colour = NA),
        plot.margin = ggplot2::margin(0.75,0.2,1.25,0, unit = "cm"),
        panel.border = element_rect(colour = "black", size = 1),
        legend.position = "none", 
        axis.title.y = element_text(vjust = 1, size = 22)) +
  labs(x = "", 
       y = "") + 
  stat_compare_means(comparisons = list(c("Generalist", "Specialist")), 
                     symnum.args = symnum.args, size = 8, label.y = 0.5)

kruskal.test(Avg_popFst_data$adjusted_fst, Avg_popFst_data$behavior)

NFst_sp <- NFst_sp +
  annotation_custom(
    grob = ggplotGrob(Fst_class_inset),
    xmin = 4,  # Adjust position
    xmax = 8.6,
    ymin = 0.4,
    ymax = 1.09
  )

##############
# dxy
##############
NDxy_sp <- ggplot(Avg_popDxy_data, aes(x = species, y = adjusted_Dxy, fill = guild)) +
  geom_point(data = Avg_popDxy_data[Avg_popDxy_data$species == "AKY",], shape = 21, size = 3) +
  geom_boxplot(data = Avg_popDxy_data[Avg_popDxy_data$species %in% c("CRP", "CLK", "AKA", "PRD", "MEL"),], lwd = 2, outlier.size = 3, outlier.color = NA) +
  geom_boxplot(data = Avg_popDxy_data[Avg_popDxy_data$species %in% c("CRP", "CLK", "AKA", "PRD", "MEL"),], fatten = 3, col = "white", outlier.color = NA) +
  geom_point(data = Avg_popDxy_data[Avg_popDxy_data$species == "POL",], shape = 21, size = 3) +
  geom_boxplot(data = Avg_popDxy_data[Avg_popDxy_data$species == "SAN",], lwd = 2, outlier.size = 3, outlier.color = NA) +
  geom_boxplot(data = Avg_popDxy_data[Avg_popDxy_data$species == "SAN",], fatten = 3, col = "white", outlier.color = NA) +
  geom_jitter(data = Avg_popDxy_data[Avg_popDxy_data$n_points > 1,], fill = "white", col = "black", shape = 21, size = 3) +
  scale_fill_manual("Host category", values = guildcols) +
  scale_x_discrete(
    breaks = levels(Avg_popDxy_data$species),
    labels = sapply(levels(Avg_popDxy_data$species), function(i) sub("Amphiprion_", "A. ", data$species[data$abbrev ==i][1]))
  ) +
  theme_light() + 
  coord_cartesian(ylim = c(0, 1)) +
  theme(text = element_text(size = 20, color = 'black'), 
        axis.text.y = element_text(size = 20, color = 'black'),
        axis.ticks.x = element_line(color = "black", linewidth = 2),
        axis.text.x = element_text(size = 20, color = 'black', face = "italic", angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey85"), 
        panel.grid.minor.y = element_line(color = "grey85"),
        plot.margin = ggplot2::margin(0.75,0.2,0.2,1, unit = "cm"),
        panel.border = element_rect(colour = "black", size = 1),
        legend.position = "none", 
        axis.title.y = element_text(vjust = 1, size = 22)) +
  labs(x = "", 
       y = expression("Adjusted" ~ italic(d[xy])))

Dxy_class_inset <- ggplot(Avg_popDxy_data, aes(x = behavior, y = adjusted_Dxy, fill = behavior)) +
  geom_boxplot(lwd = 2, outlier.size = 3, outlier.color = NA) +
  geom_boxplot(fatten = 3, col = "white", outlier.color = NA) +
  geom_jitter(fill = "white", col = "black", shape = 21, size = 3) +
  scale_fill_manual(values = c( "black", "grey")) +
  coord_cartesian(ylim = c(0, 0.65)) +
  theme_light() + 
  theme(text = element_text(size = 20, color = 'black'), 
        axis.text.y = element_text(size = 20, color = 'black'),
        axis.ticks.x = element_line(color = "black", linewidth = 2),
        axis.text.x = element_text(size = 20, color = 'black'),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey85"), 
        panel.grid.minor.y = element_line(color = "grey85"),
        plot.background = element_rect(fill = "transparent", colour = NA),
        plot.margin = ggplot2::margin(0.75,0.2,1.25,0, unit = "cm"),
        panel.border = element_rect(colour = "black", size = 1),
        legend.position = "none", 
        axis.title.y = element_text(vjust = 1, size = 22)) +
  labs(x = "", 
       y = "") + 
  stat_compare_means(comparisons = list(c("Generalist", "Specialist")), 
                     symnum.args = symnum.args, size = 8, label.y = 0.5)

kruskal.test(Avg_popDxy_data$adjusted_Dxy, Avg_popDxy_data$behavior)

NDxy_sp <- NDxy_sp +
  annotation_custom(
    grob = ggplotGrob(Dxy_class_inset),
    xmin = 4,  # Adjust position
    xmax = 8.6,
    ymin = 0.4,
    ymax = 1.09
  )

##############
## pi
##############

# check effect of number of samples
ggplot(Avg_popPi_data, aes(x = nsamples, y = pi, fill = category)) +
  geom_point() +
  theme_light() + 
  theme(text = element_text(size = 18, color = 'black'), 
        axis.text.y = element_text(size = 16, color = 'black'),
        axis.ticks.x = element_line(color = "black", linewidth = 2),
        axis.text.x = element_text(size = 16, color = 'black', face = "italic", angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey85"), 
        panel.grid.minor.y = element_line(color = "grey85"),
        plot.margin = ggplot2::margin(0.75,0.2,0.2,0.75, unit = "cm"),
        panel.border = element_rect(colour = "black", size = 1),
        legend.position = "none", 
        axis.title.y = element_text(vjust = 1, size = 18)) +
  labs(x = "Number of samples", 
       y = expression("Nucleotide diversity" ~ (pi)))

model <- lm(pi ~ nsamples, data = Avg_popPi_data)
summary(model) # no effect

Pi_sp <- ggplot(Avg_popPi_data, aes(x = species, y = pi, fill = category)) +
  geom_point(data = Avg_popPi_data[Avg_popPi_data$n_points == 1 & Avg_popPi_data$category == "Generalist",], shape = 21, size = 3) +
  geom_boxplot(data = Avg_popPi_data[Avg_popPi_data$n_points > 1 & Avg_popPi_data$category %in% c("Generalist", "RM specialist"),], lwd = 2, outlier.size = 3, outlier.color = NA) +
  geom_boxplot(data = Avg_popPi_data[Avg_popPi_data$n_points > 1 & Avg_popPi_data$category %in% c("Generalist", "RM specialist"),], fatten = 3, col = "white", outlier.color = NA) +
  geom_point(data = Avg_popPi_data[Avg_popPi_data$species == "EPH",], shape = 21, size = 3) +
  geom_boxplot(data = Avg_popPi_data[Avg_popPi_data$n_points > 1 & Avg_popPi_data$category %in% c("EQ specialist", "SD specialist"),], lwd = 2, outlier.size = 3, outlier.color = NA) +
  geom_boxplot(data = Avg_popPi_data[Avg_popPi_data$n_points > 1 & Avg_popPi_data$category %in% c("EQ specialist", "SD specialist"),], fatten = 3, col = "white", outlier.color = NA) +
  geom_jitter(data = Avg_popPi_data[Avg_popPi_data$n_points > 1,], fill = "white", col = "black", shape = 21, size = 3) +
  scale_fill_manual("Host category", values = guildcols) +
  scale_x_discrete(
    breaks = levels(Avg_popPi_data$species),
    labels = sapply(unique(Avg_popPi_data$species), function(i) sub("Amphiprion_", "A. ", data$species[data$abbrev ==i][1]))
  ) +
  coord_cartesian(ylim = c(0, 0.01)) + 
  theme_light() + 
  theme(text = element_text(size = 20, color = 'black'), 
        axis.text.y = element_text(size = 20, color = 'black'),
        axis.ticks.x = element_line(color = "black", linewidth = 2),
        axis.text.x = element_text(size = 20, color = 'black', face = "italic", angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey85"), 
        panel.grid.minor.y = element_line(color = "grey85"),
        plot.margin = ggplot2::margin(0.75,0.2,0.2,0.75, unit = "cm"),
        panel.border = element_rect(colour = "black", size = 1),
        legend.position = "none", 
        axis.title.y = element_text(vjust = 1, size = 22)) +
  labs(x = "", 
       y = expression("Nucleotide diversity" ~ (pi))) 

Pi_class_inset <- ggplot(Avg_popPi_data, aes(x = class, y = pi, fill = class)) +
  geom_boxplot(lwd = 2, outlier.size = 3, outlier.color = NA) +
  geom_boxplot(fatten = 3, col = "white", outlier.color = NA) +
  geom_jitter(fill = "white", col = "black", shape = 21, size = 3) +
  scale_fill_manual(values = c( "black", "grey")) +
  coord_cartesian(ylim = c(0, 0.0215)) +
  theme_light() + 
  coord_cartesian(ylim = c(0, 0.01)) + 
  theme(text = element_text(size = 20, color = 'black'), 
        axis.text.y = element_text(size = 20, color = 'black'),
        axis.ticks.x = element_line(color = "black", linewidth = 2),
        axis.text.x = element_text(size = 20, color = 'black'),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey85"), 
        panel.grid.minor.y = element_line(color = "grey85"),
        plot.background = element_rect(fill = "transparent", colour = NA),
        plot.margin = ggplot2::margin(0.75,0.2,1.25,0, unit = "cm"),
        panel.border = element_rect(colour = "black", size = 1),
        legend.position = "none", 
        axis.title.y = element_text(vjust = 1, size = 22)) +
  labs(x = "", 
       y = "")  +
  stat_compare_means(comparisons = list(c("Generalist", "Specialist")), symnum.args = symnum.args, 
                     size = 8, label.y = 0.008)

kruskal.test(Avg_popPi_data$pi, Avg_popPi_data$category)

Pi_sp <- Pi_sp +
  annotation_custom(
    grob = ggplotGrob(Pi_class_inset),
    xmin = 3.5,  # Adjust position
    xmax = 10.6,
    ymin = 0.004,
    ymax = 0.0109
  )

ggarrange(NFst_sp, NDxy_sp, Pi_sp, ncol = 3, 
          common.legend = TRUE, legend = "none", labels = letters, hjust = 0,
          font.label = list(size = 30, face = "bold", color ="black")) + bgcolor("white") + border(color = NA)
ggsave(paste0("./Figures/AvgNFST_AvgDxy_AvgPi_species_popgen_boxplot1.png"), width = 25, height = 9)

write.csv(Avg_popFst_data, "./analyses/popgen/FST_produced_dataset.csv", row.names = FALSE)
write.csv(Avg_popDxy_data, "./analyses/popgen/dxy_produced_dataset.csv", row.names = FALSE)
write.csv(Avg_popPi_data, "./analyses/popgen/Pi_produced_dataset.csv", row.names = FALSE)

###########
# Heatmaps with adjusted values
###########
pi_lim <- ceiling(max(Avg_popPi_data$pi) * 1e2) / 1e2

p_List <- list()
for (s in  species[!species %in% c("EPH", "LAT")]){
  
  pops <- popLevels[popLevels %in% data$country_ID[data$abbrev == s]]
  spname <- sub("Amphiprion_", "A. ", data$species[data$abbrev == s][1])
  
  sp_data <- merge(Avg_popDxy_data[Avg_popDxy_data$species == s, c("pop1", "pop2", "adjusted_Dxy")], 
                   Avg_popFst_data[Avg_popFst_data$species == s, c("pop1", "pop2", "adjusted_fst")], 
                   by = c("pop1", "pop2"))
  sp_data <- rbind(sp_data,
                   data.frame(pop1 = pops, 
                              pop2 = pops,
                              adjusted_Dxy = NA,
                              adjusted_fst = NA))
  
  sp_data$pop1 <- factor(sp_data$pop1, levels = pops)
  sp_data$pop2 <- factor(sp_data$pop2, levels = pops)
  
  p <- ggplot() +
    # Upper diagonal for Dxy
    geom_tile(data = sp_data,
              aes(pop1, pop2, fill = adjusted_Dxy), color = "white", linewidth = 0.5) +
    geom_text(data = sp_data,  size = 8,
              aes(pop1, pop2, label = ifelse(is.na(adjusted_Dxy), "NA", round(adjusted_Dxy, 2)), 
                  color = ifelse(is.na(adjusted_Dxy), "white", "black"))) +  
    scale_color_manual(values = c("white" = "white", "black" = "black")) +
    scale_fill_gradientn(colors = c("white", "#ED8F5F"), 
                         limits = c(0, 0.5),
                         name = expression(italic(d[xy])), 
                         guide = "none", 
                         na.value = "black") + 
    # Lower diagonal for FST
    new_scale_fill() +
    geom_tile(data = sp_data,
              aes(pop2, pop1, fill = adjusted_fst), color = "white", linewidth = 0.5) +
    geom_text(data = sp_data, size = 8,
              aes(pop2, pop1, label = ifelse(is.na(adjusted_fst), "NA", round(adjusted_fst, 2)), 
                  color = ifelse(is.na(adjusted_fst), "white", "black"))) +  
    scale_color_manual(values = c("white" = "white", "black" = "black")) +
    scale_fill_gradientn(colors = c("white", "violetred4"), 
                         limits = c(0, 0.5),
                         name = expression(F[ST]), 
                         guide = "none", 
                         na.value = "black") +
    # Diagonal for PI
    new_scale_fill() +
    geom_tile(data = Avg_popPi_data %>%
                filter(species == s),
              aes(pop, pop, fill = pi), color = "black", linewidth = 1.5) +
    geom_text(data = Avg_popPi_data %>%
                filter(species == s),  size = 8,
              aes(pop, pop, label = round(pi, 3)), color = "black") +  # Add FST values inside tiles
    scale_fill_gradientn(colors = c("white", "steelblue4"),
                         limits = c(0, pi_lim),
                         name = expression(pi),
                         guide = "none", 
                         na.value = "black") +
    scale_y_discrete(position = "right", labels = popnames[pops]) +
    scale_x_discrete(labels = popnames[pops]) +
    labs(title = spname, x = "", y = "") +
    theme_bw() +
    guides(color = "none") + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 14),
          axis.text.y = element_text(hjust = 1, size = 14),
          plot.title = element_text(size = 20, face = "italic", vjust = 1, hjust = 0.5),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_blank(),
          panel.background = element_rect(fill = "transparent"),  
          plot.margin = ggplot2::margin(0.5,0.5,0.5,0.5, unit = "cm"),
          plot.background = element_rect(fill = "white"),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16, face = "italic"))

  p_List[[s]] <- p
}

legend_dxy <- extract_legend("#ED8F5F", limits = c(0, 0.5), name = expression(italic(d[xy])), 
                             size = 20,  width = 15, height = 2)
legend_fst <- extract_legend("violetred4", limits = c(0, 0.5), name = expression(italic(F[ST])), 
                             size = 20, width = 15, height = 2)
legend_pi <- extract_legend("steelblue4", limits = c(0, pi_lim), name = expression(pi), 
                            size = 20, width = 15, height = 2)

# Arrange the legends side-by-side
grid_legends <- ggarrange(legend_dxy, legend_fst, legend_pi, ncol = 1) + bgcolor("white") + border("white")

pL <- ggarrange(p_List[["CLK"]], p_List[["CRP"]],
                p_List[["AKY"]], p_List[["MEL"]],
                p_List[["AKA"]], p_List[["PRD"]],
                p_List[["POL"]], p_List[["SAN"]],
                grid_legends, ncol = 3, nrow = 3, font.label = list(size = 26, face = "bold", color ="black"),
                labels = letters[1:8])
ggsave(paste0("./Figures/pi_dxy_Fst_adjusted.png"), plot = pL, width = 30, height = 30, limitsize = FALSE)


