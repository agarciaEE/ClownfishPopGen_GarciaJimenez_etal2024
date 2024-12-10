### EEMS plots
if(!("reemsplots2" %in% installed.packages()[,1])){
  remotes::install_github("dipetkov/reemsplots2")
}
library(reemsplots2)
library(raster)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(spatialreg)
library(spdep)

# working directory
setwd("~/Unil/Research/3_EcoPopGen/project_repository/")

path <- "./analyses/eems/"
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

data$lon[!is.na(data$lon) & data$lon < 0] = data$lon[!is.na(data$lon) & data$lon < 0]  + 360

pop_count <- data %>%
  group_by(abbrev, loc_ID) %>%
  dplyr::summarise(lon = mean(lon, na.rm = T),
                   lat = mean(lat, na.rm = T),
                   n = length(country))

# host categories
guilds <- c("Generalist", "RM specialist", "EQ specialist", "SD specialist")
# species vector
species <- c("AKA", "AKY", "CLK", "CRP", "MEL", "PRD", "POL", "SAN")
# population vector
popLevels <- c("KNY", "MYT", "MDG", "MDV", "THA", "IDN", "PHI", "TAI",
               "PAU", "PNG", "SOI", "AUS", "NCD", "FIJ", "POL")

# set guild colors
guildcols <- setNames(c("black", "orange", "red", "cyan"), guilds)

# Get world map data
world_map <- sf::st_as_sf(maps::map("world", wrap = c(0,360), plot = FALSE, fill = TRUE)) 

# Get major currents
currents <- read_sf("./raw_data/currents_maps/Major_Ocean_Currents/") %>%
  st_shift_longitude()

currents <- currents %>%
  filter(!currents$NAME %in% c("Guinea", "Norwegian") & !currents$OBJECTID %in% c(103, 104, 105)) 

# Get ITF currents
ITF <- raster("./raw_data/currents_maps/ITFcurrent.tif")
ITF[ITF>10] = NA
ITF[ITF>0] = 1
extent(ITF) <- c(84.5, 164, -29, 31.5)
poly <- rasterToPolygons(ITF)
poly <- st_as_sf(poly) %>%
  st_set_crs(crs(currents)) %>%
  st_simplify()

# log posteriors
######
chains <- 1:3
num.demes <- c(50, 200, 500)

for (s in species) {
  # Generate all combinations of chains and deme numbers
  combinations <- expand.grid(deme = num.demes, chain = chains)
  
  # Construct paths using the combinations
  mcmcpaths <- paste0(path, s, "_eemsFiles/", s, "-EEMS-ndemes", combinations$deme, "-chain", combinations$chain)
  
  plog <- plot_log_posterior(mcmcpaths)
  ggsave(paste0(path, s, "_all_log_posteriors.png"), plot = plog)
  
  plog_summ <- plog$data %>%
    group_by(path, chain, demes) %>%
    dplyr::summarise(mean_pilogl = mean(pilogl, na.rm = T),
                     sd_pilogl = sd(pilogl, na.rm = T))
  
  best_demes <- plog$data %>%
    dplyr::group_by(demes) %>%
    dplyr::summarise(mean_pilogl = mean(pilogl, na.rm = T),
                     sd_pilogl = sd(pilogl, na.rm = T)) %>%
    dplyr::arrange(desc(mean_pilogl)) %>%
    dplyr::slice(1) %>% dplyr::pull(demes)
    
  mcmcpaths <- plog_summ$path[plog_summ$demes == best_demes]

  plots <- make_eems_plots(mcmcpaths, add_grid = FALSE, col_grid = "grey25", longlat = TRUE)
  
  coords <- read.csv(paste0(path,  s, "_eemsFiles/", s, ".outer"), sep = " ", header = F)
  names(coords) <- c("LON", "LAT")
  
  guild <- data$category[data$abbrev == s][1]
  
  # Get species distribution
  dis_sf <- sf::read_sf(paste0("./analyses/sdm/", s, ".shp")) %>% st_buffer(dist = 100000)
  q_sf <- sf::st_as_sf(plots$mrates01$data, coords = c("x", "y"), crs = crs(dis_sf))
  outer_poly <- sf::st_difference(q_sf, dis_sf) %>%  sf::st_shift_longitude()
  
  theme_custom <- theme(plot.title = element_text(face = "italic", size = 16,  hjust = 0.5), legend.key.size = unit(1, "cm"),
                        panel.grid = element_blank(), panel.background = element_rect(fill = "grey92"), 
                        legend.key = element_rect(fill = "white"),
                        legend.text = element_text(size = 12), legend.title = element_text(size = 12))
  
  # plots
  m_plots[[s]] <- plots$mrates01 +
    geom_polygon(data = coords, aes(LON, LAT), fill = NA, col = "grey20", linetype = "dashed", linewidth = 1) +
    scale_x_continuous("", expand = c(0,0)) +
    scale_y_continuous("", expand = c(0,0)) +
    labs(title = gsub("_", " ", data$species[data$abbrev == s][1])) +
    geom_sf(data = currents, aes( x = NULL, y = NULL), fill = "white", col = "grey20", linewidth = 0.125, alpha = 0.75) +
    geom_sf(data = poly, aes( x = NULL, y = NULL), fill = "white", col = "grey20", linewidth = 0.125, alpha = 0.75) +
    geom_sf(data = world_map, aes( x = NULL, y = NULL), fill = "grey85", size = 1) +
    geom_point(data = pop_count[pop_count$abbrev == s,], 
               aes(x = lon, y = lat, size = n), shape = 21,  col = "grey20", fill = guildcols[guild], stroke = 1) +
    scale_size("# individuals", range = c(0.1,5), limits = c(0.1, max(pop_count$n))) +
    theme_light() +
    coord_sf(xlim = ext[1:2], ylim = ext[3:4]) +  
    theme_custom

  q_plots[[s]] <- plots$qrates01  +
    geom_sf(data = outer_poly, aes( x = NULL, y = NULL), fill ="grey92", col = "grey92") +
    scale_x_continuous("", expand = c(0,0)) +
    scale_y_continuous("", expand = c(0,0)) +
    labs(title = gsub("_", " ", data$species[data$abbrev == s][1])) +
    geom_sf(data = currents, aes( x = NULL, y = NULL), fill = "white", col = "grey20", linewidth = 0.125, alpha = 0.75) +
    geom_sf(data = poly, aes( x = NULL, y = NULL), fill = "white", col = "grey20", linewidth = 0.125, alpha = 0.75) +
    geom_sf(data = world_map, aes( x = NULL, y = NULL), fill = "grey85", size = 1) +
    geom_point(data = pop_count[pop_count$abbrev == s,], 
               aes(x = lon, y = lat, size = n), shape = 21,  col = "grey20", fill = guildcols[guild], stroke = 1) +
    scale_size("# individuals", range = c(0.1,5), limits = c(0.1, max(pop_count$n))) +
    theme_light() +
    coord_sf(xlim = ext[1:2], ylim = ext[3:4]) + 
    theme_custom

  Pm_plots[[s]] <- plots$mrates02 +
    geom_polygon(data = coords, aes(LON, LAT), fill = NA, col = "grey20", linetype = "dashed", linewidth = 1) +
    scale_x_continuous("", expand = c(0,0)) +
    scale_y_continuous("", expand = c(0,0)) +
    labs(title = gsub("_", " ", data$species[data$abbrev == s][1])) +
    geom_sf(data = currents, aes( x = NULL, y = NULL), fill = "white", col = "grey20", linewidth = 0.125, alpha = 0.75) +
    geom_sf(data = poly, aes( x = NULL, y = NULL), fill = "white", col = "grey20", linewidth = 0.125, alpha = 0.75) +
    geom_sf(data = world_map, aes( x = NULL, y = NULL), fill = "grey85", size = 1) +
    geom_point(data = pop_count[pop_count$abbrev == s,], 
               aes(x = lon, y = lat, size = n), shape = 21,  col = "grey20", fill = guildcols[guild], stroke = 1) +
    scale_size("# individuals", range = c(0.1,5), limits = c(0.1, max(pop_count$n))) +
    theme_light() +
    coord_sf(xlim = ext[1:2], ylim = ext[3:4]) +  
    theme_custom

  Pq_plots[[s]] <- plots$qrates02  +
    geom_sf(data = outer_poly, aes( x = NULL, y = NULL), fill ="grey92", col = "grey92") +
    scale_x_continuous("", expand = c(0,0)) +
    scale_y_continuous("", expand = c(0,0)) +
    labs(title = gsub("_", " ", data$species[data$abbrev == s][1])) +
    geom_sf(data = currents, aes( x = NULL, y = NULL), fill = "white", col = "grey20", linewidth = 0.125, alpha = 0.75) +
    geom_sf(data = poly, aes( x = NULL, y = NULL), fill = "white", col = "grey20", linewidth = 0.125, alpha = 0.75) +
    geom_sf(data = world_map, aes( x = NULL, y = NULL), fill = "grey85", size = 1) +
    geom_point(data = pop_count[pop_count$abbrev == s,], 
               aes(x = lon, y = lat, size = n), shape = 21,  col = "grey20", fill = guildcols[guild], stroke = 1) +
    scale_size("# individuals", range = c(0.1,5), limits = c(0.1, max(pop_count$n))) +
    theme_light() +
    coord_sf(xlim = ext[1:2], ylim = ext[3:4]) + 
    theme_custom

  theme_custom <- theme_bw() + theme(text = element_text(size = 16),  plot.margin = unit(c(1,1,1,1), "cm"), panel.border = element_rect(linewidth = 1))
  ggsave(paste0(path, s, "_best_model_obs_vs_fitted_genetic_dissimilarity_between_demes.png"),
         plot = plots$rdist01 + geom_point(size=5) + theme_custom )
  
  ggsave(paste0(path, s, "_best_model_obs_vs_fitted_genetic_dissimilarity_within_demes.png"),
         plot = plots$rdist02 + geom_point(size=5) + theme_custom)
  
  ggsave(paste0(path, s, "_best_model_eems_ibd.png"), plot = plots$rdist03 + geom_point(size=5) + theme_custom)

  xym.values <- as.data.frame(plots$mrates01$data)
  xym.Pvalues <- as.data.frame(plots$mrates02$data)
  xyq.values <- as.data.frame(plots$qrates01$data)
  xyq.Pvalues <- as.data.frame(plots$qrates02$data)
  
  xyqm_df <- data.frame(xym.values[,1:2], 
                        q = xyq.values[,3],
                        Pq = xyq.Pvalues[,3],
                        m = xym.values[,3],
                        Pm = xym.Pvalues[,3])
  
  write.csv(xyqm_df, file = paste0(path, "generated_datasets/", s, "_eems_data.csv"), row.names = FALSE)
}

dir.create("./Figures", showWarnings = FALSE)

# Plot Figures
pq <- ggarrange(Pq_plots$CLK, Pq_plots$AKY, Pq_plots$CRP,
                Pq_plots$MEL,
                Pq_plots$AKA, Pq_plots$PRD,
                Pq_plots$POL, Pq_plots$SAN,
                font.label = list(size = 20, face = "bold", color ="black"),
                nrow = 4, ncol = 2, common.legend = TRUE, labels = letters[1:10], legend = "right") +
  bgcolor("white") + border(color = NA)
ggsave(paste0("./Figures/all_species_Pq.png"), plot = pq, height = 15, width = 20)

q <- ggarrange(q_plots$CLK, q_plots$AKY, q_plots$CRP,
               q_plots$MEL,
               q_plots$AKA, q_plots$PRD,
               q_plots$POL, q_plots$SAN,
               font.label = list(size = 20, face = "bold", color ="black"),
               nrow = 4, ncol = 2, common.legend = TRUE, labels = letters[1:10], legend = "right") +
  bgcolor("white") + border(color = NA) 
ggsave(paste0("./Figures/all_species_q.png"), plot = q, height = 15, width = 20)

pm <- ggarrange(Pm_plots$CLK, Pm_plots$AKY, Pm_plots$CRP,
                Pm_plots$MEL,
                Pm_plots$AKA, Pm_plots$PRD,
                Pm_plots$POL, Pm_plots$SAN,
                font.label = list(size = 20, face = "bold", color ="black"),
                nrow = 4, ncol = 2, common.legend = TRUE, labels = letters[1:10], legend = "right") +
  bgcolor("white") + border(color = NA)
ggsave(paste0("./Figures/all_species_Pm.png"), plot = pm, height = 15, width = 20)

m <- ggarrange(m_plots$CLK, m_plots$AKY, m_plots$CRP,
               m_plots$MEL,
               m_plots$AKA, m_plots$PRD,
               m_plots$POL, m_plots$SAN,
               font.label = list(size = 20, face = "bold", color ="black"),
               nrow = 4, ncol = 2, common.legend = TRUE, labels = letters[1:10], legend = "right") +
  bgcolor("white") + border(color = NA)
ggsave(paste0("./Figures/all_species_m.png"), plot = m, height = 15, width = 20)

