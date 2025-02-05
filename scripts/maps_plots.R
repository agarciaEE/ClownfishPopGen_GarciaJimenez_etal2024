if(!("plotmaps" %in% installed.packages()[,1])){
  remotes::install_github("halasadi/plotmaps")
}
library(plotmaps)
library(sf)
library(raster)
library(tidyverse)
library(ggpubr)
library(dplyr)
library(purrr)
library(lme4)
library(lmerTest)
library(plyr)

# working directory
setwd("~/Unil/Research/3_EcoPopGen/project_repository/")

### Custom Functions
source("./scripts/custom_functions.R")

path <- "./analyses/maps/"

#####
## Load samples info
#####
data <- read.csv("./data/samples_info.csv")

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
world_lgm_map <- sf::st_read("./raw_data/world_clgm/world_cut.shp") %>%
  sf::st_set_crs("+init=epsg:4326")

# Get Continental maps only to overlay in plots
# List of countries/regions to remove
remove_list <- c("Indonesia", "Philippines", "Malaysia", "Papua New Guinea", 
                 "Solomon Islands", "Fiji", "Vanuatu", "New Caledonia", 
                 "Micronesia", "Marshall Islands", "Kiribati", "Tuvalu", 
                 "Samoa", "Tonga", "Cook Islands", "French Polynesia")

# Filter out the listed countries/regions
continental_map <- world_map[!world_map$ID %in% remove_list, ]
continental_map <- sf::st_make_valid(continental_map)
# Calculate area of each polygon
continental_map$area <- sf::st_area(continental_map)
# Set a threshold 
area_threshold <- 100000 * 1000000  # 100,000 km^2 in m^2
# Filter out small islands
continental_map <- continental_map[as.numeric(continental_map$area) > area_threshold, ]

# Get major currents
currents <- read_sf("./data/Major_Ocean_Currents/") %>%
  st_shift_longitude()

currents <- currents %>%
  filter(!currents$NAME %in% c("Guinea", "Norwegian") & !currents$OBJECTID %in% c(103, 104, 105)) 

# Get ITF currents
ITF <- raster("./data/currents_maps/ITFcurrent.tif")
ITF[ITF>10] = NA
ITF[ITF>0] = 1
extent(ITF) <- c(84.5, 164, -29, 31.5)
poly <- rasterToPolygons(ITF)
poly <- st_as_sf(poly) %>%
  st_set_crs(crs(currents)) %>%
  st_simplify()

# host category comparisons
mycomp <-  combn(levels(data$category), 2)
mycomp <- lapply(1:ncol(mycomp), function(i) mycomp[,i])
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "n.s."))

# set extent
sdm_exts <- sapply(amph_SDM, function(x) as.vector(extent(x)))
ext <- c(min(sdm_exts[1,]), max(sdm_exts[2,]), min(sdm_exts[3,]), max(sdm_exts[4,]))
ext <- c(30, 215, -35, 35)

# custom theme
theme_custom <- theme(plot.title = element_text(face = "italic", size = 16,  hjust = 0.5), legend.key.size = unit(1, "cm"),
                      panel.grid = element_blank(), panel.background = element_rect(fill = "grey92"), 
                      legend.key = element_rect(fill = "white"),
                      legend.text = element_text(size = 12), legend.title = element_text(size = 12))

# parameters
######
chains <- 1:3
num.demes <- c(50, 200, 500)
ibd_lengths <- c("0.02-1", "2-6", "6-Inf")
gT = 5 # generation Time
width = 20
height = 12
recombMap = TRUE
correct.by.area = TRUE
save.plot = TRUE
common_limits = TRUE
set.range = TRUE
mlim <- c(0,1e4)
Nlim <- c(0,1e4)
best_ndemes = TRUE

###########################
# ----------------------- #
#     log posteriors      #
# ----------------------- #
###########################
best_deme_bysp <- data.frame()
for (s in species){
  plog_sp <- list()
  for (ibd in ibd_lengths){
    
    sp_path = paste0(path, s)

    mcmcpaths <- paste0(sp_path, "/", grep("MAPS", list.dirs(sp_path, recursive = FALSE, full.names = FALSE), value = TRUE))
    mcmcpaths <- grep(paste0("-", sub("-", "_", ibd)), mcmcpaths, value = TRUE)
    
    success_runs <- sapply(lapply(mcmcpaths, list.files), length) > 4
    n_success_runs <- sum(success_runs)
    if (n_success_runs != length(mcmcpaths)){
      warning("Not all runs where successful.")
    }
    mcmcpaths <- mcmcpaths[success_runs]
    
    # get time period based on IBD segment length
    l = strsplit(as.character(ibd), "\\-")[[1]]
    l1 = as.numeric(l[1])
    l2 = as.numeric(l[2])
    Y <- ibdL2nG(l1, l2, gT)
    
    plog <- plot_trace(mcmcpaths) + 
      labs(title = paste("IBD length:", ibd, "(~", Y, "years ago)"))
    
    # add tol plot list
    plog_sp[[ibd]] <- plog
    
    # get summary of msmc logliks
    plog_summ <- plog$data %>%
      group_by(path, chain, demes) %>%
      dplyr::summarise(mean_pilogl = mean(pilogl, na.rm = T),
                       sd_pilogl = sd(pilogl, na.rm = T))
    
    # write down
    write.csv(plog_summ, file = paste0(path, "/", prefix, "log_posterior_mcmc_summary.csv"))
    
    # get best ndemes based on mean loglik across chains
    plog_bydeme_summ <- plog_summ %>%
      group_by(demes) %>%
      dplyr::summarise(mean_pilogl = mean(mean_pilogl, na.rm = T),
                       sd_pilogl = sd(sd_pilogl, na.rm = T))
    
    ndemes <- plog_bydeme_summ$demes[which.max(plog_bydeme_summ$mean_pilogl)]
    
    # add to dataset
    best_deme_bysp <- rbind(best_deme_bysp, data.frame(species = prefix, 
                                                       IBDseg = ibd,
                                                       ndemes = ndemes))
    
  }
  ggarrange(plog_sp$`0.02-1`, plog_sp$`2-6`, plog_sp$`6-Inf`, nrow = 3, common.legend = TRUE, legend = "right")
  ggsave(paste0(path, "/", prefix, "log_posterior_mcmc_.png"), height = 15, width = 10)
  
}
write.csv(best_deme_bysp, file = paste0(path, "best_deme_bysp.csv"), row.names = FALSE)

###########################
###########################
# ----------------------- #
#       PLOT MAPS         #
# ----------------------- #
###########################
best_deme_bysp <- read.csv(paste0(path, "best_deme_bysp.csv"))
N_dataset <- m_dataset <- list()
for (ibd in ibd_lengths){
  
  m_plots <- list()
  N_plots <- list()
  N_dataset[[ibd]] <- m_dataset[[ibd]]  <- list()
  
  for (prefix in species){
    
    ndemes <- best_deme_bysp %>%
      dplyr::filter(species == prefix & IBDseg == ibd) %>%
      dplyr::pull(ndemes)
    
    # get MAPS output
    #################
    sp_path = paste0(path, prefix, "/")
    
    host.category <- data$category[data$abbrev == prefix][1]
    
    # output path indicating the cM range taken
    outpath <- paste0(ibd, "_", ndemes, "demes")
    
    # mcmcpaths: choosing all three chains of best ndemes
    mcmcpaths <- paste0(sp_path, "/", grep("MAPS", list.dirs(sp_path, recursive = FALSE, full.names = FALSE), value = TRUE))
    mcmcpaths <- grep(paste0("ndemes", ndemes, "-", sub("-", "_", ibd)), mcmcpaths, value = TRUE)

    if (any(sapply(sapply(mcmcpaths, list.files), length) > 0)) {
      
      ## get MAPS plots
      p <- plot_MAPS(add.pts = FALSE, add.graph = FALSE, add.countries = FALSE,
                     longlat = TRUE, mcmcpath = mcmcpaths, correct.by.area = correct.by.area,
                     set.range = TRUE, m.limits = mlim, N.limits = Nlim, 
                     outpath = outpath, width = width, height = height, trans = "log10", 
                     save.plot = FALSE)
      
      # extract data
      mdata<- cbind(p$m$layers[[1]]$data, P = p$m_sign$layers[[1]]$data$P)
      Ndata <- cbind(p$N$layers[[1]]$data, P = p$N_sign$layers[[1]]$data$P)
      
      Ndata <- Ndata[Ndata$filter == TRUE, c("x", "y", "ss", "upper.ci", "lower.ci", "P")]
      colnames(Ndata) <- c("x", "y", "De", "upper.ci", "lower.ci", "P")
      
      mdata <- mdata[mdata$filter == TRUE, c("x", "y", "ss", "upper.ci", "lower.ci", "P")]
      colnames(mdata) <- c("x", "y", "sigma", "upper.ci", "lower.ci", "P")
      
      write.csv(Ndata, file = paste0(path, "generated_datasets/", 
                                     prefix, "_De_", ibd, "-", ndemes, ".csv"))
      write.csv(mdata, file = paste0(path, "generated_datasets/",
                                     prefix, "_sigma_", ibd, "-", ndemes, ".csv"))
      
      ### add features and custom themes to plots
      # get individual coords
      coords <- read.csv(paste0(sp_path, prefix, ".outer"), sep = " ", header = F)
      names(coords) <- c("LON", "LAT")
      
      # get distribution and outer polygon
      dis_sf <- sf::read_sf(paste0("./analyses/sdm/", s, ".shp")) %>% st_buffer(dist = 100000)
      q_sf <- st_as_sf(p$m$layers[[1]]$data, coords = c("x", "y"), crs = crs(dis_sf))
      outer_poly <-  st_difference(q_sf, dis_sf) %>%  st_shift_longitude()
      
      p$m <- p$m +
        geom_polygon(data = coords, aes(LON, LAT), fill = NA, col = "grey20", linetype = "dashed", linewidth = 1) +
        scale_x_continuous("", expand = c(0,0)) +
        scale_y_continuous("", expand = c(0,0)) +
        labs(title = gsub("_", " ", data$species[data$abbrev == prefix][1])) +
        geom_sf(data = currents, aes( x = NULL, y = NULL), fill = "white", col = "grey20", linewidth = 0.125, alpha = 0.75) +
        geom_sf(data = ITF_pol, aes( x = NULL, y = NULL), fill = "white", col = "black", linewidth = 0.125, alpha = 0.75) +
        geom_sf(data = world_map, aes( x = NULL, y = NULL), fill = "grey85", size = 1) +
        p$m$layers[[3]] + 
        geom_sf(data = continental_map, aes( x = NULL, y = NULL), fill = "grey85", size = 1) +
        geom_point(data = pop_count[pop_count$abbrev == prefix,], 
                   aes(x = lon, y = lat, size = n), shape = 21,  col = "grey20", fill = loc.cols[host.category], stroke = 1) +
        scale_size("# individuals", range = c(0.1,5), limits = c(0.1, max(pop_count$n))) +
        theme_light() +
        coord_sf(xlim = ext[1:2], ylim = ext[3:4]) +  # Replace coord_cartesian with coord_sf
        theme_custom
      if (ibd == "0.02-1") {
        p$m <- p$m + geom_sf(data = world_lgm_map, aes( x = NULL, y = NULL), fill = "grey70", col = "grey70", size = 1) +
          geom_sf(data = world_map, aes( x = NULL, y = NULL), fill = "grey85", size = 1) +
          p$m$layers[[3]] + 
          geom_point(data = pop_count[pop_count$abbrev == prefix,], 
                     aes(x = lon, y = lat, size = n), shape = 21,  col = "grey20", fill = loc.cols[host.category], stroke = 1) +
          coord_sf(xlim = ext[1:2], ylim = ext[3:4])
      }
      p$N <- p$N  +
        geom_sf(data = outer_poly, aes( x = NULL, y = NULL), fill ="grey92", col = "grey92", size = 2.5) +
        scale_x_continuous("", expand = c(0,0)) +
        scale_y_continuous("", expand = c(0,0)) +
        labs(title = gsub("_", " ", data$species[data$abbrev == prefix][1])) +
        geom_sf(data = currents, aes( x = NULL, y = NULL), fill = "white", col = "grey20", linewidth = 0.125, alpha = 0.75) +
        geom_sf(data = ITF_pol, aes( x = NULL, y = NULL), fill = "white", col = "black", linewidth = 0.125, alpha = 0.75) +
        geom_sf(data = world_map, aes( x = NULL, y = NULL), fill = "grey85", size = 1) +
        p$N$layers[[3]] + 
        geom_point(data = pop_count[pop_count$abbrev == prefix,], 
                   aes(x = lon, y = lat, size = n), shape = 21,  col = "grey20", fill = loc.cols[host.category], stroke = 1) +
        geom_sf(data = continental_map, aes( x = NULL, y = NULL), fill = "grey85", size = 1) +
        scale_size("# individuals", range = c(0.1,5), limits = c(0.1, max(pop_count$n))) +
        theme_light() +
        coord_sf(xlim = ext[1:2], ylim = ext[3:4]) +  # Replace coord_cartesian with coord_sf
        theme_custom
      if (ibd == "0.02-1") {
        p$N <- p$N  + geom_sf(data = world_lgm_map, aes( x = NULL, y = NULL), fill = "grey70", col = "grey70", size = 1) +
          geom_sf(data = world_map, aes( x = NULL, y = NULL), fill = "grey85", size = 1) +
          p$N$layers[[3]] + 
          geom_point(data = pop_count[pop_count$abbrev == prefix,], 
                     aes(x = lon, y = lat, size = n), shape = 21,  col = "grey20", fill = loc.cols[host.category], stroke = 1) +
          coord_sf(xlim = ext[1:2], ylim = ext[3:4])
      }
      p$m_sign <- p$m_sign  +
        geom_polygon(data = coords, aes(LON, LAT), fill = NA, col = "grey20", linetype = "dashed", linewidth = 1) +
        scale_x_continuous("", expand = c(0,0)) +
        scale_y_continuous("", expand = c(0,0)) +
        labs(title = gsub("_", " ", data$species[data$abbrev == prefix][1])) +
        geom_sf(data = currents, aes( x = NULL, y = NULL), fill = "white", col = "grey20", linewidth = 0.125, alpha = 0.75) +
        geom_sf(data = ITF_pol, aes( x = NULL, y = NULL), fill = "white", col = "black", linewidth = 0.125, alpha = 0.75) +
        geom_sf(data = world_map, aes( x = NULL, y = NULL), fill = "grey85", size = 1) +
        geom_point(data = pop_count[pop_count$abbrev == prefix,], 
                   aes(x = lon, y = lat, size = n), shape = 21,  col = "grey20", fill = loc.cols[host.category], stroke = 1) +
        scale_size("# individuals", range = c(0.1,5), limits = c(0.1, max(pop_count$n))) +
        theme_light() +
        coord_sf(xlim = ext[1:2], ylim = ext[3:4]) +  # Replace coord_cartesian with coord_sf
        theme_custom
      if (ibd == "0.02-1") {
        p$m_sign <- p$m_sign  + geom_sf(data = world_lgm_map, aes( x = NULL, y = NULL), fill = "grey70", col = "grey70", size = 1) +
          geom_sf(data = world_map, aes( x = NULL, y = NULL), fill = "grey85", size = 1) +
          geom_point(data = pop_count[pop_count$abbrev == prefix,], 
                     aes(x = lon, y = lat, size = n), shape = 21,  col = "grey20", fill = loc.cols[host.category], stroke = 1) +
          coord_sf(xlim = ext[1:2], ylim = ext[3:4])
      }
      p$N_sign <- p$N_sign  +
        geom_sf(data = outer_poly, aes( x = NULL, y = NULL), fill ="grey92", col = "grey92", size = 2.5) +
        scale_x_continuous("", expand = c(0,0)) +
        scale_y_continuous("", expand = c(0,0)) +
        labs(title = gsub("_", " ", data$species[data$abbrev == prefix][1])) +
        geom_sf(data = currents, aes( x = NULL, y = NULL), fill = "white", col = "grey20", linewidth = 0.125, alpha = 0.75) +
        geom_sf(data = ITF_pol, aes( x = NULL, y = NULL), fill = "white", col = "black", linewidth = 0.125, alpha = 0.75) +
        geom_sf(data = world_map, aes( x = NULL, y = NULL), fill = "grey85", size = 1) +
        geom_point(data = pop_count[pop_count$abbrev == prefix,], 
                   aes(x = lon, y = lat, size = n), shape = 21,  col = "grey20", fill = loc.cols[host.category], stroke = 1) +
        scale_size("# individuals", range = c(0.1,5), limits = c(0.1, max(pop_count$n))) +
        theme_light() +
        coord_sf(xlim = ext[1:2], ylim = ext[3:4]) +  # Replace coord_cartesian with coord_sf
        theme_custom
      if (ibd == "0.02-1") {
        p$N_sign <- p$N_sign  + geom_sf(data = world_lgm_map, aes( x = NULL, y = NULL), fill = "grey70", col = "grey70", size = 1) +
          geom_sf(data = world_map, aes( x = NULL, y = NULL), fill = "grey85", size = 1) +
          geom_point(data = pop_count[pop_count$abbrev == prefix,], 
                     aes(x = lon, y = lat, size = n), shape = 21,  col = "grey20", fill = loc.cols[host.category], stroke = 1) +
          coord_sf(xlim = ext[1:2], ylim = ext[3:4])
      }
      
      m_plots[[prefix]] <- list(m = p$m, m_sign = p$m_sign)
      N_plots[[prefix]] <- list(N = p$N, N_sign = p$N_sign)
    }
  }
  
  pN <- ggarrange(N_plots$CLK$N_sign, N_plots$AKY$N_sign, N_plots$CRP$N_sign,
                  N_plots$MEL$N_sign,
                  N_plots$AKA$N_sign, N_plots$PRD$N_sign,
                  N_plots$POL$N_sign, N_plots$SAN$N_sign,
                  font.label = list(size = 20, face = "bold", color ="black"),
                  nrow = 4, ncol = 2, common.legend = TRUE, labels = letters[1:10], legend = "right") +
    bgcolor("white") + border(color = NA)
  ggsave(paste0("./Figures/all_species_", ifelse(recombMap, "GenMap_", "NoGenMap_"), "BestNdemes_", ibd, "_PN.png"), plot = pN, height = 15, width = 20)
  
  N <- ggarrange(N_plots$CLK$N, N_plots$AKY$N, N_plots$CRP$N,
                 N_plots$MEL$N,
                 N_plots$AKA$N, N_plots$PRD$N,
                 N_plots$POL$N, N_plots$SAN$N,
                 font.label = list(size = 20, face = "bold", color ="black"),
                 nrow = 4, ncol = 2, common.legend = TRUE, labels = letters[1:10], legend = "right") +
    bgcolor("white") + border(color = NA) 
  ggsave(paste0("./Figures/all_species_", ifelse(recombMap, "GenMap_", "NoGenMap_"), "BestNdemes_", ibd, "_N.png"), plot = N, height = 15, width = 20)
  
  pm <- ggarrange(m_plots$CLK$m_sign, m_plots$AKY$m_sign, m_plots$CRP$m_sign,
                  m_plots$MEL$m_sign,
                  m_plots$AKA$m_sign, m_plots$PRD$m_sign,
                  m_plots$POL$m_sign, m_plots$SAN$m_sign,
                  font.label = list(size = 20, face = "bold", color ="black"),
                  nrow = 4, ncol = 2, common.legend = TRUE, labels = letters[1:10], legend = "right") +
    bgcolor("white") + border(color = NA)
  ggsave(paste0("./Figures/all_species_", ifelse(recombMap, "GenMap_", "NoGenMap_"), "BestNdemes_", ibd, "_Pm_abs.png"), plot = pm, height = 15, width = 20)
  
  m <- ggarrange(m_plots$CLK$m, m_plots$AKY$m, m_plots$CRP$m,
                 m_plots$MEL$m,
                 m_plots$AKA$m, m_plots$PRD$m,
                 m_plots$POL$m, m_plots$SAN$m,
                 font.label = list(size = 20, face = "bold", color ="black"),
                 nrow = 4, ncol = 2, common.legend = TRUE, labels = letters[1:10], legend = "right") +
    bgcolor("white") + border(color = NA)
  ggsave(paste0("./Figures/all_species_", ifelse(recombMap, "GenMap_", "NoGenMap_"), "BestNdemes_", ibd, "_m_abs.png"), plot = m, height = 15, width = 20)
  
}

###########################
# ----------------------- #
#      BOXPLOT MAPS       #
# ----------------------- #
###########################
best_deme_bysp <- read.csv(paste0(path, "best_deme_bysp.csv"))

N_dataset_selected <- m_dataset_selected <- list()
for (ibdseg in ibd_lengths){
  
  N_dataset_selected[[ibdseg]] <- m_dataset_selected[[ibdseg]] <- list()
  
  for (prefix in species) {
    
    ndemes <- best_deme_bysp %>%
      dplyr::filter(species == prefix & IBDseg == ibdseg) %>%
      dplyr::pull(ndemes)
    
    ndemes <- paste0(ndemes, "demes")
    
    N_dataset_selected[[ibdseg]][[prefix]] <- read.csv(paste0(path, "generated_datasets/", 
                                                              prefix, "_De_", ibdseg, "-", ndemes, ".csv"), row.names = 1)
    m_dataset_selected[[ibdseg]][[prefix]] <- read.csv(paste0(path, "generated_datasets/", 
                                                              prefix, "_sigma_", ibdseg, "-", ndemes, ".csv"), row.names = 1)
  }
  
  N_dataset_selected[[ibdseg]] <- plyr::ldply(N_dataset_selected[[ibdseg]], .id = "species")
  m_dataset_selected[[ibdseg]] <- plyr::ldply(m_dataset_selected[[ibdseg]], .id = "species")
}

N_dataset_selected <- ldply(N_dataset_selected, .id = "ibdLength")
m_dataset_selected <- ldply(m_dataset_selected, .id = "ibdLength")

# add guilds info
for (guild in guilds){
  sps <- unique(data$abbrev[data$category == guild])
  N_dataset_selected$guild[N_dataset_selected$species %in% sps] = guild
  m_dataset_selected$guild[m_dataset_selected$species %in% sps] = guild
}

# set factors
N_dataset_selected$guild <- factor(N_dataset_selected$guild, levels = guilds)
m_dataset_selected$guild <- factor(m_dataset_selected$guild, levels = guilds)

m_dataset_selected$species <- factor(m_dataset_selected$species, levels = species)
N_dataset_selected$species <- factor(N_dataset_selected$species, levels = species)

# create generaslit vs specialist factor
N_dataset_selected$behavior <- as.character(N_dataset_selected$guild)
N_dataset_selected$behavior[N_dataset_selected$behavior != "Generalist"] = "Specialist"
N_dataset_selected$behavior <- factor(N_dataset_selected$behavior, levels = c("Generalist", "Specialist"))

m_dataset_selected$behavior <- as.character(m_dataset_selected$guild)
m_dataset_selected$behavior[m_dataset_selected$behavior != "Generalist"] = "Specialist"
m_dataset_selected$behavior <- factor(m_dataset_selected$behavior, levels = c("Generalist", "Specialist"))

# add time periods
for (ibd in unique(N_dataset_selected$ibdLength)) {
  l = strsplit(as.character(ibd), "\\-")[[1]]
  l1 = as.numeric(l[1])
  l2 = as.numeric(l[2])
  ya <- ibdL2nG(l1, l2, gT)
  
  N_dataset_selected$time[N_dataset_selected$ibdLength == ibd] = ya
  m_dataset_selected$time[m_dataset_selected$ibdLength == ibd] = ya
  
}

# get marine regions
regfile <- "./raw_data/MEOW_ECOS/"
marine_regions <- sf::st_read(regfile)
marine_regions <- sf::st_transform(marine_regions, "+init=epsg:4326") # transform to WGS84

# subset marine_regions of interest
marine_regions <- marine_regions[marine_regions$REALM %in% c("Central Indo-Pacific", "Western Indo-Pacific", "Eastern Indo-Pacific",
                                                             "Temperate Australasia") | marine_regions$ECOREGION %in% c("East China Sea", "Central Kuroshio Current"),]
marine_regions <- marine_regions[!marine_regions$PROVINCE %in% "Easter Island",]
marine_regions <- marine_regions[!marine_regions$PROVINCE %in% unique(marine_regions$PROVINCE[marine_regions$REALM == "Temperate Australasia" & !marine_regions$PROVINCE %in% c("West Central Australian Shelf", "East Central Australian Shelf")]),]

plot(marine_regions["ECOREGION"])

# get delimited marine_regions representative of each species population
pop_count <- na.exclude(pop_count)
pop_count$region = NA
# Loop to assign region names to each point in pop_count
for (i in 1:nrow(pop_count)) {
  
  # Create sf point geometry for each location
  coords <- sf::st_sfc(sf::st_point(c(pop_count$lon[i], pop_count$lat[i])), crs = 4326)
  
  # Find intersecting region and assign it to pop_count$region
  ovr <- sf::st_intersects(marine_regions, coords, sparse = FALSE)
  
  # Assign the corresponding ecoregion name, if an intersection is found
  pop_count$region[i] <- if (any(ovr)) marine_regions$ECOREGION[ovr] else NA
}

world_map <- st_transform(world_map, st_crs(marine_regions))
world_lgm_map <- st_transform(world_lgm_map, st_crs(marine_regions))

marine_regions <- st_make_valid(marine_regions)
world_map <- st_make_valid(world_map)
world_lgm_map <- st_make_valid(world_lgm_map)

marine_regions <- st_difference(marine_regions, st_union(world_map))
marine_regions_lgm <- st_difference(marine_regions, st_union(world_lgm_map))

# get time points
time_points <- sort(unique(m_dataset_selected$time))

# duplicate species column for processing
N_dataset_selected$abbrev = as.character(N_dataset_selected$species)

# Process the dataset to get average values per identified regions
N_pop_summary <- N_dataset_selected %>%
  group_by(species, ibdLength) %>%
  group_modify(~ {
    
    sp <- unique(.x$abbrev)
    spname <- data$species[data$abbrev == sp][1]
    
    # Retrieve marine_regions for the specific species
    sp_regions <- pop_count$region[pop_count$abbrev == sp]
    
    # Filter for relevant polygons in marine_regions
    if (unique(.x$time) == 19125) {
      sp_reg_poly <- marine_regions_lgm %>%
        filter(ECOREGION %in% sp_regions)
    } else {
      sp_reg_poly <- marine_regions %>%
        filter(ECOREGION %in% sp_regions)
    }
    # transform subset MAPS dataset to sf object
    df_sf <- st_as_sf(
      data.frame(x = .x$x, y = .x$y, value = .x$De),
      coords = c("x", "y"), crs = 4326
    )
    
    # # get SDM map as data frame
    dis_sf <- sf::read_sf(paste0("./analyses/sdm/", s, ".shp")) %>% st_buffer(dist = 100000) %>% sf::st_shift_longitude()
    dis_sf <- st_transform(dis_sf, st_crs(marine_regions))

    # Extract values intersecting each polygon
    buffer_distance <- 1000  # 10 km buffer
    ovr_results <- map_dfr(1:nrow(sp_reg_poly), function(i) {
      
      poly <- sp_reg_poly$geometry[i]
      ecoregion_name <- sp_reg_poly$ECOREGION[i]
      
      # Get area of the polygon in km^2
      # Find intersecting points of SDM with the current polygon
      intersecting_points <- dis_sf[st_intersects(dis_sf, poly, sparse = FALSE), ]
      
      # # Create a buffer around each point
      reg_buffers <- st_buffer(intersecting_points, dist = buffer_distance)
      
      # Find intersecting points with the current polygon
      extracted_values <- st_join(reg_buffers, df_sf, join = st_intersects)
      intersecting_points <- df_sf[st_intersects(df_sf, poly, sparse = FALSE), ]
      
      # Calculate the mean "ss" value for intersecting points
      mean_De <- mean(intersecting_points$value, na.rm = TRUE)
      
      # # Calculate the area of each buffer in square kilometers and sum up
      region_area <- sum(as.numeric(st_area(poly) / 1e6))
      
      # Return as a data frame for each polygon processed
      data.frame(
        ecoregion = ecoregion_name,
        popID = substr(pop_count$loc_ID[pop_count$abbrev == sp & pop_count$region == ecoregion_name], 1, 3)[1],
        mean_De = mean_De,
        area_km2 = region_area
      )
    })
    
    bind_rows(ovr_results)
  }) %>%
  ungroup()

# add timepoints based on IBD lengths
N_pop_summary$time <- sapply(N_pop_summary$ibdLength, function(i) {
  l = strsplit(as.character(i), "\\-")[[1]]
  l1 = as.numeric(l[1])
  l2 = as.numeric(l[2])
  ibdL2nG(l1, l2, 5)})

# add guilds
N_pop_summary$guild = factor(sapply(N_pop_summary$species, function(i) 
  data$category[data$abbrev == i][1]), levels = guilds)

# add specialist and generalist as factor
N_pop_summary$behavior <- factor(ifelse(N_pop_summary$guild == "Generalist", "Generalist", "Specialist"),
                                    levels = c("Generalist", "Specialist"))

N_pop_summary$species <- factor(N_pop_summary$species, levels = species)

# duplicate species column for processing
m_dataset_selected$abbrev = as.character(m_dataset_selected$species)

# Process the dataset
m_pop_summary <- m_dataset_selected %>%
  group_by(species, ibdLength) %>%
  group_modify(~ {
    
    sp <- unique(.x$abbrev)
    
    # Retrieve marine_regions for the specific species
    sp_regions <- pop_count$region[pop_count$abbrev == sp]
    
    # Filter for relevant polygons in marine_regions
    if (unique(.x$time) == 19125) {
      sp_reg_poly <- marine_regions_lgm %>%
        filter(ECOREGION %in% sp_regions)
    } else {
      sp_reg_poly <- marine_regions %>%
        filter(ECOREGION %in% sp_regions)
    }
    
    # transform subset MAPS dataset to sf object
    df_sf <- st_as_sf(
      data.frame(x = .x$x, y = .x$y, value = .x$sigma),
      coords = c("x", "y"), crs = 4326
    )

    # Extract values intersecting each polygon
    ovr_results <- map_dfr(1:nrow(sp_reg_poly), function(i) {
      
      poly <- sp_reg_poly$geometry[i]
      ecoregion_name <- sp_reg_poly$ECOREGION[i]
      
      # Find intersecting points with the current polygon
      intersecting_points <- df_sf[st_intersects(df_sf, poly, sparse = FALSE), ]
      
      # Calculate the mean "ss" value for intersecting points
      mean_ss <- mean(intersecting_points$value, na.rm = TRUE)
      
      # # Calculate the area of each buffer in square kilometers and sum up
      region_area <- sum(as.numeric(st_area(poly) / 1e6))
      
      # Return as a data frame for each polygon processed
      data.frame(
        ecoregion = ecoregion_name,
        popID = substr(pop_count$loc_ID[pop_count$abbrev == sp & pop_count$region == ecoregion_name], 1, 3)[1],
        dispersal = mean_ss,
        area_km2 = region_area
      )
    })
    
    bind_rows(ovr_results)
  }) %>%
  ungroup()

# add timepoints based on IBD lengths
m_pop_summary$time <- sapply(m_pop_summary$ibdLength, function(i) {
  l = strsplit(as.character(i), "\\-")[[1]]
  l1 = as.numeric(l[1])
  l2 = as.numeric(l[2])
  ibdL2nG(l1, l2, 5)})

# average duplicated species/pop (if more than one popID)
m_pop_summary <- m_pop_summary %>%
  dplyr::group_by(species, time, popID) %>%
  dplyr::summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

# add guilds
m_pop_summary$guild = factor(sapply(m_pop_summary$species, function(i) data$category[data$abbrev == i][1]),
                                        levels = guilds)

# add specialist and generalist as factor
m_pop_summary$behavior <- factor(ifelse(m_pop_summary$guild == "Generalist", "Generalist", "Specialist"),
                                 levels = c("Generalist", "Specialist"))

m_pop_summary$species <- factor(m_pop_summary$species, levels = species)

## boxplots
N_behavior <- ggplot(N_pop_summary, aes(x = behavior, y = mean_De, color = behavior, fill = behavior)) +
  geom_boxplot(width = 0.5, lwd = 2, outlier.color = NA) +
  geom_boxplot(width = 0.5, fatten = 3, col = "white", outlier.color = NA) +
  geom_jitter(fill = "white", col = "black", shape = 21, size = 3) +
  scale_y_log10(expression(N[c]~"/km"^2), breaks = 10^seq(-5, 3, length.out = 5), 
                labels = parse(text = paste0("10^", seq(-5, 3, length.out = 5)))) +
  scale_x_discrete("") +
  facet_wrap(.~time, labeller = labeller(time = setNames(c("~62.5 ya", "~250 ya", "~19,125 ya"), c(62.5, 250.0, 19125.0)))) + 
  scale_fill_manual(values =  c("#009E73", "#FFD39B")) +
  scale_color_manual(values = c("grey60", "grey60")) +
  theme_linedraw() + 
  coord_cartesian(ylim = c(1e-5, 5e3)) +
  theme(text = element_text(size = 20), 
        plot.margin = ggplot2::margin(0.25, 0.25, 0, 0.25, unit = "in"),
        plot.title = element_text(hjust = 0.5),
        ggh4x.axis.ticks.length.minor = rel(-3/2),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.text.x = element_text(size = 24, face = "italic"),
        strip.background = element_rect(fill = "black"),
        legend.position = "none",
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 18, face = "bold"),
        legend.key = element_rect(fill = NA, color = NA)) +
  stat_compare_means(comparisons = list(c("Generalist", "Specialist")), symnum.args = symnum.args, size = 6)

dis_behavior <- ggplot(m_pop_summary, aes(x = behavior, y = dispersal, color = behavior, fill = behavior)) +
  geom_boxplot(width = 0.5, lwd = 2, outlier.color = NA) +
  geom_boxplot(width = 0.5, fatten = 3, col = "white", outlier.color = NA) +
  geom_jitter(fill = "white", col = "black", shape = 21, size = 3) +
  scale_y_log10(expression("Dispersal Distance (km)")) +
  scale_x_discrete("") +
  facet_wrap(.~time, labeller = labeller(time = setNames(c("~62.5 ya", "~250 ya", "~19,125 ya"), c(62.5, 250.0, 19125.0)))) + 
  scale_fill_manual(values =  c("dodgerblue4", "plum1")) +
  scale_color_manual(values = c("grey60", "grey60")) +
  theme_linedraw() + 
  coord_cartesian(ylim = c(1, 250)) +
  theme(text = element_text(size = 20), 
        plot.margin = ggplot2::margin(0.25, 0.25, 0, 0.35, unit = "in"),
        plot.title = element_text(hjust = 0.5),
        ggh4x.axis.ticks.length.minor = rel(-3/2),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.text.x = element_text(size = 24, face = "italic"),
        strip.background = element_rect(fill = "black"),
        legend.position = "none",
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 18, face = "bold"),
        legend.key = element_rect(fill = NA, color = NA)) +
  stat_compare_means(comparisons = list(c("Generalist", "Specialist")), symnum.args = symnum.args, size = 6)

ggarrange(dis_behavior, N_behavior, nrow = 2, common.legend = TRUE, legend = "none",  hjust = 0, labels = letters[1:2],
          font.label = list(size = 20, face = "bold", color ="black")) + bgcolor("white") + border(color = NA)
ggsave("./Figures/N_dis_behavior_time.png", height = 12, width = 12)


###############################################
# ------------------------------------------- #
# m and N difference between two time periods # 
# ------------------------------------------- #
###############################################
mlimits <- Nlimits <- data.frame()
l1 = "0.02-1" # select first IBD segment length (past)
l2 = "6-Inf"  # select second IBD segment length (recent)

for (prefix in species){
  
  ### population size
  df <- N_dataset_selected %>%
    filter(species == prefix) %>%
    dplyr::select(x, y, De, ibdLength)
  
  df <- spread(df, ibdLength, De)
  
  Ndiff[[prefix]] <- df[, c("x", "y", l1)]
  names(Ndiff[[prefix]])[3] <- "diff_De"
  
  # recent minus ancient
  Ndiff[[prefix]]$diff_De =  df[,l2] - df[,l1]
  
  Nlimits <- rbind(Nlimits, range(Ndiff[[prefix]]$diff_De))
  
  ### migration
  df <- m_dataset_selected%>%
    filter(species == prefix) %>%
    dplyr::select(x, y, sigma, ibdLength)
  
  df <- spread(df, ibdLength, sigma)
  
  mdiff[[prefix]] <- df[, c("x", "y", l1)]
  names(mdiff[[prefix]])[3] <- "diff_sigma"
  
  # recent minus ancient
  mdiff[[prefix]]$diff_sigma =  df[,l2] - df[,l1]
  
  mlimits <- rbind(mlimits, range(mdiff[[prefix]]$diff_sigma))
  
}
colnames(Nlimits) <- colnames(mlimits) <-  c("min", "max")

symlog_transform <- function(x) { ifelse(x == 0, 0, sign(x) * log10(abs(x))) }
symlog_inverse <- function(x) { ifelse(x == 0, 0, sign(x) * 10^abs(x)) }
round_by_magnitude <- function(x) { sign(x) * 10^ceiling(log10(abs(x))) }

m.limits <- symlog_transform(round_by_magnitude(c(min(mlimits$min), max(mlimits$max))))
m.breaks <- seq(m.limits[1], m.limits[2], 1)
m.labels <- format(symlog_inverse(m.breaks), digits = 2, scientific = TRUE)

N.limits <- symlog_transform(round_by_magnitude(c(min(Nlimits$min), max(Nlimits$max))))
N.breaks <- seq(N.limits[1], N.limits[2], 1)
N.labels <- format(symlog_inverse(N.breaks), digits = 2, scientific = TRUE)

# make plots
mdiff_plots <- Ndiff_plots <- list()
for (prefix in species){
  
  sp_path = paste0(path, prefix, "/")
  
  host.category <- data$category[data$abbrev == prefix][1]
  
  ### add features and custom themes to plots
  # get individual coords
  coords <- read.csv(paste0(sp_path, prefix, ".outer"), sep = " ", header = F)
  names(coords) <- c("LON", "LAT")
  
  # get distribution and outer polygon
  dis_sf <- sf::read_sf(paste0("./analyses/sdm/", s, ".shp")) %>% st_buffer(dist = 100000)
  q_sf <- st_as_sf(p$m$layers[[1]]$data, coords = c("x", "y"), crs = crs(dis_sf))
  outer_poly <-  st_difference(q_sf, dis_sf) %>%  st_shift_longitude()
  
  diff_N <- diff_m <- ggplot() + theme_classic() + theme(axis.line = element_blank(), 
                                                         axis.ticks = element_blank(), axis.text.y = element_blank(), 
                                                         axis.text.x = element_blank(), axis.title.x = element_blank(), 
                                                         axis.title.y = element_blank())
  
  
  diff_m <- diff_m + geom_raster(data = mdiff[[prefix]], aes(x = x, y = y, fill = symlog_transform(diff_sigma)), 
                                 alpha = 1) + 
    scale_fill_gradient2(midpoint = 0,
                         limits = m.limits,
                         breaks = m.breaks, 
                         labels = m.labels,
                         name = expression(Delta*sigma["recent-ancient"]), 
                         trans = "identity", na.value = "lightgray") +
    geom_contour(data = mdiff[[prefix]], aes(x = x, y = y, z = symlog_transform(diff_sigma)), 
                 breaks = m.breaks, color = "white") + 
    theme(legend.key.width = unit(0.75, "cm"),
          legend.text = element_text(size = 15),
          legend.key.height = unit(2.25, "cm"),
          legend.title = element_text(size = 15)) + coord_fixed() +
    geom_polygon(data = coords, aes(LON, LAT), fill = NA, col = "grey20", linetype = "dashed", linewidth = 1) +
    scale_x_continuous("", expand = c(0,0)) +
    scale_y_continuous("", expand = c(0,0)) +
    labs(title = gsub("_", " ", data$species[data$abbrev == prefix][1])) +
    geom_sf(data = currents, aes( x = NULL, y = NULL), fill = "white", col = "grey20", linewidth = 0.125, alpha = 0.75) +
    geom_sf(data = ITF_pol, aes( x = NULL, y = NULL), fill = "white", col = "black", linewidth = 0.125, alpha = 0.75) +
    geom_sf(data = world_map, aes( x = NULL, y = NULL), fill = "grey85", size = 1) +
    metR::geom_text_contour(data = mdiff[[prefix]], aes(x = x, y = y, z = diff_sigma), 
                            breaks = symlog_inverse(N.breaks), stroke = 0.2, size = 3, skip = 0, 
                            label.placer = metR::label_placer_n(n = 1,
                                                                rot_adjuster = isoband::angle_halfcircle_bottom()), 
                            check_overlap = TRUE) +
    geom_sf(data = continental_map, aes( x = NULL, y = NULL), fill = "grey85", size = 1) +
    geom_point(data = pop_count[pop_count$abbrev == prefix,], 
               aes(x = lon, y = lat, size = n), shape = 21,  col = "grey20", fill = loc.cols[host.category], stroke = 1) +
    scale_size("# individuals", range = c(0.1,5), limits = c(0.1, max(pop_count$n))) +
    theme_light() +
    coord_sf(xlim = ext[1:2], ylim = ext[3:4]) +  # Replace coord_cartesian with coord_sf
    theme_custom
  
  diff_N <- diff_N + geom_raster(data = Ndiff[[prefix]], aes(x = x, y = y, fill = symlog_transform(diff_De)), 
                                 alpha = 1) + 
    scale_fill_gradient2(midpoint = 0,
                         limits = N.limits,
                         breaks = N.breaks, 
                         labels = N.labels,
                         name = expression(Delta*D[e]["(recent-ancient)"]), 
                         trans = "identity", na.value = "lightgray") +
    geom_contour(data = Ndiff[[prefix]], aes(x = x, y = y, z = symlog_transform(diff_De)), 
                 breaks = N.breaks, color = "white") + 
    theme(legend.key.width = unit(0.75, "cm"),
          legend.text = element_text(size = 15),
          legend.key.height = unit(2.25, "cm"),
          legend.title = element_text(size = 15)) + coord_fixed() +
    geom_sf(data = outer_poly, aes( x = NULL, y = NULL), fill ="grey92", col = "grey92") +
    scale_x_continuous("", expand = c(0,0)) +
    scale_y_continuous("", expand = c(0,0)) +
    labs(title = gsub("_", " ", data$species[data$abbrev == prefix][1])) +
    geom_sf(data = currents, aes( x = NULL, y = NULL), fill = "white", col = "grey20", linewidth = 0.125, alpha = 0.75) +
    geom_sf(data = ITF_pol, aes( x = NULL, y = NULL), fill = "white", col = "black", linewidth = 0.125, alpha = 0.75) +
    geom_sf(data = world_map, aes( x = NULL, y = NULL), fill = "grey85", size = 1) +
    metR::geom_text_contour(data = Ndiff[[prefix]], aes(x = x, y = y, z = diff_De), 
                            breaks = symlog_inverse(N.breaks), stroke = 0.2, size = 3, skip = 0, 
                            label.placer = metR::label_placer_n(n = 1,
                                                                rot_adjuster = isoband::angle_halfcircle_bottom()), 
                            check_overlap = TRUE, rotate = TRUE) +
    geom_sf(data = continental_map, aes( x = NULL, y = NULL), fill = "grey85", size = 1) +
    geom_point(data = pop_count[pop_count$abbrev == prefix,], 
               aes(x = lon, y = lat, size = n), shape = 21,  col = "grey20", fill = loc.cols[host.category], stroke = 1) +
    scale_size("# individuals", range = c(0.1,5), limits = c(0.1, max(pop_count$n))) +
    theme_light() +
    coord_sf(xlim = ext[1:2], ylim = ext[3:4]) +  # Replace coord_cartesian with coord_sf
    theme_custom
  
  mdiff_plots[[prefix]] <- diff_m
  Ndiff_plots[[prefix]] <- diff_N
}

Dm <- ggarrange(mdiff_plots$CLK, mdiff_plots$AKY, mdiff_plots$CRP,
                mdiff_plots$MEL,
                mdiff_plots$AKA, mdiff_plots$PRD,
                mdiff_plots$POL, mdiff_plots$SAN,
                font.label = list(size = 20, face = "bold", color ="black"),
                nrow = 4, ncol = 2, common.legend = TRUE, labels = letters[1:10], legend = "right") +
  bgcolor("white") + border(color = NA)
ggsave(paste0("./Figures/all_species_", ifelse(recombMap, "GenMap_", "NoGenMap_"), ifelse(best_ndemes, "best_ndemes_", ndemes), "demes_delta-mrates_", l2, "-", l1, ".pdf"), 
       plot = Dm, height = 15, width = 20)

DN <- ggarrange(Ndiff_plots$CLK, Ndiff_plots$AKY, Ndiff_plots$CRP,
                Ndiff_plots$MEL,
                Ndiff_plots$AKA, Ndiff_plots$PRD,
                Ndiff_plots$POL, Ndiff_plots$SAN,
                font.label = list(size = 20, face = "bold", color ="black"),
                nrow = 4, ncol = 2, common.legend = TRUE, labels = letters[1:10], legend = "right") +
  bgcolor("white") + border(color = NA)
ggsave(paste0("./Figures/all_species_", ifelse(recombMap, "GenMap_", "NoGenMap_"), ifelse(best_ndemes, "best_ndemes_", ndemes), "demes_delta-Nrates_", l2, "-", l1, ".pdf"), 
       plot = DN, height = 15, width = 20)

