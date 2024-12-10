####################
# load tidyverse package
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(ggnewscale)
library(ggpubr)
library(ggtree)
library(vegan)

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

# remove A. polymnus samples that appear to cluster with another species (potential mislabelling or misidentification)
misidentified_samples <- c('GB_046', 'GB_048', 'GB_049')

pca_data <- data.frame()
admix_data <- data.frame()
for (prefix in species){
  
  eigenvecFile=paste0("./analyses/pca/", prefix, ".eigenvec")
  eigenvalFile=paste0("./analyses/pca/", prefix, ".eigenval")
  
  popInfo <- read.table(paste0("./raw_data/", prefix, "_popList.txt"), col.names = c("ind", "pop"))
  popInfo$pop[is.na(popInfo$pop)] = "NotDefined"
  popInfo$pop[popInfo$pop == "MVD"] = "MDV"
  
  # read in data
  pca <- read_table(eigenvecFile, col_names = F)
  eigenval <- scan(eigenvalFile)

  # sort out the pca data
  # remove redundant column
  pca <- pca[,-1]
  # set names
  names(pca)[1] <- "ind"
  names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
  
  # remove misidentified or mislabeled samples
  pca <- pca[!pca$ind %in% misidentified_samples,]
  
  nind = nrow(pca)
  cat(prefix, ":", nind)
  if (nind > 100) { 
    width = 20 
    text_size = 50
    point_size = 12.5
    blwd = 4.5
    angle = 45
  } else if (nind < 8) { 
    width = 17 
    text_size = 25
    point_size = 8.5
    blwd = 2.5
    angle = 30
  } else {
    width = 20
    text_size = 40
    point_size = 10.5
    tree.size = 2
    blwd = 2.5
    angle = 30
  }
  
  pop <- popInfo$pop[match(pca$ind, popInfo$ind)]
  pop[is.na(pop)] = "NotDefined"
  
  # remake data.frame
  pca <- as.tibble(data.frame(pca, pop))

  pop_levels <- popLevels[popLevels %in% pca$pop]
  
  # get populations lon lat coordinates
  pop_coordinates <- data.frame(
    pop = pop_levels,
    x = sapply(pop_levels, function(i) mean(data$lon[data$abbrev == prefix & data$country_ID == i], na.rm = TRUE)),
    y = sapply(pop_levels, function(i) mean(data$lat[data$abbrev == prefix & data$country_ID == i], na.rm = TRUE))
  )
  pop_coordinates <- na.exclude(pop_coordinates)
  pop_coordinates$x[pop_coordinates$x < 0] <- pop_coordinates$x[pop_coordinates$x < 0] + 360 # shift longitudes from -180,180 to 0,360
  
  # add country names
  pop_coordinates$country <- popnames[pop_coordinates$pop]
  
  pop_levels <- pop_coordinates$pop[order(pop_coordinates$x)]
  
  pca <- pca[pca$pop %in% pop_levels,]
  pca$pop <- factor(pca$pop, levels = pop_levels)
 
  popcol <- setNames(viridis::cividis(length(levels(pca$pop))), levels(pca$pop))
  guildcol <- guildcols[data$category[data$abbrev == prefix][1]]
  
  pve <- eigenval^2/sum(eigenval^2) * 100
  # plot pca
  pca_aspect_ratio <- diff(range(pca$PC2))/diff(range(pca$PC1))
  if(pca_aspect_ratio < 0.5) {pca_aspect_ratio = 0.5}
  if(pca_aspect_ratio > 2) {pca_aspect_ratio = 2}
  pca_aspect_ratio <- pca_aspect_ratio + sign(1 - pca_aspect_ratio) * 0.4

  pca.plt <- ggplot(pca, aes(PC1, PC2, fill = pop)) + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey20") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey20") +
    geom_point(size = point_size + 5, shape = 21) + 
    scale_fill_manual("Population", values = popcol) +
    labs(col = "Population") +
    theme_light() + 
    xlab(paste0("PC1 (", signif(pve[1], 3), "%)")) + 
    ylab(paste0("PC2 (", signif(pve[2], 3), "%)")) +
    theme(text = element_text(size = text_size + 14, color = "black"), 
          plot.background = element_rect(fill = "white", color = "transparent"), 
          plot.margin = margin(10, 15, 10, 15),
          panel.border = element_rect(linewidth = blwd + 2, color = "black"),
          panel.grid = element_blank(),
          legend.position = "right")

  ## Admixture plot
  admixCVerrorFile=paste0("./analyses/admixture/", prefix, "_CVerror.txt")
  CVerr <- read_table(admixCVerrorFile, col_names = c("k", "CVerr"))
  K <- CVerr$k[which.min(CVerr$CVerr)]
  
  admixFile=paste0("./analyses/admixture/", prefix, ".", K, ".Q")
  admixTbl <- read_table(admixFile, col_names = paste0("K", 1:K))
  
  # remake data.frame
  admixdf <- as.tibble(data.frame(popInfo, admixTbl))
  admixdf$pop <- factor(admixdf$pop, levels = pop_levels)
  admixdf <- admixdf[admixdf$pop %in% pop_levels,] # remove NA populations
  admixdf <- admixdf[!admixdf$ind %in% misidentified_samples,] #remove mislabeled or misidentified samples
  
  admixdf <- admixdf[order(admixdf$pop),]
  admixdf$ind <- factor(admixdf$ind, levels = unique(admixdf$ind))
  
  admixdfg <- gather(admixdf, K, prop, 3:ncol(admixdf))
  
  admxplot <- ggplot(admixdfg, aes(x = ind, y = prop, fill = K)) + 
                    geom_bar(stat = "identity", width = 0.9) +
                    scale_fill_manual(values = colorRampPalette(c(ifelse(guildcol == "black", "gray60" , colorspace::lighten(guildcol, 0.25)), 
                                                                  ifelse(guildcol == "black", "gray20", colorspace::darken(guildcol, 0.25))))(n = K)) +
                    scale_y_continuous(paste("K =", K), expand = c(0,0), breaks = seq(0,1,0.2), limits = c(0,1.1)) +
                    theme_bw() +
                    theme(legend.position = "none",
                          axis.title.x = element_blank(),
                          panel.grid = element_blank(), 
                          panel.border = element_blank(), 
                          plot.background = element_rect(fill = "white", color = "transparent"), 
                          plot.margin = margin(12, 15, 10, 15),
                          axis.ticks = element_blank(),
                          axis.title.y = element_text(size = text_size + 10, face = "bold", vjust = 2),
                          axis.text.y = element_text(size = text_size+12 ),
                          axis.text.x = element_text(size = ifelse(nind>100,30, text_size-10), angle = 45, hjust = 1, vjust = 1)) 
  
  for (p in levels(admixdf$pop)){
    idx <- range(which(admixdf$pop == p))
    admxplot <- admxplot + annotate("rect", xmin = idx[1]-0.45, xmax = idx[2]+0.45, ymin = 1.01, ymax = 1.05,
                                     fill= popcol[p])
  }

  admix_data <- rbind(admix_data, cbind(species = prefix, admixdfg))
  pca_data <- rbind(pca_data, cbind(species = prefix, pca[,c("ind", "PC1", "PC2", "pop")]))
  
  ggarrange(admxplot, pca.plt, nrow = ifelse(pca_aspect_ratio < 1, 2, 1), ncol = ifelse(pca_aspect_ratio < 1, 1, 2),
            common.legend = FALSE) +
    bgcolor("white") +
    border(color = "white") +
    theme(plot.background = element_rect(fill = "white", color = "white"),
          plot.margin = margin(20, 20, 20, 20))

  ggsave(paste0("./Figures/", prefix, "_pca_admx.png"), width = width, height = height, limitsize = FALSE)
}

# write pca and admix data
write.csv(pca_data, file = "./analyses/pca/allspecies_pca_data.csv")
write.csv(admix_data, file = "./analyses/admixture/allspecies_admixture_data.csv")
      

