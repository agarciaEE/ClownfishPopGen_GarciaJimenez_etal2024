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

#####################
### ------------- ###
#####################

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

# set factors
data$category <- factor(data$category, levels = guilds)
data$abbrev <- factor(data$category, levels = species)
data$country_ID <- factor(data$category, levels = popLevels)

# generation time
gen = 5
# mutation rate
mu = 4E-08

#####################
### ------------- ###
#####################

# load msmc reconstructions (10 individual-based boostrap reconstructions per species/population)
##################################################################################################
msmc10_path = "./analyses/msmc2/ind_bootstraps/"
msmc10_files <- list.files(msmc10_path, pattern = "final.txt$")

# read and format msmc data
msmc10_data <- lapply(msmc10_files, function(i) read.table(paste0(msmc10_path, i), sep = "\t", header = T))
names(msmc10_data) <- sub("\\.msmc2.final.txt", "", sub("all\\.", "", msmc10_files))
msmc10_data <- plyr::ldply(msmc10_data, .id = "id")

# add samples info
msmc10_data$species <- substr(msmc10_data$id, 1,3)
msmc10_data$pop <- sub("\\.[a-z][0-9]*", "", sub("[A-Z]*_", "", msmc10_data$id))
msmc10_data$bootstrap <- sub("[A-Z]*_[A-Z]*\\.", "", msmc10_data$id)

# compute time (years) and Ne
msmc10_data$left_time_boundary <- gen * msmc10_data$left_time_boundary / mu
msmc10_data$right_time_boundary <- gen * msmc10_data$right_time_boundary / mu
msmc10_data$Ne <- 1 / (msmc10_data$lambda * (2 * mu)) / 1e4

## Correct and format dataset 
#############################
# remove NA populations
msmc10_data <- msmc10_data[msmc10_data$pop != "NA",]
msmc10_data <- msmc10_data[!is.na(msmc10_data$pop),]

# correct Maldives population ID
msmc10_data$pop[msmc10_data$pop == "MVD"] = "MDV"

# add host categories (guilds)
msmc10_data$guild <- factor(sapply(msmc10_data$species, function(i) data$category[data$abbrev == i][1]), levels = levels(data$category))
# add scientific name
msmc10_data$sci_name <- sub("_", " ", sapply(msmc10_data$species, function(i) data$species[data$abbrev == i][1]))

# get middle time point between time boundaries
msmc10_data$time <- (msmc10_data$right_time_boundary + msmc10_data$left_time_boundary)/2
msmc10_data$time[!is.finite(msmc10_data$time)] = msmc10_data$left_time_boundary[!is.finite(msmc10_data$time)] # assign left time boundary as time when right time boundary is Inf

# order dataset
msmc10_data <- msmc10_data[,c("id", "guild", "sci_name", "species", "pop", "bootstrap", "time_index", "left_time_boundary", "time", "right_time_boundary", "lambda", "Ne")]

# write msmc data
write.csv(msmc10_data, file = "./analyses/msmc2/ind_bootstrap_msmc_data_raw.csv")

# load msmc reconstructions (50 genomic-window-based boostrap reconstructions per species/population)
#####################################################################################################
msmc50_path <- "./analyses/msmc2/multihet_bootstraps/"
msmc50_files <- sapply(species, function(sp) list.files(file.path(msmc50_path, sp), pattern = "final.txt$"))

# read and format msmc data
msmc50_data <- do.call(rbind, lapply(species, function(sp) { 
                    do.call(rbind, # combine msmc bootstraps
                            lapply(msmc50_files[[sp]], function(i) {
                      df <- read.table(file.path(msmc50_path, sp, i), sep = "\t", header = T) # read file
                      df$species = sp # add species name
                      df$pop = strsplit(strsplit(i, "_")[[1]][2], "\\.")[[1]][1]  # add pop id
                      df$bootstrap = strsplit(strsplit(i, "_")[[1]][3], "\\.")[[1]][1] # add bootstrap number
                      df$id = paste0(df$species, "_", df$pop, ".", df$bootstrap)  # add id with species_pop.bootstrap code
                      df
                    }))
                }))

# compute time (years) and Ne
msmc50_data$left_time_boundary <- gen * msmc50_data$left_time_boundary / mu
msmc50_data$right_time_boundary <- gen * msmc50_data$right_time_boundary / mu
msmc50_data$Ne <- (1 / msmc50_data$lambda) / (2 * mu) / 1e4

## Correct and format dataset 
#############################
# remove NA populations
msmc50_data <- msmc50_data[msmc50_data$pop != "NA",]
msmc50_data <- msmc50_data[!is.na(msmc50_data$pop),]

# correct Maldives population ID
msmc50_data$pop[msmc50_data$pop == "MVD"] = "MDV"

# add host categories (guilds)
msmc50_data$guild <- factor(sapply(msmc50_data$species, function(i) data$category[data$abbrev == i][1]), levels = levels(data$category))

# add scientific name
msmc50_data$sci_name <- sub("_", " ", sapply(msmc50_data$species, function(i) data$species[data$abbrev == i][1]))
  
# get middle time point between time boundaries
msmc50_data$time <- (msmc50_data$right_time_boundary + msmc50_data$left_time_boundary)/2
msmc50_data$time[!is.finite(msmc50_data$time)] = msmc50_data$left_time_boundary[!is.finite(msmc50_data$time)] # assign left time boundary as time when right time boundary is Inf

# order dataset
msmc50_data <- msmc50_data[,c("id", "guild", "sci_name", "species", "pop", "bootstrap", "time_index", "left_time_boundary", "time", "right_time_boundary", "lambda", "Ne")]

# write msmc data
write.csv(msmc50_data, file = "./analyses/msmc2/multihet_bootstrap_msmc_data_raw.csv")

#####################
### ------------- ###
#####################

# x-axis 
exps <- c(0,7)
xlims <- 10 ^ exps
xlabs <- sapply(exps[1]:exps[2], function(i) 10^i)
xticks <- unlist(lapply(2:length(xlabs), function(i) seq(xlabs[i-1], xlabs[i], length.out = 10)))
xtickslabs <- sapply(xticks, function(i) ifelse(i %in% xlabs, format(i, scientific = TRUE), ""))
xlab <- bquote(paste("Years (g=", .(gen), ";", ~ mu, "=", .(format(mu, scientific = TRUE)), ")"))

# y-axis
ylims <- c(0,10) 
ylab <- expression("Effective population size" ~(x10^4))

# visualize data
ggplot(msmc10_data, aes(x = time, y = Ne)) +
  geom_step(aes(group = interaction(species, pop, bootstrap), col = guild), linewidth = 1) +
  scale_x_log10(xlab,  breaks = xticks, labels = xtickslabs, expand = c(0,0)) +
  scale_y_continuous(ylab, expand = c(0,0), minor_breaks = seq(ylims[1], ylims[2], 0.25), guide = "axis_minor") + # this is added to the original code) 
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

# Some SD specialist bootstraps (specifically A. Sandaracinos) show really weird and unrealistic Ne in midpoints (potential signal of ancient hybridization??)
SAN_unrealistic_bootstraps <- unique(msmc10_data$id[msmc10_data$species == "SAN" & msmc10_data$time_index %in% 2:18 & msmc10_data$Ne > Ne10_threshold])

# A. sandaracions Australian population completely removed (Hybrid population)
# A. sandaracions Papua New Guinea population 9 out of 10 individual-based bootstraps removed (??)

msmc10_data_filtered <- msmc10_data[!msmc10_data$id %in% SAN_unrealistic_bootstraps,]

## Some time edges with unrealistic Ne values

# Define Ne threshold 
Ne10_threshold <- quantile(msmc10_data$Ne, 0.95)

# Filter based on Ne threshold to remove unrealistic estimates at the time edges
msmc10_data_filtered <- msmc10_data_filtered[msmc10_data_filtered$Ne <= Ne10_threshold,] 

# Same for msmsc50_data
SAN_unrealistic_bootstraps <- unique(msmc50_data$id[msmc50_data$species == "SAN" & msmc50_data$time_index %in% 2:18 & msmc50_data$Ne > Ne10_threshold])
msmc50_data_filtered <- msmc50_data[!msmc50_data$id %in% SAN_unrealistic_bootstraps,]
msmc50_data_filtered <- msmc50_data_filtered[msmc50_data_filtered$Ne <= Ne10_threshold,] 

# visualize data
ggplot(msmc10_data_filtered, aes(x = time, y = Ne)) +
  geom_step(aes(group = interaction(species, pop, bootstrap), col = guild), linewidth = 1) +
  scale_x_log10(xlab,  breaks = xticks, labels = xtickslabs, expand = c(0,0)) +
  scale_y_continuous(ylab, expand = c(0,0), minor_breaks = seq(ylims[1], ylims[2], 0.25), guide = "axis_minor") + # this is added to the original code) 
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

# write filtered msmc data
write.csv(msmc10_data_filtered, file = "./analyses/msmc2/ind_bootstrap_msmc_data_filtered.csv")
write.csv(msmc50_data_filtered, file = "./analyses/msmc2/multihet_bootstrap_msmc_data_filtered.csv")

