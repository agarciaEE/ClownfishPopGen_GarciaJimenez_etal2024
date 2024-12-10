###########################
# ----------------------- #
# TEST RECOMBINATION MAPS #
# ----------------------- #
###########################
# working directory
setwd("~/Unil/Research/3_EcoPopGen/project_repository/")

### Custom Functions
source("./scripts/custom_functions.R")


ndemes = 50
prefix = "CLK"
host.category <- "Generalist"
width = 20
height = 12
ibd_lengths <- c("0.02-1", "2-6", "6-Inf")
set.range = TRUE
mlim = c(0,1e3)
Nlim = c(0,1e3)

path <- "./analyses/maps/CLK_testrecombmaps"

RecombMap_list <- list()
RecombMap_plots <-  list()
for (ibd in ibd_lengths){
  

  
  RecombMap_list[[ibd]] <- list()
  RecombMap_plots[[ibd]] <- list()
  
  # get mcmc file paths
  mcmcpaths <- paste0(path, "/", grep("MAPS", list.dirs(path, recursive = FALSE, full.names = FALSE), value = TRUE))
  mcmcpaths <- grep(paste0("ndemes", ndemes, "-", sub("-", "_", ibd), "params-clarkii"), mcmcpaths, value = TRUE)

  RecombMap_plots[[ibd]] <- list()
  for (m in mcmcpaths){
    if (length(list.files(m)) == 1){
      next
    }

    rcm <- strsplit(basename(m), "_")[[1]][3]
    
    # output path indicating the cM range taken
    outpath <- paste0(rcm, "RecombMaps_", ibd, "_", ndemes, "demes")
    
    ## get MAPS plots
    p <- plot_MAPS(add.pts = FALSE, add.graph = FALSE, add.countries = FALSE,
                   longlat = TRUE, mcmcpath = m, correct.by.area = TRUE,
                   set.range = set.range, m.limits = mlim, N.limits = Nlim,
                   outpath = outpath, height = height, width = width, trans = "identity")
    
    RecombMap_list[[ibd]][[rcm]] <- list(m = summary(c(p$m_summary$avg)),
                                         N = summary(c(p$N_summary$avg)))
    RecombMap_plots[[ibd]][[rcm]] <- p 
    
  }
}

RecombMap_mdf <- plyr::ldply(lapply(RecombMap_list, function(i) 
  plyr::ldply(lapply(i, function(r) 
    r[["m"]]), .id = "recombmap")), .id = "ibdsegment")
RecombMap_Ndf <- plyr::ldply(lapply(RecombMap_list, function(i) 
  plyr::ldply(lapply(i, function(r) 
    r[["N"]]), .id = "recombmap")), .id = "ibdsegment")

par(mfrow = c(1,2))
boxplot(Mean ~ ibdsegment, data = RecombMap_Ndf, ylab = "dispersal distance")
boxplot(Mean ~ ibdsegment, data = RecombMap_mdf, ylab = "population density")

dir.create("./Figures", showWarnings = FALSE)

ggarrange(ggarrange(RecombMap_plots$`0.02-1`$`IDN-A`$m, RecombMap_plots$`0.02-1`$`IDN-B`$m, RecombMap_plots$`0.02-1`$`PNG-A`$m, RecombMap_plots$`0.02-1`$`PNG-B`$m,
                    nrow = 1, ncol = 4, common.legend = TRUE, legend = "right"),
          ggarrange(RecombMap_plots$`2-6`$`IDN-A`$m, RecombMap_plots$`2-6`$`IDN-B`$m, RecombMap_plots$`2-6`$`PNG-A`$m, RecombMap_plots$`2-6`$`PNG-B`$m,
                    nrow = 1, ncol = 4, common.legend = TRUE, legend = "right"),
          ggarrange(RecombMap_plots$`6-Inf`$`IDN-A`$m, RecombMap_plots$`6-Inf`$`IDN-B`$m, RecombMap_plots$`6-Inf`$`PNG-A`$m, RecombMap_plots$`6-Inf`$`PNG-B`$m,
                    nrow = 1, ncol = 4, common.legend = TRUE, legend = "right"),
          font.label = list(size = 20, face = "bold", color ="black"),
          nrow = 3, ncol = 1, common.legend = FALSE, labels = letters[1:3]) +
  bgcolor("white") + border(color = NA)
ggsave(paste0("./Figures/clarkii_RecombMap_test_sigma.png"),  height = 15, width = 20)

ggarrange(ggarrange(RecombMap_plots$`0.02-1`$`IDN-A`$N, RecombMap_plots$`0.02-1`$`IDN-B`$N, RecombMap_plots$`0.02-1`$`PNG-A`$N, RecombMap_plots$`0.02-1`$`PNG-B`$N,
                    nrow = 1, ncol = 4, common.legend = TRUE, legend = "right"),
          ggarrange(RecombMap_plots$`2-6`$`IDN-A`$N, RecombMap_plots$`2-6`$`IDN-B`$N, RecombMap_plots$`2-6`$`PNG-A`$N, RecombMap_plots$`2-6`$`PNG-B`$N,
                    nrow = 1, ncol = 4, common.legend = TRUE, legend = "right"),
          ggarrange(RecombMap_plots$`6-Inf`$`IDN-A`$N, RecombMap_plots$`6-Inf`$`IDN-B`$N, RecombMap_plots$`6-Inf`$`PNG-A`$N, RecombMap_plots$`6-Inf`$`PNG-B`$N,
                    nrow = 1, ncol = 4, common.legend = TRUE, legend = "right"),
          font.label = list(size = 20, face = "bold", color ="black"),
          nrow = 3, ncol = 1, common.legend = FALSE, labels = letters[1:3]) +
  bgcolor("white") + border(color = NA)
ggsave(paste0("./Figures/clarkii_RecombMap_test_De.png"),  height = 15, width = 20)

ggarrange(ggarrange(RecombMap_plots$`0.02-1`$`IDN-A`$m_sign, RecombMap_plots$`0.02-1`$`IDN-B`$m_sign, RecombMap_plots$`0.02-1`$`PNG-A`$m_sign, RecombMap_plots$`0.02-1`$`PNG-B`$m_sign,
                    nrow = 1, ncol = 4, common.legend = TRUE, legend = "right"),
          ggarrange(RecombMap_plots$`2-6`$`IDN-A`$m_sign, RecombMap_plots$`2-6`$`IDN-B`$m_sign, RecombMap_plots$`2-6`$`PNG-A`$m_sign, RecombMap_plots$`2-6`$`PNG-B`$m_sign,
                    nrow = 1, ncol = 4, common.legend = TRUE, legend = "right"),
          ggarrange(RecombMap_plots$`6-Inf`$`IDN-A`$m_sign, RecombMap_plots$`6-Inf`$`IDN-B`$m_sign, RecombMap_plots$`6-Inf`$`PNG-A`$m_sign, RecombMap_plots$`6-Inf`$`PNG-B`$m_sign,
                    nrow = 1, ncol = 4, common.legend = TRUE, legend = "right"),
          font.label = list(size = 20, face = "bold", color ="black"),
          nrow = 3, ncol = 1, common.legend = FALSE, labels = letters[1:3]) +
  bgcolor("white") + border(color = NA)
ggsave(paste0("./Figures/clarkii_RecombMap_test_Psigma.png"),  height = 15, width = 20)

ggarrange(ggarrange(RecombMap_plots$`0.02-1`$`IDN-A`$N_sign, RecombMap_plots$`0.02-1`$`IDN-B`$N_sign, RecombMap_plots$`0.02-1`$`PNG-A`$N_sign, RecombMap_plots$`0.02-1`$`PNG-B`$N_sign,
                    nrow = 1, ncol = 4, common.legend = TRUE, legend = "right"),
          ggarrange(RecombMap_plots$`2-6`$`IDN-A`$N_sign, RecombMap_plots$`2-6`$`IDN-B`$N_sign, RecombMap_plots$`2-6`$`PNG-A`$N_sign, RecombMap_plots$`2-6`$`PNG-B`$N_sign,
                    nrow = 1, ncol = 4, common.legend = TRUE, legend = "right"),
          ggarrange(RecombMap_plots$`6-Inf`$`IDN-A`$N_sign, RecombMap_plots$`6-Inf`$`IDN-B`$N_sign, RecombMap_plots$`6-Inf`$`PNG-A`$N_sign, RecombMap_plots$`6-Inf`$`PNG-B`$N_sign,
                    nrow = 1, ncol = 4, common.legend = TRUE, legend = "right"),
          font.label = list(size = 20, face = "bold", color ="black"),
          nrow = 3, ncol = 1, common.legend = FALSE, labels = letters[1:3]) +
  bgcolor("white") + border(color = NA)
ggsave(paste0("./Figures/clarkii_RecombMap_test_PDe.png"), height = 15, width = 20)
