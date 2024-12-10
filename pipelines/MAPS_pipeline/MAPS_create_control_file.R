#!/usr/bin/Rscript
# example run: Rscript MAPS_create_control_file.R -e /Users/agarciaj/Unil/Research/3_EcoPopGen/Results/Aclarkii_ref/EEMS/CLK_eemsFiles/ -i /Users/agarciaj/Unil/Research/3_EcoPopGen/Results/Aclarkii_ref/MAPS/CLK_MAPSFiles -p CLK -o FALSE -n NULL -d 500 -g 1000 -l 0.1 -u Inf -m 2000000 -b 1000000 -t 1000 

# Read in the arguments
library("optparse")
option_list = list(
  make_option(c("-e", "--eemspath"), type="character", default=NULL, 
              help="Directory where EEMS files are located (Default: current directory)", metavar="character"),
  make_option(c("-i", "--datapath"), type="character", default=NULL, 
              help="Directory where EEMS files are located (Default: current directory)", metavar="character"),
  make_option(c("-p", "--prefix"), type="character", default="", 
              help="Common prefix in all bed files.", metavar="character"),
  make_option(c("-o", "--overwrite"), type="logical", default=FALSE, 
              help="Logical whether to overwrite results in the existing output directory or create a new MCMC chain (Default: FALSE).", metavar="logical"),
  make_option(c("-n", "--nsamples"), type="integer", default=-1, 
              help="Number of samples (default: -1, it will be computed)", metavar="integer"),
  make_option(c("-d", "--ndemes"), type="integer", default=200, 
              help="Number of demes (default: 200)", metavar="integer"),
  make_option(c("-g", "--genomesize"), type="integer", default=1000, 
              help="Approximate genome length (default: 1000 cM)", metavar="integer"),
  make_option(c("-c", "--chain"), type="integer", default=1, 
              help="Number of the mcmc chain (Default: 1).", metavar="integer"),
  make_option(c("-l", "--lower"), type="numeric", default=0.1, 
              help="Lower PSC segment boundary (default: 0.1).", metavar="numeric"),
  make_option(c("-u", "--upper"), type="numeric", default=Inf, 
              help="Upper PSC segment boundary (default: Inf).", metavar="numeric"),
  make_option(c("-m", "--mcmciter"), type="integer", default=2000000, 
              help="Number of mcmc iterations (Default: 2000000)", metavar="integer"),
  make_option(c("-b", "--mcmcburn"), type="integer", default=1000000, 
              help="Burn-in size (Default: 1000000).", metavar="integer"),
  make_option(c("-t", "--mcmcthin"), type="integer", default=999, 
              help="Thining size", metavar="integer")
) 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Check that all required arguments are provided
if (is.null(opt[["datapath"]])){
  print_help(opt_parser)
  stop("Please provide the data path", call.=FALSE)
}

invisible(sapply(1:length(opt), function(i) assign(names(opt)[i], opt[[i]], envir = globalenv())))

system(paste("mkdir -p", datapath))

cat("\nCreating control file...\n")

inputFilesExt <- c(".coord", ".sims", ".outer", ".ipmap", paste0("_", ndemes,".edges"), paste0("_", ndemes,".demes"))

files <- list.files(datapath, pattern = prefix)
missing_files <- sapply(inputFilesExt, function(i) !any(stringr::str_detect(files, i)))
missing_inputFilesExt <- names(missing_files[missing_files])

if (any(missing_files)){
  warning("Missing the following input files in datapath folder: ", 
          paste(paste0(prefix, missing_inputFilesExt, collapse = ", ")))
  if (!is.null(eemspath)){
    cat("Copying from eems data folder...")
    if (".coord" %in% missing_inputFilesExt) {
      coord <- read.table(file.path(eemspath, paste0(prefix, ".coord")), col.names = c("x", "y"))
      locs <- paste0(coord$x, " ", coord$y)
      locs <- c(locs, locs)
      write.table(data.frame(locs), file = file.path(datapath, paste0(prefix, ".coord")), 
                  quote=FALSE, sep = " ", row.names = FALSE, col.names=FALSE)
    }
    if (".ipmap" %in% missing_inputFilesExt) {
      ipmap <- read.table(file.path(eemspath, paste0(prefix, "-EEMS-ndemes", ndemes, "-chain1"), "ipmap.txt"), col.names = c("deme"))
      ipmap <- c(ipmap$deme, ipmap$deme)
      write.table(data.frame(ipmap), file = file.path(datapath, paste0(prefix,  "_", ndemes, ".ipmap")), 
                  quote=FALSE, sep = " ", row.names = FALSE, col.names=FALSE)
    }
    if (any(paste0("_", ndemes, c(".edges", ".demes")) %in% missing_inputFilesExt)) {
      for (i in c("edges", "demes")){
        system(paste("cp", paste0(file.path(eemspath, paste0(prefix, "-EEMS-ndemes", ndemes, "-chain1")), "/", paste0(i, ".txt"), collapse = " "), file.path(datapath, paste0(prefix, "_", ndemes, ".", i))))
      }
    }
    if (".outer" %in% missing_inputFilesExt) {
      system(paste("cp", paste0(file.path(eemspath, paste0(prefix, "-EEMS-ndemes", ndemes, "-chain1")), "/", paste0("outer.txt"), collapse = " "), file.path(datapath, paste0(prefix, ".outer"))))
    }
    else {
      message("Please provide missing files to the input data folder.")
    }
  } else {
    message("Please provide missing files to the input data folder.")
  }
}

dataPATH <- file.path(datapath, prefix)
outputPATH <- file.path(datapath, paste0(prefix, "-MAPS-ndemes", ndemes, "-chain", chain))

if (dir.exists(outputPATH) & !overwrite){
  warning("Output directory already exists, creating a new mcmc chain...")
  chain = chain + 1
  outputPATH <- file.path(datapath, paste0(prefix, "-MAPS-ndemes", ndemes, "-chain", chain))
}

dir.create(outputPATH, recursive = T, showWarnings = F)

if (nsamples == -1){
  cat("Counting number of individuals...")
  if (file.exists(file.path(datapath, paste0(prefix, ".ipmap")))) {
    nsamples <- nrow(read.table(file.path(datapath, paste0(prefix, ".ipmap"))))
  } else if (file.exists(file.path(datapath, paste0(prefix, ".coord")))) {
    nsamples <- nrow(read.table(file.path(datapath, paste0(prefix, ".coord"))))
  } else if (file.exists(file.path(datapath, paste0(prefix, ".sims")))) {
    nsamples <- nrow(read.table(file.path(datapath, paste0(prefix, ".sims"))))
  } else {
    stop("Missing means to obtain the number of samples.")
  }
}

controlFile <- list(datapath = dataPATH,
                    gridpath = dataPATH,
                    mcmcpath = outputPATH,
                    nIndiv = as.character(format(nsamples, scientific = F)),
                    nDemes = as.character(format(ndemes, scientific = F)),
                    genomeSize = as.character(format(genomesize, scientific = F)),
                    numMCMCIter = as.character(format(mcmciter, scientific = F)),
                    numBurnIter = as.character(format(mcmcburn, scientific = F)),
                    numThinIter = as.character(format(mcmcthin, scientific = F)),
                    lowerBound =  as.character(format(lower, scientific = F)),
                    upperBound =  as.character(format(upper, scientific = F))
)

# Convert the list to a dataframe
controlFile <- data.frame(names(controlFile), unlist(controlFile))
# Write the dataframe to a CSV file
write.table(controlFile, quote = F, sep = " = ",  col.names = F, 
            file =  file.path(datapath, paste0(prefix, "_ndemes", ndemes, "_params-chain", chain, ".ini")), row.names = FALSE)

