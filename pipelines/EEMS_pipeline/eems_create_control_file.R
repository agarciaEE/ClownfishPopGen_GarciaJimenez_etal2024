#!/usr/bin/Rscript

# Read in the arguments
library("optparse")
option_list = list(
  make_option(c("-i", "--datapath"), type="character", default=NULL, 
              help="Directory where bed files are located (Default: current directory)", metavar="character"),
  make_option(c("-p", "--prefix"), type="character", default="", 
              help="Common prefix in all bed files.", metavar="character"),
  make_option(c("-o", "--overwrite"), type="logical", default=FALSE, 
              help="Logical whether to overwrite results in the existing output directory or create a new MCMC chain (Default: FALSE).", metavar="logical"),
  make_option(c("-N", "--nsamples"), type="integer", default=NA, 
              help="Number of samples (default: NULL, it will be computed)", metavar="integer"),
  make_option(c("-D", "--ndemes"), type="integer", default=200, 
              help="Number of demes (default: 200)", metavar="integer"),
  make_option(c("-S", "--nsites"), type="integer", default=NA, 
              help="Number of sites (default: NULL, it will be computed)", metavar="integer"),
  make_option(c("-C", "--chain"), type="integer", default=1, 
              help="Number of the mcmc chain (Default: 1).", metavar="integer"),
  make_option(c("-P", "--diploid"), type="logical", default=TRUE, 
              help="Logical whether inidividuals are diploid (default: true).", metavar="logical"),
  make_option(c("-M", "--mcmciter"), type="integer", default=2000000, 
              help="Number of mcmc iterations (Default: 2000000)", metavar="integer"),
  make_option(c("-B", "--mcmcburn"), type="integer", default=1000000, 
              help="Burn-in size (Default: 1000000).", metavar="integer"),
  make_option(c("-T", "--mcmcthin"), type="integer", default=999, 
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

cat("\nCreating control file...\n")

files <- list.files(datapath, pattern = prefix)
missing_files <- sapply(c(".coord", ".diffs", ".outer"), function(i) !any(stringr::str_detect(files, i)))

if (any(missing_files)){
  warning("Missing the following input files in datapath folder: ", 
          paste(paste0(prefix, names(missing_files[missing_files])), collapse = ", "))
}

if (is.na(nsamples)){
  if (file.exists(file.path(datapath, paste0(prefix, ".order")))) {
    nsamples <- nrow(read.table(file.path(datapath, paste0(prefix, ".order"))))
  } else if (file.exists(file.path(datapath, paste0(prefix, ".coord")))) {
    nsamples <- nrow(read.table(file.path(datapath, paste0(prefix, ".coord"))))
  } else {
    stop("Missing means to obtain the number of samples.")
  }
}

if (is.na(nsites)){
  cat("\nObtaining number of sites...\n")
  if (file.exists(file.path(datapath, paste0(prefix, ".bim")))) {
    nsites <- as.integer(system(paste("wc -l", file.path(datapath, paste0(prefix, ".bim")), "| awk '{print $1}'"), intern = TRUE))
  } else {
    stop("Missing means to obtain the number of sites")
  }
}

dataPATH <- file.path(datapath, prefix)
outputPATH <- file.path(datapath, paste0(prefix, "-EEMS-ndemes", ndemes, "-chain", chain))

if (dir.exists(outputPATH) & !overwrite){
  warning("Output directory already exists, creating a new mcmc chain...")
  chain = chain + 1
  outputPATH <- file.path(datapath, paste0(prefix, "-EEMS-ndemes", ndemes, "-chain", chain))
}

dir.create(outputPATH, recursive = T, showWarnings = F)

controlFile <- list(datapath = dataPATH,
                    mcmcpath = outputPATH,
                    nIndiv = as.character(format(nsamples, scientific = F)),
                    nSites = as.character(format(nsites, scientific = F)),
                    nDemes = as.character(format(ndemes, scientific = F)),
                    diploid = ifelse(diploid, "true", "false"),
                    numMCMCIter = as.character(format(mcmciter, scientific = F)),
                    numBurnIter = as.character(format(mcmcburn, scientific = F)),
                    numThinIter = as.character(format(mcmcthin, scientific = F))
)

# Convert the list to a dataframe
controlFile <- data.frame(names(controlFile), unlist(controlFile))
# Write the dataframe to a CSV file
write.table(controlFile, quote = F, sep = " = ",  col.names = F, 
            file =  file.path(datapath, paste0(prefix, "_ndemes", ndemes, "_params-chain", chain, ".ini")), row.names = FALSE)
