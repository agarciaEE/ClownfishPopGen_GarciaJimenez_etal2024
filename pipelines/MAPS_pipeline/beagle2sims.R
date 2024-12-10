# Install optparse if not installed
if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse")
}

library(optparse)

# Define command line arguments
option_list <- list(
  make_option(c("-i", "--ibdFile"), type = "character", default = NA,
              help = "Path to the IBD file", metavar = "character"),
  make_option(c("-p", "--prefix"), type = "character", default = NA,
              help = "Prefix of the .fam file", metavar = "character"),
  make_option(c("-d", "--dir"), type = "character", default = "./",
              help = "Working directory", metavar = "character"),
  make_option(c("-l", "--lowerBnd"), type = "numeric", default = 2,
              help = "Lower bound value (default is 2)", metavar = "numeric"),
  make_option(c("-u", "--upperBnd"), type = "numeric", default = Inf,
              help = "Upper bound value (default is Inf)", metavar = "numeric")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Parse command line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if required arguments are provided
if (is.na(opt$ibdFile) || is.na(opt$prefix)) {
  print_help(opt_parser)
  stop("Error: The first two arguments (ibdFile and prefix) are required.")
}

# Assign arguments to variables
ibdFile <- opt$ibdFile
prefix <- opt$prefix
path <- opt$dir
lowerBnd <- opt$lowerBnd
upperBnd <- opt$upperBnd

# Output the values for debugging purposes
cat("Script parameters:\n")
cat("ibdFile:", ibdFile, "\n")
cat("prefix:", prefix, "\n")
cat("working directory:", path, "\n")
cat("lowerBnd:", lowerBnd, "\n")
cat("upperBnd:", upperBnd, "\n")

## read in the meta data
meta_info <- read.table(file.path(path, paste0(prefix, ".fam")), sep = " ", as.is = TRUE, header = FALSE, col.names = c("fid", "iid", "PID", "MID", "sex", "phenotype"))
meta_info <- meta_info[,1:2]

ids <- paste0(meta_info$iid, "_", meta_info$fid)
print(ids)

# BEAGLE outputs haplotype data, we treat
# the two different chromosomes of each individual as independent
ids <- c(paste0(ids, "_1"), paste0(ids, "_2"))

n <- length(ids) 

# set up the matrix
ibd_summary <- matrix(nrow = n, ncol = n, 0)
rownames(ibd_summary) <- ids
colnames(ibd_summary) <- ids

## read in the PSC calls. Here, we read in only chromosome 1
## Here, the file was generated using the Snakemake pipeline found here: https://github.com/halasadi/ibd_data_pipeline
## In this pipeline, both the .ibd and .hbd files are utilized (see the relabelIBD rule in the Snakemake file)

ibdFile <- strsplit(ibdFile, ",")[[1]]

multifiles = FALSE
if (length(ibdFile) > 1){
  ibd_data <- lapply(ibdFile, read.table, header = TRUE, stringsAsFactors = FALSE)
  # output filename prefix
  outPrefix = prefix
} else {
  ibd_data <- list(read.table(ibdFile, header = TRUE, stringsAsFactors = FALSE))
  # output filename prefix
  outPrefix = strsplit(basename(ibdFile), "\\.")[[1]][1]
}

for (i in 1:length(ibdFile)){
  # compute the lengths of the PSC segments
  lengths <- ibd_data[[i]]$end-ibd_data[[i]]$start
  
  cat("-", ibdFile[i],"\n")
  cat("Lengths of PSC segments:\n")
  print(lengths)
  
  # Check if all PSC segments are shorter than the selected lower boundary
  if (all(lengths <= lowerBnd)) {
    message("All PSC segments are shorter than the selected lower boundary.")
    #file.create(file.path(path, paste0(outPrefix, ".maps.", lowerBnd, "_", upperBnd, ".sims")), showWarnings = FALSE)
    
    # Check if all PSC segments are longer than the selected upper boundary
  } else if (all(lengths >= upperBnd)) {
    message("All PSC segments are longer than the selected upper boundary.")
    #file.create(file.path(path, paste0(outPrefix, ".maps.", lowerBnd, "_", upperBnd, ".sims")), showWarnings = FALSE)
    
    # Handle cases where segments are within the specified boundaries
  } else {
    # Highlight PSC segments greater than the lower bound and less than the upper bound
    selected_inds <- which(lengths > lowerBnd & lengths < upperBnd)
    
    # Corrected the typo in 'length'
    cat("Number of PSC segments passing the filter:", length(selected_inds), "\n")
    
    # Ensure ibd_summary is initialized before incrementing
    if (!exists("ibd_summary")) {
      ibd_summary <- matrix(0, nrow = max(ids), ncol = max(ids))
      rownames(ibd_summary) <- colnames(ibd_summary) <- as.character(ids)
    }
    
    # if any ibd lenght passed the filter:
    if (length(selected_inds) > 0){
      for (j in 1:length(selected_inds)) {
        id1 <- ibd_data[[i]]$id1[selected_inds[j]]
        id2 <- ibd_data[[i]]$id2[selected_inds[j]]
        if (id1 %in% ids & id2 %in% ids) {
          ibd_summary[id1, id2] <- ibd_summary[id1, id2] + 1
          ibd_summary[id2, id1] <- ibd_summary[id1, id2]
        }
      }
    }
  }
}

# Write similarity matrix to a .sims file
write.table(ibd_summary, file = file.path(path, paste0(outPrefix, ".maps.", lowerBnd, "_", upperBnd, ".sims")),
            quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
