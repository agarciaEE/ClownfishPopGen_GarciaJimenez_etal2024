#!/usr/bin/Rscript
# example run:  Rscript Rscripts/eems_generate_input.R --infofile data/popgen_dataset.csv --beddir eems --bedprefix _AperRef_Filter1 --outputdir eems --overwrite TRUE

# Read in the arguments
library("optparse")
option_list = list(
  make_option(c("-f", "--infofile"), type="character", default=NULL, 
              help="File name of the dataset with all the information needed (coordinates, pop ID, sample ID.", metavar="character"),
  make_option(c("-p", "--bedprefix"), type="character", default="", 
              help="Common prefix in all bed files.", metavar="character"),
  make_option(c("-d", "--beddir"), type="character", default=".", 
              help="Directory where bed files are located (Default: current directory).", metavar="character"),
  make_option(c("-o", "--outputdir"), type="character", default=".", 
              help="Directory to output files (Default: current directory).", metavar="character"),
  make_option(c("-N", "--ndemes"), type="integer", default=200, 
              help="Number of demes (default: 200)", metavar="integer"),
  make_option(c("-C", "--chain"), type="integer", default=1, 
              help="Number of the mcmc chain (Default: 1).", metavar="integer"),
  make_option(c("-D", "--diploid"), type="logical", default=TRUE, 
              help="Logical whether inidividuals are diploid (default: true).", metavar="logical"),
  make_option(c("-M", "--mcmciter"), type="integer", default=2000000, 
              help="Number of mcmc iterations (Default: 2000000)", metavar="integer"),
  make_option(c("-B", "--mcmcburn"), type="integer", default=1000000, 
              help="Burn-in size (Default: 1000000).", metavar="integer"),
  make_option(c("-T", "--mcmcthin"), type="integer", default=999, 
              help="Thining size", metavar="integer"),
  make_option(c("-I", "--indcol"), type="integer", default=1, 
              help="Column number containing individual ID information (Default: 1).", metavar="integer"),
  make_option(c("-S", "--spscol"), type="integer", default=2, 
              help="Column number containing species ID information (Default: 2).", metavar="integer"),
  make_option(c("-A", "--abrcol"), type="integer", default=3, 
              help="Column number containing species abbreviation information (Default: 2).", metavar="integer"),
  make_option(c("-P", "--popcol"), type="integer", default=9, 
              help="Column number containing population ID information (Default: 9).", metavar="integer"),
  make_option(c("-X", "--loncol"), type="integer", default=10, 
              help="Column number containing latitud information (Default: 10).", metavar="integer"),
  make_option(c("-Y", "--latcol"), type="integer", default=11, 
              help="Column number containing longitud information (Default: 11).", metavar="integer"),
  make_option(c("-t", "--translon"), type="logical", default=TRUE, 
              help="Logical whether translocate negative longitud coordinates to ceter on the Indo-Pacific (Default: TRUE).", metavar="logical"),
  make_option(c("-O", "--overwrite"), type="logical", default=FALSE, 
              help="Logical whether to overwrite results in the existing output directory or create a new MCMC chain (Default: FALSE).", metavar="logical"),
  make_option(c("-r", "--regshfile"), type="character", default="/Users/agarciaj/Unil/Research/3_EcoPopGen/MEOW_ECOS/meow_ecos.shp", 
              help="file path to marine regions shape file (Default: /Users/agarciaj/Unil/Research/3_EcoPopGen/MEOW_ECOS/meow_ecos.shp).", metavar="character"),
  make_option(c("-s", "--distshdir"), type="character", default="/Users/agarciaj/Unil/Research/3_EcoPopGen/ENMs_shapefiles/", 
              help="file path to directory containing species distribution shape files (Default: /Users/agarciaj/Unil/Research/3_EcoPopGen/ENMs_shapefiles/).", metavar="character"),
  make_option(c("-k", "--plinkpath"), type="character", default="/work/FAC/FBM/DBC/nsalamin/software/nsalamin/clownfishes/agarciaj/plink", 
              help="file path to plink software (Default: /Users/agarciaj/Softwares/plink/plink).", metavar="character"),
  make_option(c("-b", "--bed2diffspath"), type="character", default="/work/FAC/FBM/DBC/nsalamin/software/nsalamin/clownfishes/agarciaj/eems/bed2diffs/src-wout-openmp/bed2diffs_v1", 
              help="File path to bed2diffs software (Default: /Users/agarciaj/Softwares/eems/bed2diffs/src-wout-openmp/bed2diffs_v1)", metavar="character"),
  make_option(c("-g", "--makeouterpath"), type="character", default="/work/FAC/FBM/DBC/nsalamin/clownfish/agarciaj/EEMS_pipeline/eems_generate_outerfile.R", 
              help="File path to generate_outerfile.R script (Default: /Users/agarciaj/Unil/Research/3_EcoPopGen/Rscripts/generate_outerfile.R)", metavar="character")
  ) 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
cat(names(opt))
# Check that all required arguments are provided
if (is.null(opt[["infofile"]])){
  print_help(opt_parser)
  stop("Please provide the info file", call.=FALSE)
}

invisible(sapply(1:length(opt), function(i) assign(names(opt)[i], opt[[i]], envir = globalenv())))

dir.create(outputdir, recursive = T, showWarnings = F)
data <- read.csv(infofile)

for (i in unique(data[,spscol])) {

  abrev <- data[data[,spscol] == i, abrcol][1]
  cat("\nPreparing files for", i, "...")
  outputdir_i <- file.path(outputdir, paste0(abrev, "_eemsFiles"))
  dir.create(outputdir_i, recursive = T, showWarnings = F)

  subset <- data[data[,spscol] == i,]
  famFILE <- paste0(abrev, bedprefix, ".fam")
  bimFILE <- paste0(abrev, bedprefix, ".bim")
  bedFILE <- paste0(abrev, bedprefix, ".bed")

  if (file.exists(file.path(beddir, famFILE))){
    cat("\nCopying bed files to", outputdir_i, "folder...")
    system(paste("cp", file.path(beddir, bedFILE), file.path(beddir, bimFILE), file.path(beddir, famFILE), outputdir_i))
    cat("\nReading fam file...")
    fam <- read.table(file.path(beddir, famFILE))
    samples <- fam[,1]
    cat("\nSamples:", paste(samples, collapse = ", "))
    subset <- subset[subset[,indcol] %in% samples,]
    subset[,indcol]  <- factor(subset[,indcol], levels = samples)
    subset <- subset[order(subset[,indcol]),]

    nsamples <- length(samples)
    nsites <- as.integer(system(paste("wc -l", file.path(beddir, bimFILE), "| awk '{print $1}'"), intern = TRUE))
    cat("\nNumber of samples: ", nsamples, "\nNumber of sites: ", nsites, "\n")
    
    ## write .coords
    coords <- as.data.frame(subset[, c(latcol, loncol)])
    if (translon){
      coords[,2][coords[,2] < 0] = coords[,2][coords[,2] < 0] + 360
    }
    
    if (any(is.na(coords))){
      warning("NA coordintates found. removing individuals with no coordinates...")
      remove_samples <- samples[is.na(coords[,1])]
      coords <- coords[!is.na(coords[,1]),]
      out_excluded <- file.path(outputdir_i, paste0(abrev, bedprefix, "_excluded_indv.txt"))
      write.table(cbind(remove_samples, remove_samples), file = out_excluded, append = FALSE, quote = F, col.names = F, row.names = F)
      system(paste0(plinkpath, " --bfile ", file.path(outputdir_i, paste0(abrev, bedprefix)), " --remove ", out_excluded,  " --make-bed --out ", file.path(outputdir_i, paste0(abrev, bedprefix, "_wcoords"))))
      system(paste("mv", file.path(outputdir_i, paste0(abrev, bedprefix, "_wcoords.bed")), file.path(outputdir_i, paste0(abrev, bedprefix, ".bed"))))
      system(paste("mv", file.path(outputdir_i, paste0(abrev, bedprefix, "_wcoords.bim")), file.path(outputdir_i, paste0(abrev, bedprefix, ".bim"))))
      system(paste("mv", file.path(outputdir_i, paste0(abrev, bedprefix, "_wcoords.fam")), file.path(outputdir_i, paste0(abrev, bedprefix, ".fam"))))

      nsites <- as.integer(system(paste("wc -l", file.path(outputdir_i, paste0(abrev, bedprefix)), "| awk '{print $1}'"), intern = TRUE))
    }
    cat("\nWritting .coord file...")
    write.table(coords, quote = F, col.names = F, row.names = F , file = file.path(outputdir_i, paste0(abrev, bedprefix, ".coord")))
    
    ## write .diffs
    if (!file.exists(file.path(outputdir_i, paste0(abrev, bedprefix, ".diffs")))) {
      cat("\nGenerating .diffs file...")
      system(paste0(bed2diffspath, " --bfile ", file.path(outputdir_i, paste0(abrev, bedprefix))))
    }

    ## write .outer
    if (!file.exists(file.path(outputdir_i, paste0(abrev, bedprefix, ".outer")))) {
      cat("\nGenerating .outer file...")
      outerfile <- file.path(outputdir_i, paste0(abrev, bedprefix))
      system(paste("Rscript", makeouterpath, file.path(distshdir,paste0(abrev, ".shp")), regshfile, outerfile, TRUE), intern = TRUE)
    }
    
    ## write control file
    cat("\nCreating control file...")
    dataPATH <- file.path(outputdir_i, paste0(abrev, bedprefix))
    outputPATH <- file.path(outputdir_i, paste0(abrev, "-EEMS-ndemes", ndemes, "-chain", chain))
    if (dir.exists(outputPATH) & !overwrite){
      warning("Output directory already exists, creating a new mcmc chain...")
      chain = chain + 1
      outputPATH <- file.path(dataPATH, paste0(bedprefix, "-EEMS-ndemes", ndemes, "-chain", chain))
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
              file =  file.path(outputdir_i, paste0(abrev, "_params-chain", chain, ".ini")), row.names = FALSE)

  } else {
    cat("\nNo BED files were found for", i, "\n")
  }
}
