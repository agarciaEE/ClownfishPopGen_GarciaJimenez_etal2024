#!/usr/bin/bash
# example run: sbatch -J AKA_vcfStats_report 5B_produce_vcfStatistics_report.sh . SNPcall_pipeline vcfStats AKA_report AKA all all AKA_popList.txt 0.25

#SBATCH --time=06:00:00
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# ARGUMENTS:
########################################
workindir=${1-.} 			# Working directory
RscriptFolder=${2-SNPcall_pipeline}     # Folder containing the Rscript 5B_output_vcfStats.R
inputfileDir=${3-vcfStats} 		# intput and output file directory
outputDir=${4-AKA_report}		# Output directory to create
PREFIX=${5-AKA} 			# prefix to search target files and to add to output files
CHR=${6-all}				# Chromosome number to analyse. Default "all" which will analyse all byChr vcf files or in case they are all togethe rin one vcf, it will analyse that one.
BOOTSTRAPS=${7-all} 			# number of bootstraps to use in the report
POPS=${8-}				# Populations to be analyzed. Comma separated or text file with two columns (1srt sample name and 2nd pop name)
PROP_CHR_TH=${9-0.25}			# proportion of chromomes that need to be outlier to be considered outlier
##########################################

##############################

echo "Using following arguments:

Working directory:	"$workindir"
Input directory (if relative path, from "$workindir" :	"$inputfileDir"
Files prefix:	"$PREFIX"
Chromosome list:        "$CHR"
Number of subsets:     "$BOOTSTRAPS"
Populations:	"$POPS"
Outlier detection threshold:  "$PROP_CHR_TH

echo "
If settings are wrong, please look at instructions and arguments to pass

"
####################
#### load module ###
####################

module load gcc/10.4.0 r/4.2.1
#module load miniconda3/4.10.3
#conda activate ATLAS

cd "$workindir"

outputDir="$inputfileDir"/"$outputDir"

mkdir -p "$outputDir"

echo "Producing vcf statistics report and saving them in $outputDir folder..."

Rscript "$RscriptFolder"/5B_output_vcfStats.R "$inputfileDir" "$outputDir" "$PREFIX" "$BOOTSTRAPS" "$CHR" "$POPS" "$PROP_CHR_TH"

echo "Reports have been produced successfully."


