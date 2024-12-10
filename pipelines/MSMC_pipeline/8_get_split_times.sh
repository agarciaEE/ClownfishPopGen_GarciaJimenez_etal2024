#!/usr/bin/bash --login
# example run: sbatch 8_get_split_times.sh PREFIX . MSMC/Combined_results MSMC/Combined_results PREFIX_msmc 5 2 4e-8 5

#SBATCH --job-name=getSplitTimes
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=00:20:00

##############################
# Argument to pass:
prefix=${1-AKA}   		# Prefix
workinDir=${2-.}                 # Working directory. Detault: current directory (.)
inputDir=${3-MSMC/Combine_results} # Directory where the combined cross and pop msmc2 files are stored. Default: MSMC/Combined_results
outputDir=${4-MSMC/Combined_results}  # Output directory
outputfile=${5-msmc_result}	# Output filename
dropRecent=${6-5}		# Number of entries to drop from recent estimates
dropPast=${7-2}			# Number of entries to drop from past estimates
mu=${8-4e-8}			# mutation ratr
gen=${9-5}			# generation time
############################################################

echo "Using following arguments:

Input file:         							"$inputFile"
Working directory:							"$workinDir"
Combined msmc2 results directory (relative to "$workinDir"):    	"$inputDir"
Output directory:							"$outputDir"
Output filename:							"$outputfile"
Dropped recent entries:							"$dropRecent"
Dropped past entries:							"$dropPast"
Mutation rate:								"$mu"
Generation Time:							"$gen

echo "
If settings are wrong, please look at instructions and arguments to pass

"

####################
# LOAD MODULES
####################

module load gcc miniconda3

conda activate MSMC

cd $workinDir

mkdir -p $outputDir

#####################
# BEGINNNIN OF SCRIPT
#####################

msmcfiles=$(ls "$inputDir"/"$prefix"*)

python MSMC_pipeline/get_split_times.py "$msmcfiles" --out "$outputDir"/"$outputfile" \
			--mu "$mu" --gen "$gen" --drop_recent "$dropRecent" --drop_past "$dropPast"

echo "Script done."
