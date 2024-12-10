#!/usr/bin/bash  --login
# example run: sbatch 7.combine_CROSS_POP.sh Crossfiles.txt . MSMC/CROSS_results MSMC/POP_results MSMC/Combined_results

#SBATCH --job-name=combine_MSMC
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=02:00:00

##############################
# Argument to pass:
inputFile=${1-Crossfiles.txt}   # File containing a list of CROSS-population msmc2 results
workinDir=${2-.}                 # Working directory. Detault: current directory (.)
inputCROSSDir=${3-MSMC/CROSS_results} # Directory where the cross-population msmc2 files are stored. Default: MSMC/inputCROSS
inputPOPDir=${4-MSMC/POP_results}   # Directory where the populations msmc2 files are stored. Default: MSMC/inputPOP
outputDir=${5-MSMC/Combined_results}  # Output directory
############################################################

echo "Using following arguments:

Input file:         							"$inputFile"
Working directory:							"$workinDir"
Cross-population msmc2 results directory (relative to "$workinDir"):    "$inputCROSSDir"
Population msmc2 results directory (relative to "$workinDir"):        	"$inputPOPDir"
Output directory:							"$outputDir

echo "
If settings are wrong, please look at instructions and arguments to pass

"

####################
# LOAD MODULES
####################

module load gcc python/3.9.13
MSMCTOOLS=/work/FAC/FBM/DBC/nsalamin/software/nsalamin/clownfishes/agarciaj/msmc-tools/

cd $workinDir

mkdir -p $outputDir

#####################
# BEGINNNIN OF SCRIPT
#####################

while read -r file; do

	echo "Working on file $file..."

	# get species prefix
        spPrefix=$(echo "$file" | cut -d '_' -f1)

	# get run name
	runName=$(echo "$file" | cut -d '.' -f3 | cut -d '.' -f1)

        # get chromosome info
        chr=$(echo "$file" | cut -d '.' -f2 | cut -d '.' -f1)

	# get populations
	pop1=$(echo "$file" | cut -d '_' -f2 | cut -d '.' -f1 | cut -d '-' -f1)
	pop2=$(echo "$file" | cut -d '_' -f2 | cut -d '.' -f1 | cut -d '-' -f2)

	# get each population msmc2 file
	msmc_pop1=$(ls "$inputPOPDir"/"$spPrefix"_"$pop1"."$chr"."$runName".*.final.txt)
	msmc_pop2=$(ls "$inputPOPDir"/"$spPrefix"_"$pop2"."$chr"."$runName".*.final.txt)
        msmc_cross="$inputCROSSDir"/"$file"

	all_files=true
	if [ ! -f "$msmc_pop1" ]; then
		echo "$pop1 population msmc2 file not found."
		all_files=false
	fi
        if [ ! -f "$msmc_pop2" ]; then
                echo "$pop2 population msmc2 file not found."
                all_files=false
        fi
        if [ ! -f "$msmc_cross" ]; then
                echo "$pop1-$pop2 cross-population msmc2 file not found."
                all_files=false
        fi

	# combine
	if "$all_files"; then
		out_filename="$outputDir"/"$spPrefix"_"$pop1"-"$pop2"."$chr"."$runName".msmc2.final.combined.txt
		$MSMCTOOLS/combineCrossCoal.py "$msmc_cross" "$msmc_pop1" "$msmc_pop2" > "$out_filename"
		echo "MSMC2 outputs combined in $out_filename."
	else
		echo "Cross-population msmc2 $spPrefix between $pop1 and $pop2 failed."
	fi

done < "$inputFile"
