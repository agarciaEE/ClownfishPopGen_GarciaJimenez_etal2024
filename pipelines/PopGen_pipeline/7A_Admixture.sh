#!/usr/bin/bash --login
# example run: sbatch -J AKA_Admixture 7A_Admixture.sh AKA.bed . Plink/Raw Admixture AKA "1-6" 10 200

#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --time=3-00:00:00
#SBATCH --output=%x_%j_%a.out
#SBATCH --error=%x_%j_%a.err

# ARGUMENTS:
###############
BEDFile=${1-AKA.bed}   		              	# source BED file
workindir=${2-.}                                # Working directory
inputfileDir=${3-Plink}                    # Input vcf file directory
outputfileDir=${4-Admixture}                        # Output file directory
outPrefix=${5-AKA}                              # Prefix of output file
nK=${6-"1-6"}					# Number of K to run Admixture
CV=${7-10}					# X-fold Cross-Validation
nB=${8-200}					# Number of bootstraps
###############

####################
#### load module ###
####################

cd "$workindir"

echo "Using following parameters :

VCF file :                              "$vcfFile"
Working directory :                     "$workindir"
Input file directory (relative to $workindir) :               "$inputfileDir"
Output file directory :                    "$outputfileDir"
Output file prefix :                    "$outPrefix"
Number of K :				"$nK"
CV fold :				"$CV"
Number of bootstraps :			"$nB

echo "
If settings are wrong, please look at instructions and arguments to pass

"

############################################################
## Relounch the script as array with the number of jobs corresponding to the number of bootstraps
if [[ "$SLURM_ARRAY_TASK_ID" == "" ]]; then

        # Relaunch this script as an array
        exec sbatch -J "$SLURM_JOB_NAME" --array=$nK $0 $BEDFile $workindir $inputfileDir $outputfileDir $outPrefix $nK $CV
fi
############################################################
K=$SLURM_ARRAY_TASK_ID

mkdir -p "$outputfileDir"

module load gcc/10.4.0 miniconda3/4.10.3 r/4.2.1
conda activate PopGen

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
awk '{$1="0";print $0}' "$inputfileDir"/"${BEDFile%.bed}".bim > "$inputfileDir"/"${BEDFile%.bed}".bim.tmp
mv "$inputfileDir"/"${BEDFile%.bed}".bim.tmp "$inputfileDir"/"${BEDFile%.bed}".bim

# generate .map file with plink
plink --bfile "$inputfileDir"/"${BEDFile%.bed}" --recode --out "$inputfileDir"/"${BEDFile%.bed}"

admixture --cv="$CV" -B"$nB" -j"$SLURM_CPUS_PER_TASK" "$inputfileDir"/"$BEDFile" "$K" > "$outputfileDir"/"$outPrefix"_K"$K"_admixture_log.out
