#!/usr/bin/bash
# example run: sbatch 1C_Concat_PhasedVCF.sh samplesList.txt . MajorMinor/Phased 8

#SBATCH --job-name=Concat_Phased
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=%x_%j_%a.out
#SBATCH --error=%x_%j_%a.err
#SBATCH --time=08:00:00

##############################
# Argument to pass:
samplesList=${1-samplesList.txt}		# File with samples names
workindir=${2-.}                                # Working directory. Default: Current Directory
inputfileDir=${3-Phased}                        # Name of the directory with the input mapping files (should be in working directory)
threads=${4-8}					# Number of threads
##############################

echo "Using following arguments:

List of samples:	"$samplesList"
Working directory:	"$workindir"
VCF file directory (if relative path, from "$workindir" :  "$inputfileDir"
threads:	"$threads
echo "
If settings are wrong, please look at instructions and arguments to pass

"

############################################################
## Relounch the script as array with the number of jobs corresponding to the number of files
nline=$(wc -l $samplesList | awk '{print $1}')

if [[ "$SLURM_ARRAY_TASK_ID" == "" ]]; then
     # Relaunch this script as an array
     exec sbatch --array=1-$nline $0 $samplesList $workindir $inputfileDir $threads
fi

############################################################

####################
# LOAD MODULES
####################

module load gcc/10.4.0 shapeit4/4.1.3 bcftools/1.15.1

####################
# BEGINNNIN OF SCRIPT
####################

# Get name of samples for array
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p $samplesList)

# Print the name of the sample
echo "Working on:        "$sample

# Change to working directory
cd $workindir

mkdir -p tmp

ls -1v "$inputfileDir"/"$sample"_[0-9]*_phased.vcf > tmp/"$sample"_files.txt

bcftools concat -f tmp/"$sample"_files.txt -o "$inputfileDir"/"$sample"_phased.vcf --threads "$threads"

bcftools view -Oz -o "$inputfileDir"/"$sample"_phased.vcf.gz "$inputfileDir"/"$sample"_phased.vcf

bcftools index "$inputfileDir"/"$sample"_phased.vcf.gz

echo "End of script."
