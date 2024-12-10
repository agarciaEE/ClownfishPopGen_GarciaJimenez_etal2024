#!/usr/bin/bash  --login

#SBATCH --job-name=GLF
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH --output=outscripts/GLF_%a.out
#SBATCH --error=outscripts//GLF_%a.err
#SBATCH --time=24:00:00

##############################
# Argument to pass:
SampleList=${1-AllSamples.txt}		# Name of the file with sample names. Default: AllSamples.txt
workindir=${2-.}				# Working directory. Default: Current Directory
inputfileDir=${3-FilteredMapping}			# Name of the directory with the input Filtered Mapping files (should be in workindir). Default: FilteredMapping
windsize=${4-100000}				# Window size for GLF computation. Default:100 Kb, as it reduces computation time and memory. ATLAS default would be 1Mb
AtlasExe=${5-atlas}	#Atlas executable via conda
##############################

############################################################
## Relounch the script as array with the number of jobs corresponding to the number of files
nline=$(wc -l $SampleList | awk '{print $1}')


if [[ "$SLURM_ARRAY_TASK_ID" == "" ]]; then
     # Relaunch this script as an array
     exec sbatch --array=1-$nline $0 $SampleList $workindir $inputfileDir $windsize $AtlasExe
fi

############################################################

echo "Using following arguments: 

List of Samples:	"$SampleList"
Working directory:	"$workindir"
Mapping Files directory (if relative path, from "$workindir" :	"$inputfileDir"
Window size for ATLAS:	"$windsize"
Atlas executable: 	"$AtlasExe

echo "
If settings are wrong, please look at instructions and arguments to pass

"
####################
# LOAD MODULES
####################

module load gcc/10.4.0 bwa/0.7.17 python/3.9.13 samtools/1.15.1 picard/2.26.2 bamtools/2.5.2 r/4.2.1
module load miniconda3/4.10.3
conda activate ATLAS


####################
# BEGINNNIN OF SCRIPT
####################

# Get name of samples for array
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SampleList)

# Print the name of the sample
echo "Working on sample:	"$sample

# Change to working directory
cd $workindir

####################
# GENERATE Genotype likelihood files
####################

mkdir -p GLF

$AtlasExe task=GLF bam="./FilteredMapping/"$sample".BWA.Apercula.Sort.Filt_mergedReads.bam" window="$windsize" out="./GLF/"$sample".BWA.Apercula.Sort.Filt_mergedReads"

####################
# END
####################

echo "script for sample "$sample" has finished"










