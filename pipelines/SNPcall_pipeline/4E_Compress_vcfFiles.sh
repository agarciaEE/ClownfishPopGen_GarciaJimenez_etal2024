#!/usr/bin/bash
# example run: sbatch 4E_Compress_vcfFiles.sh vcfFiles.txt . . .

#SBATCH --job-name=vcf_compress
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=30G
#SBATCH --output=%x_%a.out
#SBATCH --error=%x_%a.err
#SBATCH --time=32:00:00

##############################
# Argument to pass:
vcfFilesList=${1-vcfFiles.txt}		# filename containing the list of files
workindir=${2-.}				# Working directory. Default: Current Directory
inputfileDir=${3-.}			# Name of the directory with the input Filtered Mapping files (should be in workindir). Default: FilteredMapping
outputfileDir=${4-.}
##############################

echo "Using following arguments:

Filename with vcf files:	"$vcfFilesList"
Working directory:	"$workindir"
Vcf Files directory (if relative path, from "$workindir" :	"$inputfileDir"
Output vcf file directory:	"$outputfileDir

echo "
If settings are wrong, please look at instructions and arguments to pass

"

############################################################
## Relounch the script as array with the number of jobs corresponding to the number of files
nline=$(wc -l $vcfFilesList | awk '{print $1}')

if [[ "$SLURM_ARRAY_TASK_ID" == "" ]]; then
     # Relaunch this script as an array
     exec sbatch --array=1-$nline $0 $vcfFilesList $workindir $inputfileDir $outputfileDir
fi

############################################################

####################
# LOAD MODULES
####################

module load gcc/10.4.0 bcftools/1.15.1

####################
# BEGINNNING OF SCRIPT
####################

# Get name of samples for array
vcfFile=$(sed -n ${SLURM_ARRAY_TASK_ID}p $vcfFilesList)

# Print the name of the sample
echo "Working on vcf file:        "$vcfFile

# Change to working directory
cd $workindir

if ! file $vcfFile | grep -q compressed ; then

	bcftools view -Oz -o $inputfileDir/$vcfFile.gz $outputfileDir/$vcfFile

	bcftools index $outputfileDir/$vcfFile.gz

else

	echo "File is already compressed."

        bcftools index $inputfileDir/$vcfFile

fi

####################
# END
####################

echo "File have been compressed and indexed successfully."
