#!/usr/bin/bash
# example run: sbatch 2B_Mapping_generating_sorted_BAM.sh samplesnames.txt . TrimmedReads AperculaReference/AperGenome remove Apercula

#SBATCH --job-name=BWA
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --output=%x_%a.out
#SBATCH --error=%x_%a.err
#SBATCH --time=32:00:00

##############################
# Argument to pass:
SampleList=${1-samplesnames.txt}                  # Name of the file with sample names. Default: AllSamples.txt
workindir=${2-.}                                        # Working directory. Default: Current Directory
inputfileDir=${3-TrimmedReads}                  # Name of the directory with the input fastq.gz files (should be in workindir). Default: TrimmedReads
ReferenceGenome=${4-AperculaReference/AperGenome} # Name (and location) of the index of the reference genome
TEMP=${5-remove}                               # Do remove SAM and unsorted bam files, by default. If something else, do not remove them!
suffix=${6-Apercula}				# Suffix to add with infor about the reference genome used.
##############################


############################################################
## Relounch the script as array with the number of jobs corresponding to the number of files
nline=$(wc -l $SampleList | awk '{print $1}')


if [[ "$SLURM_ARRAY_TASK_ID" == "" ]]; then
     # Relaunch this script as an array
     exec sbatch --array=1-$nline $0 $SampleList $workindir $inputfileDir $ReferenceGenome $removeTEMP $suffix
fi

############################################################

echo "Using following arguments: 

List of Samples:        "$SampleList"
Working directory:	"$workindir"
Trimmed data directory (if relative path, from "$workindir" :   "$inputfileDir"
Reference genome (index):	"$ReferenceGenome"
Remove or not 'temporary' sam and unsorted bam files:   "$removeTEMP

echo "
If settings are wrong, please look at instructions and arguments to pass

"
####################
# LOAD MODULES
####################

module load gcc/10.4.0 bwa/0.7.17 python/3.9.13 samtools/1.15.1 picard/2.26.2 bamtools/2.5.2 r/4.2.1

####################
# BEGINNNIN OF SCRIPT
####################


# Get name of samples for array
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SampleList)

# Print the name of the sample
echo "Working on sample:        "$sample

# Change to working directory
cd $workindir

####################
# Run BWA
####################

mkdir -p MappingOutput
mkdir -p summary_mapping

if [[ "$TEMP" = "pipe" ]]
then
	bwa mem -M -t 2 $ReferenceGenome $inputfileDir"/"$sample"_R1.Trimmed.Paired.fastq.gz" $inputfileDir"/"$sample"_R2.Trimmed.Paired.fastq.gz" \
        | samtools view -Sb -@ 2 | samtools sort -@ 2 > "./MappingOutput/"$sample".BWA."$suffix".Sorted.bam"
else
	bwa mem -M -t 2 $ReferenceGenome $inputfileDir"/"$sample"_R1.Trimmed.Paired.fastq.gz" $inputfileDir"/"$sample"_R2.Trimmed.Paired.fastq.gz" > "./MappingOutput/"$sample".BWA."$suffix".sam"

	####################
	# Convert to BAM
	####################
	echo "Converting to BAM..."
	samtools view -Sb -@ 2 "./MappingOutput/"$sample".BWA.Apercula.sam" > "./MappingOutput/"$sample".BWA."$suffix".bam"

	####################
	# Sort BAM (by coordinate)
	####################
	echo "Sorting..."
	samtools sort -@ 2 "./MappingOutput/"$sample".BWA.Apercula.bam" > "./MappingOutput/"$sample".BWA."$suffix".Sorted.bam"
fi

####################
# Index BAM
####################
samtools index "./MappingOutput/"$sample".BWA."$suffix".Sorted.bam"

####################
# Generate additional mapping statistics
####################

bamtools  stats -in "./MappingOutput/"$sample".BWA."$suffix".Sorted.bam"  > "summary_mapping/"$sample"MappingStats.txt"
picard CollectInsertSizeMetrics I="./MappingOutput/"$sample".BWA."$suffix".Sorted.bam" O="summary_mapping/"$sample".InsertMetrics.txt" H="summary_mapping/"$sample".InsertPlot.pdf"

####################
# Remove sam and unsorted bam
####################

if [[ "$TEMP" == "remove" ]]
then
  echo "Removing sam and unsorted bam files."
  rm "./MappingOutput/"$sample".BWA."$suffix".sam"
  rm "./MappingOutput/"$sample".BWA."$suffix".bam"

fi

####################
# END
####################

echo "script for sample "$sample" has finished"




