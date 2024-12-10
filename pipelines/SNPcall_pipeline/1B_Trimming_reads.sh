#!/usr/bin/bash
#example run: sbatch SNPcall_pipeline/1B_Trimming_reads.sh samplenames.txt . rawData adapters/adapters.fa 20 50 0

#SBATCH --job-name=QCTrim
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --output=%x_%a.out
#SBATCH --error=%x_%a.err
#SBATCH --time=24:00:00

##############################
# Argument to pass:
SampleList=${1-samplenames.txt}          # Name of the file with sample names. Default: AllSamples.txt
workindir=${2-.}                                # Working directory. Default: Current Directory
inputfileDir=${3-rawData}                       # Name of the directory with the input fastq.gz files (should be in workindir). Default: RawData
                                                # The fastq files should be in the format $sample_R1.fastq.gz and $sample_R2.fastq.gz
adapterfiles=${4-adapters/adapters.fa}  # Path and name of the files where adapter sequences are. Default: adapters/adapters.fa
minQuality=${5-20}                              # Minimum quality score for trimming. Default: 20
minLength=${6-50}                               # Minimum length of reads: Default: 50
headcrop=${7-0}
##############################

##############################
## Relounch the script as array with the number of jobs corresponding to the number of files
#nline=$(wc -w $1 | awk '{print $1}') # if samples in list of elements
nline=$(wc -l $SampleList | awk '{print $1}') # if samples in lines


if [[ "$SLURM_ARRAY_TASK_ID" == "" ]]; then
     # Relaunch this script as an array
     exec sbatch --array=1-$nline $0 $SampleList $workindir $inputfileDir $adapterfiles $minQuality $minLength $headcrop
fi
##############################

echo "Using following arguments: 

List of Samples:        "$SampleList"
Working directory:	"$workindir"
Raw data directory (if relative path, from "$workindir" :	"$inputfileDir"
Adapter files:  "$adapterfiles"
Minimum quality of bases for trimming:  "$minQuality"
Minimum length of reads:        "$minLength"
Number of basis cut beginning of reads  "$headcrop

echo "
If settings are wrong, please look at instructions and arguments to pass
"

module load gcc/10.4.0 fastqc/0.11.9 trimmomatic/0.39

# Get name of samples for array
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SampleList) #if SamplesList is store by line
#sample=$(cat $SampleList | awk -v i=$SLURM_ARRAY_TASK_ID '{print $i}')

# Print the name of the sample
echo "Working on sample:        "$sample

# Change to working directory
cd $workindir

####################
# Run quality report with fastqc for each read file 
####################

# make directory fore fastqc raw data if does not exist
mkdir -p fastQC_RawData

# generate fastqc report
fastqc -o fastQC_RawData $inputfileDir"/"$sample"_R1.fastq.gz"
fastqc -o fastQC_RawData $inputfileDir"/"$sample"_R2.fastq.gz"

####################
# Run the trimming and adapter check with trimmomatic
####################

# make directory fore Trimmed reads and Summary of trimming if does not exist
mkdir -p TrimmedReads
mkdir -p Trimming_Summary

trimmomatic PE -phred33 -summary "Trimming_Summary/Summary_"$sample".txt" $inputfileDir"/"$sample"_R1.fastq.gz" $inputfileDir"/"$sample"_R2.fastq.gz" "TrimmedReads/"$sample"_R1.Trimmed.Paired.fastq.gz" "TrimmedReads/"$sample"_R1.Trimmed.Unpaired.fastq.gz" "TrimmedReads/"$sample"_R2.Trimmed.Paired.fastq.gz" "TrimmedReads/"$sample"_R2.Trimmed.Unpaired.fastq.gz" ILLUMINACLIP:$adapterfiles:2:30:10 SLIDINGWINDOW:4:15 LEADING:$minQuality TRAILING:$minQuality MINLEN:$minLength HEADCROP:$headcrop


####################
# Run quality report with fastqc for each Trimmed read file 
####################
mkdir -p fastQC_TrimmedData
fastqc -o fastQC_TrimmedData "TrimmedReads/"$sample"_R1.Trimmed.Paired.fastq.gz"
fastqc -o fastQC_TrimmedData "TrimmedReads/"$sample"_R2.Trimmed.Paired.fastq.gz"
fastqc -o fastQC_TrimmedData "TrimmedReads/"$sample"_R1.Trimmed.Unpaired.fastq.gz"
fastqc -o fastQC_TrimmedData "TrimmedReads/"$sample"_R2.Trimmed.Unpaired.fastq.gz"


# END
echo "script for sample "$sample" has finished"




