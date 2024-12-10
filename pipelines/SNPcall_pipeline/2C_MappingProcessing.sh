#!/usr/bin/bash  --login
# example run: sbatch 2C_MappingProcessing.sh samplesnames.txt . MappingOutput 30 atlas Apercula

#SBATCH --job-name=Processing
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=%x_%a.out
#SBATCH --error=%x_%a.err
#SBATCH --time=24:00:00

##############################
# Argument to pass:
SampleList=${1-samplesnames.txt}			# Name of the file with sample names. Default: AllSamples.txt
workindir=${2-.}					# Working directory. Default: Current Directory
inputfileDir=${3-MappingOutput}			# Name of the directory with the input mapping files (should be in workindir). Default: MappingOutput
QualityThreshold=${4-30}				# Mapping quality threshold. Default 30
AtlasExe=${5-atlas}	#Atlas executable. Default
suffix=${6-Apercula}                            # Suffix to add with infor about the reference genome used.
##############################

############################################################
## Relounch the script as array with the number of jobs corresponding to the number of files
nline=$(wc -l $SampleList | awk '{print $1}')


if [[ "$SLURM_ARRAY_TASK_ID" == "" ]]; then
     # Relaunch this script as an array
     exec sbatch --array=1-$nline $0 $SampleList $workindir $inputfileDir $QualityThreshold $AtlasExe $suffix
fi

############################################################

echo "Using following arguments:

List of Samples:	"$SampleList"
Suffix: 		"$suffix"
Working directory:	"$workindir"
Mapping Files directory (if relative path, from "$workindir" :	"$inputfileDir"
Mapping quality threshold:	"$QualityThreshold"
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
# Generate file with read name to later extract read group
####################

mkdir -p temp

samtools view "./"$inputfileDir"/"$sample".BWA."$suffix".Sorted.bam" | head -n 1 > "./temp/"$sample".header.txt"

####################
# Get information for Readgroups
####################

#string="$(cut -d '	' -f1 "./temp/"$sample".header.txt")"  ## doesn't work
string="$(cat "./temp/"$sample".header.txt" | sed -r 's/\s*//g')"
flowcell="$(cut -d':' -f1 <<<"$string")"
lane="$(cut -d':' -f2 <<<"$string")"
barcode="$(cut -d':' -f3 <<<"$string")"
RGID=$flowcell"."$lane"."$barcode
RGLB=$sample".Lib1"
RGPL="Illumina"
RGPU=$RGID"."$RGLB
RGSM=$sample

####################
# Add Read Groups
####################
mkdir -p FilteredMapping

echo $RGSM
echo $RGPU
echo $RGPL
echo $RGID
echo $RGLB

#java -Xmx8g  -jar /dcsrsoft/spack/arolle/v1.0/spack/opt/spack/linux-rhel8-zen2/gcc-10.4.0/picard-2.26.2-uymkctyol3ii7b3uzmphfm5gqu3djacm/bin/picard.jar AddOrReplaceReadGroups I="./"$inputfileDir"/"$sample".BWA.Apercula.Sorted.bam" O="./FilteredMapping/"$sample".BWA.Apercula.Sorted.RG.bam" RGID="$RGID" RGLB="$RGLB" RGPL="$RGPL" RGPU="$RGPU" RGSM="$RGSM" 
# does not work (killed reason:out-of-memory error)
#picard AddOrReplaceReadGroups I="./"$inputfileDir"/"$sample".BWA.Apercula.Sorted.bam" O="./FilteredMapping/"$sample".BWA.Apercula.Sorted.RG.bam" RGID="$RGID" RGLB="$RGLB" RGPL="$RGPL" RGPU="$RGPU" RGSM="$RGSM"
samtools addreplacerg -r ID:$RGID -r LB:$RGLB -r PL:$RGPL -r PU:$RGPU -r SM:$RGSM -o "./FilteredMapping/"$sample".BWA."$suffix".Sorted.RG.bam" "./"$inputfileDir"/"$sample".BWA."$suffix".Sorted.bam"

####################
# Index the bam file with read group
####################

samtools index "./FilteredMapping/"$sample".BWA."$suffix".Sorted.RG.bam"

####################
# Filter BAM
# Keep only files with FLAGS 2 (read mapped in proper pair) and remove files with  FLAGS 256 (not a primary alignment)
####################

samtools view -b -f 2 -F 256 -q $QualityThreshold -@ 2 "./FilteredMapping/"$sample".BWA."$suffix".Sorted.RG.bam" > "./FilteredMapping/"$sample".BWA."$suffix".Sort.Filt.bam"

####################
# Index the filtered bam
####################

samtools index "./FilteredMapping/"$sample".BWA."$suffix".Sort.Filt.bam"

####################
# Run Atlas readOverlap
# Generate statistics on overlapping read pairs with ATLAS.
####################

mkdir -p ProcessingStatistics

$AtlasExe task=readOverlap bam="./FilteredMapping/"$sample".BWA."$suffix".Sort.Filt.bam" out="./ProcessingStatistics/"$sample".BWA."$suffix".Sort.Filt"

####################
# Run atlas_readMerge
#  Generate a new mapping file by setting the quality of one of the overlapping segment paired reads to zero with atlas.
####################

$AtlasExe task=mergeReads bam="./FilteredMapping/"$sample".BWA."$suffix".Sort.Filt.bam" out="./FilteredMapping/"$sample".BWA."$suffix".Sort.Filt"

mkdir -p ProcessingStatistics

####################
# Generate Statistics with bamtools
####################

bamtools stats -in "./FilteredMapping/"$sample".BWA."$suffix".Sort.Filt_mergedReads.bam" > "./ProcessingStatistics/"$sample".BWA."$suffix".Sort.Filt_mergedReads.Stats"

####################
# Generate additional mapping statistics (insert size)
####################

picard CollectInsertSizeMetrics I="./FilteredMapping/"$sample".BWA."$suffix".Sort.Filt_mergedReads.bam" O="./ProcessingStatistics/"$sample".Filt.MergedReads.InsertMetrics.txt" H="./ProcessingStatistics/"$sample".BWA."$suffix".Filt.MergedReads.InsertPlot.pdf"

####################
# Generate additional mapping statistics wth ATLAS (insert size)
####################

$AtlasExe task=BAMDiagnostics bam="./FilteredMapping/"$sample".BWA."$suffix".Sort.Filt_mergedReads.bam" out="./ProcessingStatistics/"$sample".BWA."$suffix".Sort.Filt_mergedReads"

####################
# END
####################

echo "script for sample "$sample" has finished"
echo "If you wanna run ATLAS createDepthMask, assessSoftClipping or pileup, please look at the end of the 2C_MappingProcessing.sh script"

#######################################################################################
#######################################################################################

# following task are not performed by default, as they are computationally intensive and the output is not essential, as they are mostly for statistics purposes.
#A filter on depth will be done on the VCF in later steps. Softclipping is not considered in callers. If you still want to run the following task, please run it separately

#######################################################################################
#######################################################################################

####################
# Run atlas_depth_mask
# Create a depth mask in BED format
# If still want to run, please check and set the parameters minDepthForMasking and maxDepthForMask

#minDepthForMasking=0
#maxDepthForMask=40
#$AtlasExe task=createDepthMask minDepthForMask=$minDepthForMasking maxDepthForMask=$maxDepthForMasking bam="./FilteredMapping/"$sample".BWA.Apercula.Sort.Filt_mergedReads.bam" out="./FilteredMapping/"$sample".BWA.Apercula.Sort.Filt_mergedReads"

####################
# Run Assess soft clipping
# Assess soft clipping bases
# Soft-clipped bases are not considered by callers

#$AtlasExe task=assessSoftClipping bam="./FilteredMapping/"$sample".BWA.Apercula.Sort.Filt_mergedReads.bam" out="./ProcessingStatistics/"$sample".BWA.Apercula.Sort.Filt_mergedReads"

####################
# Generate additional mapping statistics wth ATLAS (pileup)
# !!!!! By default, we do not run it as much memory and computation !!!!!

#$AtlasExe task=pileup bam="./FilteredMapping/"$sample"BWA.Apercula.Sort.Filt_mergedReads.bam" out="./ProcessingStatistics/"$sample"BWA.Apercula.Sort.Filt_mergedReads"







