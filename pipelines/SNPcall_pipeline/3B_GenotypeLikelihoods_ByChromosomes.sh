#!/usr/bin/bash --login
#example run: sbatch 3B_GenotypeLikelihoods_ByChromosomes.sh AllSamples.txt . FilteredMapping 100000 atlas Chromosomes.txt Apercula

#SBATCH --job-name=GLF
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH --output=GLFchrom_%a.out
#SBATCH --error=GLFchrom_%a.err
#SBATCH --time=32:00:00

##############################
# Argument to pass:
SampleList=${1-AllSamples.txt}		# Name of the file with sample names. Default: AllSamples.txt
workindir=${2-.}				# Working directory. Default: Current Directory
inputfileDir=${3-FilteredMapping}			# Name of the directory with the input Filtered Mapping files (should be in workindir). Default: FilteredMapping
windsize=${4-100000}				# Window size for GLF computation. Default:100 Kb, as it reduces computation time and memory. ATLAS default would be 1Mb
AtlasExe=${5-atlas}	#Atlas executable
FileChromosomes=${6-Chromosomes.txt}   # File with chromosomes names
suffix=${7-Apercula}			# reference genome suffix
##############################


############################################################
## Relounch the script as array with the number of jobs corresponding to the number of files
nline=$(wc -l $SampleList | awk '{print $1}')


if [[ "$SLURM_ARRAY_TASK_ID" == "" ]]; then
     # Relaunch this script as an array
     exec sbatch --array=1-$nline $0 $SampleList $workindir $inputfileDir $windsize $AtlasExe $FileChromosomes $suffix
fi

############################################################

echo "Using following arguments:

List of Samples:	"$SampleList"
Working directory:	"$workindir"
Mapping Files directory (if relative path, from "$workindir" :	"$inputfileDir"
Window size for ATLAS:	"$windsize"
Atlas executable: 	"$AtlasExe"
Chromosome File:	"$FileChromosomes"
Reference suffix: "$suffix


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
# BEGINNNING OF SCRIPT
####################

# Get name of samples for array
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SampleList)

# Print the name of the sample
echo "Working on sample:	"$sample

# Change to working directory
cd $workindir

####################
# SPLIT BAM by Chromosome
####################

mkdir -p "./$inputfileDir/"$sample

# not necessary as all bam finally processed together and do not differ in chromosome coding
#samtools view -h ./FilteredMapping/"$sample".BWA.Apercula.Sort.Filt_mergedReads.bam | sed -e '/^@SQ/s/SN\:/SN\:chr/' -e '/^[^@]/s/\t/\tchr/2' | samtools view -bS - > ./FilteredMapping/"$sample"/"$sample".BWA.Apercula.Sort.Filt_mergedReads.bam
#cp ./FilteredMapping/"$sample"/"$sample".BWA.Apercula.Sort.Filt_mergedReads.bam ./FilteredMapping/"$sample"/"$sample".BWA.Apercula.Sort.Filt_mergedReads.NewHead.bam
#bamtools split -in ./FilteredMapping/"$sample"/"$sample".BWA.Apercula.Sort.Filt_mergedReads.NewHead.bam -reference

cp "$inputfileDir"/"$sample".BWA."$suffix".Sort.Filt_mergedReads.bam "$inputfileDir"/"$sample"/"$sample".BWA."$suffix".Sort.Filt_mergedReads.bam
bamtools split -in "$inputfileDir"/"$sample"/"$sample".BWA."$suffix".Sort.Filt_mergedReads.bam -reference -refPrefix Chr

# Depending on the reference the extracted chromosome reference will have already 'chr'. In order to run the GLF, we need to remove it
for i in $(ls "$inputfileDir"/"$sample"/*Chr*)
do
mv "$i" "$(echo $i | sed 's/Chrchr/Chr/g')"
done

####################
# Create index files and Generate Genotype likelihood files
####################

cat $FileChromosomes | while read line
do
  samtools index ./FilteredMapping/"$sample"/"$sample".BWA."$suffix".Sort.Filt_mergedReads."$line".bam
  echo "GLF for chom "$line
  mkdir -p "./GLF/"$line
  $AtlasExe task=GLF bam=./FilteredMapping/"$sample"/"$sample".BWA."$suffix".Sort.Filt_mergedReads."$line".bam window="$windsize" out=./GLF/"$line"/"$sample"."$line".BWA."$suffix".Sort.Filt_mergedReads
done

####################
# END
####################

echo "script for sample "$sample" has finished"
