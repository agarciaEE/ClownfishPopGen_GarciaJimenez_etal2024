#!/usr/bin/bash --login
# example run: sbatch -J AKA_Majorminor 4B_Run_MajorMinor_ByChromosomes.sh  ALL_POL_Samples.txt . GLF AKA ATLAS BWA.Apercula.Sort.Filt_mergedReads.glf.gz atlas Chromosomes.txt

#SBATCH --array=1-24
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --output=%x_ByChrom_%a.out
#SBATCH --error=%x_ByChrom_%a.err
#SBATCH --time=72:00:00

##############################
# Argument to pass:
SampleList=${1-ALL_POL_Samples.txt}			# Name of the file with sample names. Default: AllSamples.txt
workindir=${2-.}				# Working directory. Default: Current Directory
GLFDirectory=${3-GLF}					# Directory with GLF files
PrefixOut=${4-AKA}				# Prefix to add in output files
OutfileSuffix=${5-ATLAS}		# Output Suffix
GLFSuffix=${6-BWA.Apercula.Sort.Filt_mergedReads.glf.gz}   #Suffix for mapping files
AtlasExe=${7-atlas}	#Atlas executable
FileChromosomes=${8-Chromosomes.txt}  #Files with Chromosomes names

##############################

echo "Using following arguments:

List of Samples:	"$SampleList"
Working directory:	"$workindir"
GLF data directory (if relative path, from "$workindir" :	"$GLFDirectory"
GLF files suffix:	"$GLFSuffix"
GLF file prefix: 	"GLFFilePrefixOut"
Output file Suffix:	"OutfileSuffix"
Atlas executable: 	"$AtlasExe"
Chromosome File:  "$FileChromosomes

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

# Get name of Chromosomes for array
chr=$(sed -n ${SLURM_ARRAY_TASK_ID}p $FileChromosomes)

# Print the name of the chromosome
echo "Working on chromosome:	"$chr

####################
# Create a file with the name of the GLF files, depending on the samples and for each chromosome
####################

cd $workindir

# Create GLF Files name for each chromosome
GLFFileName="$PrefixOut"_GLFFiles."$chr".txt
echo $GLFFileName
rm -f $GLFFileName

while IFS= read -r line
do
	linename="$workindir"/"$GLFDirectory"/"$chr"/"$line"."$chr"."$GLFSuffix"
	echo "$linename"
	echo "$linename" >> $GLFFileName
done < $SampleList

####################
# Run Minor/Major to get vcf
####################

mkdir -p MajorMinor

OutfileName="$PrefixOut"_"$chr"_"$OutfileSuffix"
$AtlasExe task=majorMinor glf="$GLFFileName" out=MajorMinor/"$OutfileName"








