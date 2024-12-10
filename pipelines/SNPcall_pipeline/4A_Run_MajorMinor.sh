#!/usr/bin/bash  --login

#SBATCH --job-name=MajorMinor
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --output=outscripts/MajorMinor_%a.out
#SBATCH --error=outscripts/MajorMinor_%a.err
#SBATCH --time=72:00:00

##############################
# Argument to pass: 
SampleList=${1-AllSamples.txt}			# Name of the file with sample names. Default: AllSamples.txt
workindir=${2-.}				# Working directory. Default: Current Directory
GLFDirectory=${3-GLF}					# Directory with GLF files
GLFFileName=${4-GLFFiles.txt}				# GLF file name
OutfileName=${5-ATLAS_majorMinor.vcf.gz}		# OutputFileName
GLFSuffix=${6-.BWA.Apercula.Sort.Filt_mergedReads.glf.gz}
AtlasExe=${7-atlas}	#Atlas executable

##############################



echo "Using following arguments: 

List of Samples:	"$SampleList"
Working directory:	"$workindir"
GLF data directory (if relative path, from "$workindir" :	"$GLFDirectory"
GLF files suffix:	"$GLFSuffix"
GLF file name: 	"$GLFFileName"
Output file name:	"$OutfileName"
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

####################
# Create a file with the name of the GLF files, depending on the samples
####################

cd $workindir

while IFS= read -r line
do
	linename=$workindir"/"$GLFDirectory"/"$line$GLFSuffix
	echo "$linename" >> $GLFFileName
done < $SampleList

####################
# Run Minor/Major to get vcf
####################

$AtlasExe task=majorMinor glf=$GLFFileName out=$OutfileName













