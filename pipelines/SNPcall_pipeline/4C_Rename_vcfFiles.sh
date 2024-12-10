#!/usr/bin/bash
# example run: sbatch 4C_Rename_vcfFiles.sh . MajorMinor vcfFiles.txt ATLAS

#SBATCH --job-name=rename_vcf
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=32:00:00

##############################
# Argument to pass:
workindir=${1-.}				# Working directory. Default: Current Directory
inputDir=${2-MajorMinor}			# Name of the directory with the input MajorMinor files (should be in workindir). Default: MajorMinor
vcfFiles=${3-POL_vcfFiles.txt}   # File with vcf filenames to concatenate
match=${4-ATLAS} # matching pattern to rename the vcfFiles after renaming samples ID. File will have inserted the string "renamed_" to differentiate from the original one.
##############################


############################################################
## Relounch the script as array with the number of jobs corresponding to the number of files
nline=$(cat $vcfFiles | wc -l)

if [[ "$SLURM_ARRAY_TASK_ID" == "" ]]; then
     # Relaunch this script as an array
     exec sbatch --array=1-$nline $0 $workindir $inputDir $vcfFiles $match
fi

############################################################

echo "Using following arguments:

List of vcf files:	"$vcfFiles"
Working directory:	"$workindir"
Vcf Files directory (if relative path, from "$workindir" :	"$inputDir"
match: "$match

echo "
If settings are wrong, please look at instructions and arguments to pass

"
####################
# LOAD MODULES
####################

module load gcc/10.4.0 bwa/0.7.17 python/3.9.13 vcftools/0.1.14 bcftools/1.15.1 samtools/1.15.1 picard/2.26.2 bamtools/2.5.2 r/4.2.1

####################
# BEGINNNING OF SCRIPT
####################

# Get name of samples for array
line=$(sed -n ${SLURM_ARRAY_TASK_ID}p $vcfFiles)

# Print the name of the sample
echo "Working on file:        "$vcfFile

# Change to working directory
cd $workindir

# delete "renamed_$vcfFiles" if exist as we are gonna append the renamed vcfFiles to it
rm renamed_$vcfFiles

###
# Decompress vcf file
###
vcfline=$(echo $line | sed 's/.gz//')
echo "Uncompressing $line to $vcfline..."
bcftools view $line -O v -o $vcfline

###
# Rename sample id
###
filename="${vcfline%%${match}*}renamed_${match}${vcfline##*${match}}"
bcftools query -l $vcfline | cut -d '/' -f4 | cut -d '.' -f 1 > $inputDir/samplesID.txt
echo "Renaming samples ID as:"
cat $inputDir/samplesID.txt
echo "New filename: $filename"
bcftools reheader -s $inputDir/samplesID.txt -o $filename $vcfline
echo "Removing $vcfline..."
rm $vcfline

###
# Compress vcf file back
###
echo "Compressing "$filename" to "$filename".gz..."
bcftools view $filename -O z -o $filename.gz

###
# Add new filename to vcfFiles list
###
echo "Appending "$filename" to renamed_"$vcfFiles

echo $filename.gz >> $inputDir/renamed_$vcfFiles
echo "Removing $filename..."
rm $filename

####################
# END
####################

echo "Files have been renamed with success."
