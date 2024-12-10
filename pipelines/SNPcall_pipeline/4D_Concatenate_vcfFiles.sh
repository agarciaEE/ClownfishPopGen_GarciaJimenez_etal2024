#!/usr/bin/bash --login
# example run: sbatch -J AKA_concat_vcf 4D_Concatenate_vcfFiles.sh . MajorMinor vcfFiles.txt . AKA_ATLAS_majorMinor.vcf ATLAS false

#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=60G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=3-00:00:00

##############################
# Argument to pass:
workindir=${1-.}				# Working directory. Default: Current Directory
inputDir=${2-POL_MajorMinor}			# Name of the directory with the input Filtered Mapping files (should be in workindir). Default: FilteredMapping
vcfFiles=${3-renamed_POL_vcfFiles.txt}   # File with vcf filenames to concatenate
outputDir=${4-MajorMinor}
outputfile=${5-POL_ATLAS_MajorMinor.vcf}
match=${6-ATLAS} # matching pattern to rename the vcfFiles after renaming samples ID. File will have inserted the string "renamed_" to differentiate from the original one.
rename_vcf=${7-false}
##############################

echo "Using following arguments:

List of vcf files:	"$vcfFiles"
Working directory:	"$workindir"
Vcf Files directory (if relative path, from "$workindir" :	"$inputDir"
Output vcf file directory:	"$outputDir"
Vcf filename: "$outputfile

echo "
If settings are wrong, please look at instructions and arguments to pass

"
####################
# LOAD MODULES
####################

module load gcc/10.4.0 bwa/0.7.17 python/3.9.13 vcftools/0.1.14 bcftools/1.15.1 samtools/1.15.1 picard/2.26.2 bamtools/2.5.2 r/4.2.1
module load miniconda3

conda activate ATLAS

####################
# BEGINNNING OF SCRIPT
####################

# Change to working directory
cd $workindir

mkdir -p $outputDir

PREFIX=$(echo $outputfile | cut -d '_' -f 1)

if [ "$rename_vcf" = "true" ]; then
	# delete "renamed_$vcfFiles" if exist as we are gonna append the renamed vcfFiles to it
	rm renamed_$vcfFiles

	# for each vcf file: decompress, rename and compress again to be able to concatenate them as they have the chromosome number as part of the sample id
	cat $vcfFiles | sort -V | while read line
	do
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
		bcftools query -l $vcfline | cut -d '/' -f4 | cut -d '.' -f 1 > $inputDir/"$PREFIX"_samplesID.txt
		echo "Renaming samples ID as:"
		cat $inputDir/"$PREFIX"_samplesID.txt
		echo "New filename: $filename"
		bcftools reheader -s $inputDir/"$PREFIX"_samplesID.txt -o $filename $vcfline
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
	done
	vcfFiles=$inputDir/renamed_$vcfFiles
fi

###
# Concatenate vcf Files
###

#bcftools concat -f $vcfFiles -O z -o $outputDir/$outputfile
#vcf-concat -f $vcfFiles | pbgzip -c > $outputDir/$outputfile

# Compression with pbgzip seems to give issues to be read by bcftools downstream
# I added this alternative in case files are not compressed properly

vcf-concat -f $vcfFiles > $outputDir/$outputfile

bcftools view -Oz -o $outputDir/$outputfile.gz $outputDir/$outputfile

bcftools index $outputDir/$outputfile.gz

####################
# END
####################

echo "Files have been concatenated with success."
