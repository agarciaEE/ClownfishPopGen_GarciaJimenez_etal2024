#!/usr/bin/bash  --login
# example run: sbatch -J AKA_vcfStats 5A_vcf_Statistics.sh AKA_sample.vcf.gz . MajorMinor 0.012 AKA 1 popList.txt

#SBATCH --time=24:00:00
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=60G
#SBATCH --output=%x_%j_%a.out
#SBATCH --error=%x_%j_%a.err

# ARGUMENTS:
########################################
vcfFile=${1-AKA_sample.vcf.gz} 	# Vcf file
workindir=${2-.} 		# Working directory
inputfileDir=${3-MajorMinor}	# Input directory containing the vcf file
SamplingProportion=${4-0.012} 	# proportion of variants we want to retain. Default = 0.012
outPrefix=${5-AKA} 		# prefix to add to output files
Nbootstraps=${6-1} 		# number of bootstraps
popListFile=${7-popList.txt}		# File containing sample ID in the first columns and pop ID in the 2nd to carry AF, het and HWE test per population
##########################################

echo "Using following arguments:

VCF file:        "$vcfFile"
Working directory:	"$workindir"
Input file directory:	"$inputfileDir"
Population info file:	"$popListFile"
Sampling proportion (r):        "$SamplingProportion"
Number of bootstraps:	"$Nbootstraps"
Output PREFIX:  "$outPrefix

echo "
If settings are wrong, please look at instructions and arguments to pass

"

####################
#### load module ###
####################

module load gcc/10.4.0 r/4.2.1 vcftools/0.1.14 bcftools/1.15.1 bzip2/1.0.8 samtools/1.15.1

cd $workindir
mkdir -p "$inputfileDir"/vcfStats
mkdir -p temp # write temporary vcf files in scratch


############################################################
## Relounch the script as array with the number of jobs corresponding to the number of bootstraps
if [[ "$SLURM_ARRAY_TASK_ID" == "" ]]; then

	if [ ! -f "$inputfileDir"/vcfStats/"$outPrefix"_nVariants.txt ]; then
		# Get number of variants in the vcf file
		echo "Calculating number of variants..."
		nSNP=$(bcftools view -H $inputfileDir/$vcfFile | wc -l)

		echo "Number of variants: "$nSNP
		echo $nSNP > "$inputfileDir"/vcfStats/"$outPrefix"_nVariants.txt
	fi
	bcftools query -l $inputfileDir/$vcfFile > "$inputfileDir"/vcfStats/"$outPrefix"_samplesID.txt
	# Relaunch this script as an array
     	exec sbatch -J "$SLURM_JOB_NAME" --array=1-$Nbootstraps $0 $vcfFile $workindir $inputfileDir $SamplingProportion $outPrefix $Nbootstraps $popListFile
fi
############################################################
BOOTSTRAP_NUM=$SLURM_ARRAY_TASK_ID

####################
#### load module ###
####################

module load gcc/10.4.0 r/4.2.1 vcftools/0.1.14 bcftools/1.15.1 bzip2/1.0.8 samtools/1.15.1
module load miniconda3/4.10.3
conda activate ATLAS

##########################
##### Vcf Statistics #####
##########################

if [[ "$SamplingProportion" -lt 1 ]]; then
# randomly sample variants of the vcf file
SUBSET_VCF=temp/$(echo $vcfFile | sed "s/.vcf.gz/_subset$BOOTSTRAP_NUM.vcf.gz/g")

if [ ! -f "$SUBSET_VCF" ]; then
	echo "Sampling VCF subset num. $BOOTSTRAP_NUM as $SUBSET_VCF..."
	bcftools view "$inputfileDir/$vcfFile" | vcfrandomsample -r "$SamplingProportion" | bcftools view -Oz -o "$SUBSET_VCF"  # requires vcflib (not installed in curnagl) and uses python 3.1 (downgrades newer version of phyton).
	# compress vcf
	#pbgzip "$SUBSET_VCF"
	# adding the compressed extension to the variable SUBSET_VCF
	#SUBSET_VCF="$SUBSET_VCF".gz
	# index vcf
	bcftools index "$SUBSET_VCF"
fi

else

SUBSET_VCF="$inputfileDir"/"$vcfFile"

fi

# create output file name according to the prefix and the boostrap number.
outFile="$inputfileDir"/vcfStats/"$outPrefix"_b"$BOOTSTRAP_NUM"

# Calculate mean depth per individual
############################
echo "Calculating mean depth per individual..."
vcftools --gzvcf "$SUBSET_VCF" --depth --out "$outFile"

# Calculate mean depth per site
############################
echo "Calculating mean depth per site..."
vcftools --gzvcf "$SUBSET_VCF" --site-mean-depth --out "$outFile"

# Calculate site quality
############################
echo "Calculating site quality..."
vcftools --gzvcf "$SUBSET_VCF" --site-quality --out "$outFile"

# Calculate proportion of missing data per individual
############################
echo "Calculating proportion of missing data per individual..."
vcftools --gzvcf "$SUBSET_VCF" --missing-indv --out "$outFile"

# Calculate proportion of missing data per site
############################
echo "Calculating proportion of missing data per site..."
vcftools --gzvcf "$SUBSET_VCF" --missing-site --out "$outFile"

if [ -f "$popListFile" ];then

	echo "Subsetting by population..."

	spsDir=$(echo "$outPrefix" | sed 's/_Chr[0-9]*//g')

	bash SNPcall_pipeline/0_split_samples_by_pop.sh "$popListFile" "$spsDir"

	popFilesDir="$spsDir"_popFiles

	for file in "${popFilesDir}"/*.txt; do

		# Extract the population name from the text file's name
  		population=$(basename "$file" .txt | sed 's/_Chr[0-9]//g' | cut -d '_' -f 2)

		# Create a subset VCF file for each population
		pop_subset_vcf="temp/"$outPrefix"_"$population"_subset"$BOOTSTRAP_NUM".vcf.gz"

  		# Extract sample IDs from the text file
  		#sample_ids=$(cat "$file")

		echo "Subsetting $vcfFile selecting samples from $population population..."

		# Use vcftools to subset the vcf file by sample ids
		vcf-subset -c "$file" "$SUBSET_VCF" | bcftools view -Oz -o "$pop_subset_vcf"
  		# Use 'bcftools view' to subset the VCF file by sample IDs
  		#bcftools view -Oz -o "$pop_subset_vcf" -S "$file" "$SUBSET_VCF" # gives ERROR

  		# Index the subset VCF file
  		bcftools index "$pop_subset_vcf"

  		echo "Subset VCF file '$pop_subset_vcf' has been created for population '$population'."

		popOutFile="$inputfileDir"/vcfStats/"$outPrefix"_"$population"_b"$BOOTSTRAP_NUM"

		# Calculate allele frequency
		############################
		echo "Calculating allele frequency for $population population..."
		vcftools --gzvcf "$pop_subset_vcf" --freq2 --out "$popOutFile" --max-alleles 2

		# Calculate heterozygosity and inbreeding coefficient per individual
		####################################################################
		echo "Calculating heterozygosity and inbreeding coefficient per individual in $population population..."
		vcftools --gzvcf "$pop_subset_vcf" --het --out "$popOutFile"

                # Perform HWE test
                ##################
                echo "Assessing Hardy-Weinberg Equilibrium per site in $population population..."
                vcftools --gzvcf "$pop_subset_vcf" --hardy --out "$popOutFile"

		# Remove pop vcf file
		#####################
#		echo "Removing population subset vcf file..."
#		rm "$pop_subset_vcf"
	done

	echo "Subset VCF files have been created for all populations."

else

# Calculate allele frequency
############################
echo "Calculating allele frequency..."
vcftools --gzvcf "$SUBSET_VCF" --freq2 --out "$outFile" --max-alleles 2

# Calculate heterozygosity and inbreeding coefficient per individual
############################
echo "Calculating heterozygosity and inbreeding coefficient per individual..."
vcftools --gzvcf "$SUBSET_VCF" --het --out "$outFile"

# Perform HWE test
##################
echo "Assessing Hardy-Weinberg Equilibrium per site in $population population..."
vcftools --gzvcf "$SUBSET_VCF" --hardy --out "$outFile"

fi

# Remove subset vcf to reduce unnecessary memory usage
######################################################
#echo "Removing subset vcf file..."
#rm "$SUBSET_VCF"

echo "All analyses have successfully finished."


