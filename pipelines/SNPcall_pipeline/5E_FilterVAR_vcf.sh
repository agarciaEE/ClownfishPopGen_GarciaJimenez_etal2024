#!/usr/bin/bash --login
# example run: sbatch -J AKA_3rdFilterVCF 5E_Filter3_vcf.sh AKA_samples.vcf.gz . MajorMinor/Filter1 MajorMinor/Filter2 vcfStats/AKA_report/AKA.prunne.in.txt keep vcfStats/AKA_report/AKA_remove_ind.txt

#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=30G
#SBATCH --time=32:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# ARGUMENTS:
###############
VCF_IN=${1-AKA_samples.vcf.gz}                  # source vcf file
workindir=${2-.}                                # Working directory
inputfileDir=${3-MajorMinor}                    # Input vcf file directory
outputfileDir=${4-MajorMinor/Filtered}          # Output filtered vcf file directory
outPrefix=${5-AKA_samples_LDprunned}		# output file prefix
VariantsFile=${6-}   				# File containing the variants to keep or exclude (CHR\t$POS)
type=${7-keep}		        		# Specify wheter to 'keep' or 'exclude' the variants of the variants file (Default: "keep")
removeIndFilename=${8-}				# File containing the samples to remove
###############

echo "Using following filtering parameters :

VCF file :            "$VCF_IN"
Working directory:	"$workindir"
Inputfile directory:	"$inputfileDir"
Output directory:	"$outputfileDir"
Output file prefix:	"$outPrefix"
Variants file:	"$VariantsFile"
Filter type:	"$type"
Remove individuals file:	"$removeIndFilename

echo "
If settings are wrong, please look at instructions and arguments to pass

"

####################
#### load module ###
####################

module load gcc/10.4.0 r/4.2.1 vcftools/0.1.14 bcftools/1.15.1 bzip2/1.0.8 samtools/1.15.1
module load miniconda3
conda activate ATLAS

cd $workindir

mkdir -p $outputfileDir

if [ -f "$VariantsFile" ];then

   if [ "$type" == "exclude" ];then
	echo "Using $VariantsFile file to filter variants out."

  	vcftools --gzvcf "$inputfileDir/$VCF_IN" --exclude-positions "$VariantsFile" --recode --recode-INFO-all --out "$outputfileDir"/"$outPrefix"

   elif [ "$type" == "keep" ];then
        echo "Using $VariantsFile file to subset VCF file."

        vcftools --gzvcf "$inputfileDir/$VCF_IN" --positions "$VariantsFile" --recode --recode-INFO-all --out "$outputfileDir"/"$outPrefix"
   fi

else

   echo "File $VariantsFile not found."
   exit 1

fi

if [ -n "$removeIndFilename" ]; then

   echo "Removing samples specified in $removeIndFilename..."
   vcftools --vcf "$outputfileDir"/"$outPrefix".recode.vcf --remove "$removeIndFilename" --recode --recode-INFO-all --out "$outputfileDir"/"$outPrefix"

fi

bcftools view -Oz -o "$outputfileDir"/"$outPrefix".vcf.gz "$outputfileDir"/"$outPrefix".recode.vcf

bcftools index "$outputfileDir"/"$outPrefix".vcf.gz

rm "$outputfileDir"/"$outPrefix".recode.vcf

echo "Calculating the remaining number of variants..."

echo "Filter3" "$(bcftools view -H "$outputfileDir/$outPrefix.vcf.gz" | wc -l)" >> "$outputfileDir"/"$outPrefix"_nvariants.txt

echo "Filtering has finished."
