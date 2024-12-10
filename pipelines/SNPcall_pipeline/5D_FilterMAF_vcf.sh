#!/usr/bin/bash --login
# example run: sbatch -J AKA_2ndFilterVCF 5D_Filter2_vcf.sh AKA_samples.vcf.gz . MajorMinor/Filter1 MajorMinor/Filter2 vcfStats/AKA_report/AKA_filtering_parameters.txt vcfStats/AKA_report/AKA_remove_ind.txt

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
outputfileDir=${4-MajorMinor/Filtered}                   # Output filtered vcf file directory
filterFilename=${5-}   # Filename containing the parameters to filter
removeIndFilename=${6-}        # Filename containing the name of individuals to remove from the vcf file
###############

####################
#### load module ###
####################

module load gcc/10.4.0 r/4.2.1 vcftools/0.1.14 bcftools/1.15.1 bzip2/1.0.8 samtools/1.15.1
module load miniconda3
conda activate ATLAS

cd $workindir

mkdir -p $outputfileDir

if [ -n "$filterFilename" ]
then
  echo "Using filtering parameters provided in file $filterFilename"
  MAF=$(cat "$filterFilename" | grep "MAF" | awk '{print $2}')
fi

if [ ! -n "$MAF" ];then MAF=0.02; fi

echo "Using following filtering parameters :

Minimum Allele Frequency threshold :       	"$MAF"
Removing individuals specified file :		"$removeIndFilename

echo "
If settings are wrong, please look at instructions and arguments to pass

"

# If filtering parameters file has multiple MAF by population
#MAF=$(cat $filterFilename | grep "MAF" | awk '{print $2}' | awk '{ sum += $1 } END { printf(sum/NR) }' RS=" ")

# perform the filtering with vcftools

echo "Filtering vcf file following the given parameters..."

suffix=_maf"$MAF"

outPrefix=$(echo $VCF_IN | sed "s/.vcf.gz/$suffix/g")

if [ -n "$removeIndFilename" ]; then
        vcftools --gzvcf "$inputfileDir/$VCF_IN" --maf "$MAF" --remove "$removeIndFilename" --recode --recode-INFO-all --out "$outputfileDir"/"$outPrefix"
else
	vcftools --gzvcf "$inputfileDir/$VCF_IN" --maf "$MAF" --recode --recode-INFO-all --out "$outputfileDir"/"$outPrefix"
fi

bcftools view -Oz -o "$outputfileDir"/"$outPrefix".vcf.gz "$outputfileDir"/"$outPrefix".recode.vcf

bcftools index "$outputfileDir"/"$outPrefix".vcf.gz

rm "$outputfileDir"/"$outPrefix".recode.vcf

echo "Calculating the remaining number of variants..."

echo "Filter2" "$(bcftools view -H "$outputfileDir/$outPrefix.vcf.gz" | wc -l)" >> "$outputfileDir"/"$outPrefix"_nvariants.txt

echo "Filtering has finished."
