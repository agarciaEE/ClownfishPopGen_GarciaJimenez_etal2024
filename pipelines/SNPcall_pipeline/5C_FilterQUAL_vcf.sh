#!/usr/bin/bash --login
# example run: sbatch -J AKA_1stFilterVCF 5B_Filter_vcf.sh AKA_samples.vcf.gz . . MajorMinor/1_filter vcfStats/AKA_report/AKA_filtering_parameters.txt vcfStats/AKA_report/AKA_remove_ind.txt

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
  MIN_QUALITY=$(cat "$filterFilename" | grep "MIN_QUAL" | awk '{print $2}')
  MIN_DEPTH=$(cat "$filterFilename" | grep "MIN_DEPTH" | awk '{print $2}')
  MAX_DEPTH=$(cat "$filterFilename" | grep "MAX_DEPTH" | awk '{print $2}')
  MISSING_THRESHOLD=$(cat "$filterFilename" | grep "MISS_THRESHOLD" | awk '{print $2}')
fi

if [ ! -n "$MIN_QUALITY" ];then MIN_QUALITY=30; fi
if [ ! -n "$MIN_DEPTH" ];then MIN_DEPTH=5; fi
if [ ! -n "$MAX_DEPTH" ];then MAX_DEPTH=50; fi
if [ ! -n "$MISSING_THRESHOLD" ];then MISSING_THRESHOLD=0.9; fi

echo "Using following filtering parameters :

Minimum quality threshold :       	"$MIN_QUALITY"
Minimum depth threshold	:		"$MIN_DEPTH"
Maximum depth threshold	:		"$MAX_DEPTH"
Missing data proportion threshold :		"$MISSING_THRESHOLD"
Removing individuals specified file :		"$removeIndFilename

echo "
If settings are wrong, please look at instructions and arguments to pass

"

# perform the filtering with vcftools

echo "Filtering vcf file following the given parameters..."

# create file prefix where is gonna be adding all filtering info as suffixes along the way
missVal=$(echo "100 - $MISSING_THRESHOLD * 100" | bc)
missVal=$(echo ${missVal%.*})

suffix=_q"$MIN_QUALITY"dp"$MIN_DEPTH"_"$MAX_DEPTH"miss"$missVal"
outPrefix=$(echo $VCF_IN | sed "s/.vcf.gz/$suffix/g")

echo "Performing first filter to $VCF_IN and creating filtered file: $outPrefix.recode.vcf.gz..."

if [ -n "$removeIndFilename" ];
then
	vcftools --gzvcf "$inputfileDir/$VCF_IN" --remove-indels \
		--remove "$removeIndFilename" \
        	--max-missing "$MISSING_THRESHOLD" \
        	--minQ "$MIN_QUALITY" \
        	--min-meanDP "$MIN_DEPTH" --max-meanDP "$MAX_DEPTH" \
        	--minDP "$MIN_DEPTH" --maxDP "$MAX_DEPTH" --recode --recode-INFO-all --out "$outputfileDir"/"$outPrefix"

else
	vcftools --gzvcf "$inputfileDir/$VCF_IN" \
		--max-missing "$MISSING_THRESHOLD" --minQ "$MIN_QUALITY" \
		--min-meanDP "$MIN_DEPTH" --max-meanDP "$MAX_DEPTH" \
		--minDP "$MIN_DEPTH" --maxDP "$MAX_DEPTH" --remove-indels \
		--recode --recode-INFO-all --out "$outputfileDir"/"$outPrefix"
fi

bcftools view -Oz -o "$outputfileDir"/"$outPrefix".vcf.gz "$outputfileDir"/"$outPrefix".recode.vcf

bcftools index "$outputfileDir"/"$outPrefix".vcf.gz

#pbgzip -c "$outputfileDir"/"$outPrefix".recode.vcf > "$outputfileDir"/"$outPrefix".vcf.gz

rm "$outputfileDir"/"$outPrefix".recode.vcf

echo "Calculating the remaining number of variants..."

echo "Filter1" "$(bcftools view -H "$outputfileDir/$outPrefix.vcf.gz" | wc -l)" >> "$outputfileDir"/"$outPrefix"_nvariants.txt

echo "Analysis has finished."

