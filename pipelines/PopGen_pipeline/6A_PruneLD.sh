#!/usr/bin/bash --login
# example run: sbatch -J AKA_pruneLD 6A_PruneLD.sh AKA.vcf.gz . MajorMinor MajorMinor AKA 50 10 0.1

#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=30G
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# ARGUMENTS:
###############
vcfFile=${1-AKA_samples.vcf.gz}                 # source vcf file
workindir=${2-.}                                # Working directory
inputfileDir=${3-MajorMinor}                    # Input vcf file directory
outputfileDir=${4-MajorMinor}                   # Output file directory
outPrefix=${5-AKA}				# Prefix of output file
WindowSize=${6-50}   				# Window Size (Kb) of the analysis
SlidingWindow=${7-10}        			# Size of sliding window to compute LD (bp)
LDthreshold=${8-0.1}				# r2 threshold. Threshold of linkage to tolerate. variants with r2 greater than threshold will be prunned out.
###############

####################
#### load module ###
####################

#module load gcc/10.4.0 r/4.2.1 vcftools/0.1.14 bcftools/1.15.1 bzip2/1.0.8 samtools/1.15.1 plink-ng
module load gcc/10.4.0 miniconda3/4.10.3 r/4.2.1
conda activate PopGen

cd "$workindir"

echo "Using following parameters :

VCF file :             			"$vcfFile"
Working directory :               	"$workindir"
Input file directory (relative to $workindir) :               "$inputfileDir"
Output file directory :                    "$outputfileDir"
Output file prefix :             	"$outPrefix"
Window Size :           		"$WindowSize"
Sliding Window Size :           	"$SlidingWindow"
Linkage Disequilibrium threshold :      "$LDthreshold

echo "
If settings are wrong, please look at instructions and arguments to pass

"

mkdir -p "$outputfileDir"

#plink2 --vcf "$inputfileDir"/"$vcfFile" --double-id --allow-extra-chr \
#	--set-missing-var-ids @:# \
#	--indep-pairwise "$WindowSize" "$SlidingWindow" "$LDthreshold" --out "$outputfileDir"/"$outPrefix"
# plink2 gives an error when less than 50 samples

plink --vcf "$inputfileDir"/"$vcfFile" --double-id --chr-set 24 no-xy no-mt \
	 --set-missing-var-ids @:# \
	--indep-pairwise "$WindowSize" "$SlidingWindow" "$LDthreshold" --out "$outputfileDir"/"$outPrefix"



