#!/usr/bin/bash --login
# example run: sbatch -J AKA_plinkPCA 6A_PlinkPCA.sh AKA.vcf.gz . MajorMinor Plink AKA AKA.prune.in AKA_popList.txt

#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=30G
#SBATCH --time=05:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# ARGUMENTS:
###############
vcfFile=${1-AKA_samples.vcf.gz}                 # source vcf file
workindir=${2-.}                                # Working directory
inputfileDir=${3-MajorMinor}                    # Input vcf file directory
outputfileDir=${4-Plink}	                # Output file directory
outPrefix=${5-AKA}				# Prefix of output file
PruneInFile=${6-AKA.prune.in} 			# Prunned file
popFile=${7-AKA_popList.txt}			# File with population info (no headers) (1st column: samplenames; 2nd column: population)
###############

####################
#### load module ###
####################

module load gcc/10.4.0 r/4.2.1 vcftools/0.1.14 bcftools/1.15.1 bzip2/1.0.8 samtools/1.15.1 plink-ng
#module load gcc/10.4.0 miniconda3/4.10.3 r/4.2.1
#conda activate PopGen
PLINK=/work/FAC/FBM/DBC/nsalamin/software/nsalamin/clownfishes/agarciaj/plink

cd "$workindir"

echo "Using following parameters :

VCF file :             			"$vcfFile"
Working directory :               	"$workindir"
Input file directory (relative to $workindir) :               "$inputfileDir"
Output file directory :                    "$outputfileDir"
Output file prefix :             	"$outPrefix"
Prunned file :      "$PruneInFile"
Population information file :		"$popFile

echo "
If settings are wrong, please look at instructions and arguments to pass

"

mkdir -p "$outputfileDir"

#plink2 --vcf "$inputfileDir"/"$vcfFile" --double-id --allow-extra-chr \
#	--set-missing-var-ids @:# \
#	--extract "$outputfileDir"/"$PruneInFile" \
#	--make-bed --pca --out "$outputfileDir"/"$outPrefix"
# plink2 gives and error when less than 50 samples

echo "Carrying out Plink PCA..."

# --allow-extra-chr \
$PLINK --vcf "$inputfileDir"/"$vcfFile" --double-id --chr-set 24 no-xy no-mt \
        --set-missing-var-ids @:# \
        --extract "$outputfileDir"/"$PruneInFile" \
        --make-bed --pca --out "$outputfileDir"/"$outPrefix"

eigenvecFile="$outputfileDir"/"$outPrefix".eigenvec
eigenvalFile="$outputfileDir"/"$outPrefix".eigenval

echo "Ploting PCA.."

Rscript SNPcall_pipeline/6B_plotPlinkPCA.R "$eigenvecFile" "$eigenvalFile" "$popFile" "$outputfileDir" "$outPrefix"

echo "Script has finished."
