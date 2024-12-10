#!/usr/bin/env bash
# example run: sbatch 2_generate_globalmask.sh AperculaReference/Amphiprion_percula.Nemo_v1.dna_sm.primary_assembly.AllChrom.fa Aperc_reference

#SBATCH --job-name=make_reference_masks
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --output=%x_%a.out
#SBATCH --error=%x_%a.err
#SBATCH --time=24:00:00

##############################
# Argument to pass:
GENOME=${1-AclarkiiReference/A.clarkii_FinalAssembly.fasta} 	# reference genome to make masks
outfile_PREFIX=${2-AclarkiiGenome} 				# Prefix of the output file masks
outputDir=${3-MSMC/refmasks}					# output directory to create if doesn't exist
k=${4-35}							# k-mer
##############################

####################
# LOAD MODULES
####################
module load gcc python/3.9.13

msmctoolsPATH=/work/FAC/FBM/DBC/nsalamin/clownfish/agarciaj/msmc/msmc-tools

####################
# BEGINNNIN OF SCRIPT
####################

mkdir -p $outputDir

echo "Starting extraction of overlapping ${k}-mer subsequences"

$msmctoolsPATH/splitfa $GENOME $k | split -l 20000000
cat x* >> $outputDir/${outfile_PREFIX}_split.$k

echo "Aligning ${k}-mer reads to the genome with BWA, then converting to sam file"

# the genome needs to be indexed prior to this step-- if it has not already been indexed, run:
if [ -f "${GENOME}.bwt" ]; then
        echo "$GENOME already indexed"
else
        echo "indexing $GENOME"
        bwa index $GENOME
fi

echo "aligning reads to genome with BWA and converting to sam"
bwa aln -t 8 -R 1000000 -O 3 -E 3 ${GENOME} $outputDir/${outfile_PREFIX}_split.${k} > $outputDir/${outfile_PREFIX}_split.${k}.sai
bwa samse -f $outputDir/${outfile_PREFIX}_split.${k}.sam $GENOME $outputDir/${outfile_PREFIX}_split.${k}.sai $outputDir/${outfile_PREFIX}_split.${k}

echo "reads aligned, starting to generate rawMask"
gen_raw_mask.pl $outputDir/${outfile_PREFIX}_split.${k}.sam > $outputDir/${outfile_PREFIX}_rawMask.${k}.fa

echo "raw mask created as ${prefix}_rawMask.35.fa, now generating final mask with stringency r=50%"
gen_mask -l ${k} -r 0.5 $outputDir/${outfile_PREFIX}_rawMask.${k}.fa > $outputDir/${outfile_PREFIX}_mask.${k}.50.fa

echo "all done! final mask saved as ${outfile_PREFIX}_mask.${k}.50.fa"

$msmctoolsPATH/makeMappabilityMask.py  $outputDir/${outfile_PREFIX}_mask.${k}.50.fa $outputDir/$outfile_PREFIX

####################
# END
####################

echo "All masks have been created."



