#!/usr/bin/bash
# example run: sbatch 2A_Run_BWA_RefIndexing.sh AperculaReference Amphiprion_percula.Nemo_v1.dna_sm.primary_assembly.AllChrom.fa AperGenome

#SBATCH --job-name=BWA_index
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=1:00:00

# Use this script only if the Reference indexing (for BWA) has not been performed yet

module load gcc/10.4.0 bwa/0.7.17

pathtoReference=${1:-AperculaReference}
RefGenom=${2:-Amphiprion_percula.Nemo_v1.dna_sm.primary_assembly.AllChrom.fa}
RefOutputName=${3:-AperGenome}
 
echo "Using following input:
Path to reference genome:		"$pathtoReference"
Reference Genome Name:			"$RefGenom"
Reference output name for indexing: 	"$RefOutputName
echo "If settings are wrong, please look at instructions"

bwa index $pathtoReference"/"$RefGenom -p $pathtoReference"/"$RefOutputName


