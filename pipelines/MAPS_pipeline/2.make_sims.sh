#!/bin/bash

#SBATCH --job-name=make_sims_array
#SBATCH --output=make_sims_%A_%a.log
#SBATCH --error=make_sims_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-7 # Index for 8 species

# Load necessary modules
########################
module load gcc R # Ensure R is available on your system

# Define species array
######################
species=("AKA" "CLK" "AKY" "PRD" "CRP" "MEL" "SAN" "POL")

# Retrieve the species for this task
####################################
prefix=${species[$SLURM_ARRAY_TASK_ID]}
ibd_pipeline_dir="${prefix}_MAPSFiles/ibd_pipeline/"
lowerBnd=$1
upperBnd=$2

# Log the species being processed
#################################
echo "Processing species: $prefix"

# Run Rscript command
#####################
Rscript beagle2sims.R \
  -i "${ibd_pipeline_dir}finalqc_chr*.ibd" \
  -p "${prefix}" \
  -d "${ibd_pipeline_dir}" \
  -l "${lowerBnd}" \
  -u "${upperBnd}"

# Completion message
####################
echo "Sims file creation complete for ${prefix}."
