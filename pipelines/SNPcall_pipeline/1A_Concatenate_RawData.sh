#!/usr/bin/bash

#SBATCH --array=1-5
#SBATCH --job-name=CatLA
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --output=CatRawData_%a.out
#SBATCH --error=CatRawData_%a.err
#SBATCH --time=1:00:00

sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p RerunSamples.txt)

cd rawData

cat $sample"_L1_R1_"* $sample"_L2_R1_"* $sample"_L3_R1_"* $sample"_L4_R1_"* > $sample"_R1.fastq.gz"
cat $sample"_L1_R2_"* $sample"_L2_R2_"* $sample"_L3_R2_"* $sample"_L4_R2_"* > $sample"_R2.fastq.gz"

rm $sample"_L1_R1_"* $sample"_L2_R1_"* $sample"_L3_R1_"* $sample"_L4_R1_"* $sample"_L1_R2_"* $sample"_L2_R2_"* $sample"_L3_R2_"* $sample"_L4_R2_"*
