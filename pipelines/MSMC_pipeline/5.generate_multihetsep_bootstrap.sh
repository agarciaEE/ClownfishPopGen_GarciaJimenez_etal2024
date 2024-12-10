#!/usr/bin/env bash
# example run: sbatch -J AKA_multihet_bootstraps 5.generate_multihet_bootstrap.sh AKA . /work/FAC/FBM/DBC/nsalamin/clownfish/agarciaj/MSMC/inputPOP2/ /work/FAC/FBM/DBC/nsalamin/clownfish/agarciaj/MSMC/inputPOP2/bootstraps/ AKA_popList.txt 50 2000000 10

#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --output=%x_%a.out
#SBATCH --error=%x_%a.err
#SBATCH --time=06:00:00

##############################
# Argument to pass:
PREFIX=${1-AKA}         # Prefix. Default: AKA
workindir=${2-.}	# Working directory
inputDir=${3-/work/FAC/FBM/DBC/nsalamin/clownfish/agarciaj/MSMC/inputPOP2}
outputDir=${4-/work/FAC/FBM/DBC/nsalamin/clownfish/agarciaj/MSMC/inputPOP2/het_bootstraps}          # output directory.
POPFILE=${5-AKA_popList.txt}                      # PopList
nbootstraps=${6-50}	# number of bootstraps
chunk_size=${7-20000000} # size of bootstrap chunks
chunk_per_chr=${8-10}	# number of chunks to put inone chromosome in the bootstrap
##############################

POPS=(cat $POPFILE | awk '{print $2}' | sort | uniq)
npops=${#POPS[@]}

############################################################
## Relounch the script as array with the number of jobs corresponding to the number of populations
if [[ "$SLURM_ARRAY_TASK_ID" == "" ]] && [ -n "$POPFILE" ]; then
	# Relaunch this script as an array
     	exec sbatch -J "$SLURM_JOB_NAME" --array=1-$npops $0 $PREFIX $workinDir $inputDir $outputDir $POPFILE $nbootstraps $chunk_size $chunk_per_chr
fi
############################################################
POP=${POPS[$(echo "$SLURM_ARRAY_TASK_ID - 1" | bc)]}

####################
# LOAD MODULES
####################

module load gcc python/3.9.13
MSMCTOOLS=/work/FAC/FBM/DBC/nsalamin/software/nsalamin/clownfishes/agarciaj/msmc-tools/

echo "Using following arguments:

Prefix:                	"$PREFIX"
Population:		"$POP"
Working directory:	"$workindir"
Input directory:	"$inputDir"
Output directory:       "$outputDir"
Population info file:   "$POPFILE"
Number of bootstraps: 	"$nbootstraps"
Bootstrap chunks size:	"$chunk_size"
Number of chunks per chromosome:	"$chunk_per_chr
	
echo "
If settings are wrong, please look at instructions and arguments to pass

"
mkdir -p $outputDir

#input for the bootstrapping
BS_INPUT=("$inputDir"/"$PREFIX"_"$POP".chr*.multihet.txt)
nchr=${#BS_INPUT[@]}

#output from the bootstrapping 
BS_OUTPUT=${outputDIR}/${PREDIX}.${POP}.bootstrap

echo "generating bootstraps for ${PREFIX}.${POP}"
$MSMCTOOLS/multihetsep_bootstrap.py -n $nbootstraps -s $chunk_size --chunks_per_chromosome $chunk_per_chr --nr_chromosomes $nchr $BS_OUTPUT $BS_INPUT

#echo "multihetsep bootstraps for ${PREFIX}.${POP} done."
