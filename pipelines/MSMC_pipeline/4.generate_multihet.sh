#!/usr/bin/env bash
# example run1: sbatch -J AKA_multihet 3_generate_multihet.sh ALL_samples.txt MSMC/masks MSMC/ref_masks MSMC/vcf_phased SCAFFOLDS.txt MSMC/inputIND AKA
# example run2: sbatch -J AKA_multihet 3_generate_multihet.sh ALL_popList.txt MSMC/masks MSMC/ref_masks MSMC/vcf_phased SCAFFOLDS.txt MSMC/inputPOP AKA

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
SampleList=${1-AKA_samples.txt} 	# Name of the file with sample names. Default: AKA_samples.txt
MASKDir=${2-MSMC/masks} 		# Directory where individual mask files are stored. Default : MSMC/masks
refmasksDir=${3-MSMC/ref_masks}		# Direcotry where reference genome masks are stored. Default: MSMC/ref_masks
PHASEDVCFDir=${4-MSMC/vcf_phased} 	# Directory where phased vcf files are stored. Default: MSMC/vcf_phased
CHR=${5-SCAFFOLDS.txt} 			# .txt file with a list of chromosomes or a specific chromosome.
OUTPUTDir=${6-MSMC/inputFiles} 		# output directory. Default: MSMC/inputFiles
OUTPrefix=${7-AKA} 			# Output file prefix in case individuals are to be combined. Default: AKA
##############################

if [ -f "$CHR" ] &&  [[ "$SLURM_ARRAY_TASK_ID" == "" ]]; then
	nchr=$(cat $CHR | wc -l)
	echo "Submitting array with $nchr jobs."
	# Relaunch this script as an array
	exec sbatch -J "$SLURM_JOB_NAME" --array=1-$nchr $0 $SampleList $MASKDir $refmasksDir $PHASEDVCFDir $CHR $OUTPUTDir $OUTPrefix
fi

####################
# LOAD MODULES
####################

module load gcc python/3.9.13
MSMCTOOLS=/work/FAC/FBM/DBC/nsalamin/software/nsalamin/clownfishes/agarciaj/msmc-tools/

####################
# BEGINNNIN OF SCRIPT
####################

if [ -n "$SLURM_ARRAY_TASK_ID" ]; then
	CHR=$(sed -n ${SLURM_ARRAY_TASK_ID}p "$CHR")
	echo "Working on chromosome $CHR..."
fi

REFMASK=$(ls "$refmasksDir"/*_"$CHR".mask.bed.gz)

echo "Using following arguments:

List of Samples:        	"$SampleList"
Individual masks directory:     "$MASKDir"
Reference Mask Directory:	"$refmasksDir"
Phased VCF Directory:		"$PHASEDVCFDir"
Chromosome:   			"$CHR"
Reference mask:			"$REFMASK"
Output Dir:			"$OUTPUTDir"
Prefix output file:		"$OUTPrefix

echo "
If settings are wrong, please look at instructions and arguments to pass

"
mkdir -p $OUTPUTDir


if [ $(awk 'NR==1{print NF}' "$SampleList") -ge 2 ]; then

	echo "Grouping individuals using second column information in the file..."

	# Get group IDs
	pops=$(cat "$SampleList" | awk '{print $2}' | sort | uniq)
	for p in $pops; do
		POP_INDS=$(cat "$SampleList" | grep $p | awk '{print $1}')
		nInd=$(echo $POP_INDS | wc -w)
		echo "Subsetting $nInd individuals from population $p..."
		masks=()
		vcfs=()
		missingIND=()
		index=0
		# remove index file if exists
		rm -f "$OUTPUTDir"/"$OUTPrefix"_"$p"."$CHR".multihet.index

		for IND in $POP_INDS; do
			if [ -f "$MASKDir/ind_mask.$IND.$CHR.bed.gz" ] && [ -f  "$PHASEDVCFDir/$IND.$CHR.whatshap_phased.vcf.gz" ]; then
			  masks+=("--mask $MASKDir/ind_mask.$IND.$CHR.bed.gz ")
        		  vcfs+=("$PHASEDVCFDir/$IND.$CHR.whatshap_phased.vcf.gz ")
			  echo "$CHR $IND ${index},$(echo "${index}+1" | bc)" >> "$OUTPUTDir"/"$OUTPrefix"_"$p"."$CHR".multihet.index
			  index=$(echo "$index+2" | bc)
			else
			  echo "ERROR: Mask and/or phased vcf of individual $IND not found. Individual removed"
			  missingIND+=("$IND")
			  nInd=$(echo "$nInd-1" | bc)
			fi
        	done
		# Run multihetsep
		echo "Running multihetsep with $nInd individuals in population $p..."
	        eval "$MSMCTOOLS"/generate_multihetsep.py --chr "$CHR" --mask "$REFMASK" "${masks[@]}" "${vcfs[@]}" > "$OUTPUTDir"/"$OUTPrefix"_"$p"."$CHR".multihet.txt
		echo "Done. MSMC2 input file with name $OUTPrefix_$p.$CHR.multihet.txt generated. Samples indexes can be found in the .index extension of the same filename."
	done
else
        masks=()
        vcfs=()
        missingIND=()
	index=0
        nInd=$(cat "$SampleList" | wc -l)
	echo "Generating multihet file with $nInd individuals..."

	# remove index file if exists
	rm -f "$OUTPUTDir"/"$OUTPrefix"."$CHR".multihet.index

	while IFS= read -r IND; do
		if [ -f "$MASKDir/ind_mask.$IND.$CHR.bed.gz" ] && [ -f  "$PHASEDVCFDir/$IND.$CHR.whatshap_phased.vcf.gz" ]; then
                	masks+=("--mask $MASKDir/ind_mask.$IND.$CHR.bed.gz ")
                	vcfs+=("$PHASEDVCFDir/$IND.$CHR.whatshap_phased.vcf.gz ")
                        echo "$CHR $IND ${index},$(echo "${index}+1" | bc)" >> "$OUTPUTDir"/"$OUTPrefix"."$CHR".multihet.index
                        index=$(echo "$index+2" | bc)
                else
                	echo "ERROR: Mask and/or phased vcf of individual $IND not found. Individual removed"
                	missingIND+=("$IND")
                        nInd=$(echo "$nInd-1" | bc)
                fi
	done < "$SampleList"
        # Run multihetsep
        echo "Running generate_multihetsep.py with $nInd individuals..."
        eval "$MSMCTOOLS"/generate_multihetsep.py --chr "$CHR" --mask "$REFMASK" "${masks[@]}" "${vcfs[@]}" > "$OUTPUTDir"/"$OUTPrefix"."$CHR".multihet.txt
        echo "Done."
fi

####################
# END
####################

echo "Multihet files have been generated succesfully"



