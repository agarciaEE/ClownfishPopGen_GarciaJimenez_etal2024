#!/usr/bin/bash --login
# example run: sbatch -J PhasingVCF 1B_PhasingVCF.sh ALL_samples.txt . MSMC/vcf FilteredMapping MSMC/phasedVCF whatshap_phased SCAFFOLDS.txt AclarkiiReference/A.clarkii_FinalAssembly.fasta

#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --output=%x_%j_%a.out
#SBATCH --error=%x_%j_%a.err
#SBATCH --time=10:00:00

##############################
# Argument to pass:
sampleList=${1-ALL_samples.txt}        				# File with samples names
workindir=${2-.}                                		# Working directory. Default: Current Directory
inputVCFDir=${3-MSMC/vcf}                 			# Name of the directory with the input VCF files relative to working directory.
inputBAMDir=${4-FilteredMapping}				# Name of the direcotry with the input BAM files relative to working directory.
outputDir=${5-MSMC/phasedVCF}                       		# Name of the ourput directory (to create if non-existent). Default: pahsedVCF
suffix=${6-whatshap_phased}     				# Suffix to add to output file. Default "whatshap_phased".
Chromosomes=${7-SCAFFOLDS.txt}			                # File with the list of chromosomes or scaffolds to run.
RefGen=${8-AclarkiiReference/A.clarkii_FinalAssembly.fasta}	# Reference genome
##############################

echo "Using following arguments:

List of samples:        "$sampleList"
Working directory:	"$workindir"
VCF file directory (if relative path, from "$workindir" :  "$inputVCFDir"
BAM file directory (if relative path, from "$workindir" :  "$inputBAMDir"
Output file directory (if relative path, from "$workindir" :  "$outputDir"
Output file suffix:     "$suffix"
Reference Genome:	"$RefGen"
List of chromosomes:    "$Chromosomes

echo "
If settings are wrong, please look at instructions and arguments to pass

"

############################################################
## Relounch the script as array with the number of jobs corresponding to the number of files
nline=$(cat "$sampleList" | wc -l)

if [[ "$SLURM_ARRAY_TASK_ID" == "" ]]; then
     	# Relaunch this script as an array
	exec sbatch -J "$SLURM_JOB_NAME" --array=1-$nline $0 $sampleList $workindir $inputVCFDir $inputBAMDir $outputDir $suffix $Chromosomes $RefGen
fi

############################################################

####################
# LOAD MODULES
####################

module load gcc/10.4.0 shapeit4/4.1.3 bcftools/1.15.1 miniconda3/4.10.3
conda activate Phasing

####################
# BEGINNNIN OF SCRIPT
####################

# Get chromosome number
IND=$(sed -n ${SLURM_ARRAY_TASK_ID}p "$sampleList")

# Print the name of the sample
echo "Working on:        "$IND

# Change to working directory
cd $workindir

mkdir -p $outputDir/stats/

BAMFILE=$(ls -1 "$inputBAMDir"/"$IND".*.bam | grep -oE '[^/]+$')

# Check if BAM file exists
if [ -z "$BAMFILE" ]; then
	echo "Error: BAM file not found."
    	exit 1
# Check if the variable has more than one element
elif [[ "$BAMFILE" == *" "* ]]; then
	echo "More than one BAM file for the same individual found."
fi

if [[  "$Chromosomes" == "all" ]]; then

	VCFFILE=$(ls -1 "$inputVCFDir"/"$IND".*.vcf.gz | grep -oE '[^/]+$')

	if [ -z "$VCFFILE" ]; then
        	echo "Error: VCF file not found."
        	exit 1
	# Check if the variable has more than one element
	elif [[ "$VCFFILE" == *" "* ]]; then
        	echo "More than one VCF file for the same individual found."
	fi

	# Check if vcf file is properly indexed
	if [ ! -f "$inputVCFDir"/"$VCFFILE".csi ];then
        	echo "Indexing vcf file..."
        	bcftools index "$inputfileDir"/"$VCFFILE"
	fi

	# Phasing...
	echo "Phasing individual $IND..."
        whatshap phase --reference $RefGen --ignore-read-groups -o ${outputDir}/${IND}.${suffix}.vcf ${inputVCFDir}/$VCFFILE ${inputBAMDir}/$BAMFILE

	# Compress with bgzip
	bgzip ${outputDir}/${IND}.${suffix}.vcf

	# Collect statistics...
	echo "Collecting statistics..."
	whatshap stats --tsv=${outputDir}/stats/${IND}.${suffix}.stats.tsv ${outputDir}/${IND}.${suffix}.vcf.gz

	bcftools index ${outputDir}/${IND}.${suffix}.vcf.gz

else
        #VCFFILE=$(ls -1 "$inputVCFDir"/"$IND"*.vcf.gz)
        #if [ -z "$VCFFILE" ]; then
        #        echo "Error: VCF file not found."
        #        exit 1
        # Check if the variable has more than one element
        #elif [[ "$VCFFILE" == *" "* ]]; then
        #        echo "More than one VCF file for the same individual found."
        #fi

	for CHR in $(cat $Chromosomes)
        	do
		echo "Phasing individual $IND chromosome $CHR..."
	        CHR_VCFFILE=$(ls -1 "$inputVCFDir"/"$IND"."$CHR".vcf.gz | grep -oE '[^/]+$')

		# Check if vcf file is found.
        	if [ -z "$CHR_VCFFILE" ]; then
			echo "Error: VCF file not found."
                	exit 1
        	# Check if the variable has more than one element
        	elif [[ "$CHR_VCFFILE" == *" "* ]]; then
                	echo "More than one VCF file for the same individual found."
        	fi

	        # Check if vcf file is properly indexed
        	if [ ! -f "$inputVCFDir"/"$CHR_VCFFILE".csi ];then
                	echo "Indexing vcf file..."
                	bcftools index "$inputfileDir"/"$CHR_VCFFILE"
        	fi

		# Phasing...
                whatshap phase --reference $RefGen --ignore-read-groups -o ${outputDir}/${IND}.${CHR}.${suffix}.vcf ${inputVCFDir}/$CHR_VCFFILE ${inputBAMDir}/$BAMFILE
	        #shapeit4 --input ${inputVCFDir}/$VCFFILE --region "$CHR" --output ${outputDir}/${IND}.${CHR}.${suffix}.vcf --thread "$SLURM_CPUS_PER_TASK"

	        # Compress with bgzip
        	bgzip ${outputDir}/${IND}.${CHR}.${suffix}.vcf

                # Collect statistics...
        	echo "Collecting statistics..."
		whatshap stats --tsv=${outputDir}/stats/${IND}.${CHR}.${suffix}.stats.tsv ${outputDir}/${IND}.${CHR}.${suffix}.vcf.gz

		#bcftools view -Oz -o ${outputDir}/${IND}.${CHR}.${suffix}.vcf.gz ${outputDir}/${IND}.${CHR}.${suffix}.vcf
		bcftools index ${outputDir}/${IND}.${CHR}.${suffix}.vcf.gz
	done
fi

echo "End of script."
