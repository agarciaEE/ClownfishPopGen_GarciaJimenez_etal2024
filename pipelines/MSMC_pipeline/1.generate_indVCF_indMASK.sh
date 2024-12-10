#!/usr/bin/bash
# example run: sbatch -J genVCF_MASK 1A.generate_indVCF_indMASK.sh samplesnames.txt . AclarkiiRef/A.clarkii_FinalAssembly.fasta MSMC FilteredMapping split .BWA.Aclarkii.Sort.Filt_mergedReads.bam

#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --output=%x_%a.out
#SBATCH --error=%x_%a.err
#SBATCH --time=1-00:00:00

##############################
# Argument to pass:
SampleList=${1-samplesnames.txt} 		# Name of the file with sample names. Default: samplesnames.txt
workindir=${2-.}				# Working directory
GENOME=${3-AclarkiiRef/A.clarkii_FinalAssembly.fasta}	# Directory where bam files are stored. Default : FilteredMapping
OUTDIR=${4-MSMC} 				# Directory for storing output files
BAMDIR=${5-FilteredMapping} 			# Directory with BAM files5
CHR=${6-all}		 			# File with chromosome info
BAMSUFFIX=${7-.BWA.Aclarkii.Sort.Filt_mergedReads.bam}	# BAM file suffix
############################################################

echo "Using following arguments:

List of Samples:        				"$SampleList"
Working directory:					"$workindir"
Output directory (relative to working directory):	"$OUTDIR"
Reference genome:       				"$GENOME"
BAM files dirctory:     				"$BAMDIR"
BAM file suffix:		                  	"$BAMSUFFIX"
Chromosome:                                             "$CHR

echo "
If settings are wrong, please look at instructions and arguments to pass

"

############################################################
## Relounch the script as array with the number of jobs corresponding to the number of files
nline=$(wc -l $SampleList | awk '{print $1}')

if [[ "$SLURM_ARRAY_TASK_ID" == "" ]]; then
     	# Relaunch this script as an array
	exec sbatch -J "$SLURM_JOB_NAME" --array=1-$nline $0 $SampleList $workindir $GENOME $OUTDIR $BAMDIR $CHR $BAMSUFFIX
fi

############################################################

####################
# LOAD MODULES
####################

module load gcc/10.4.0 bwa/0.7.17 python/3.9.13 samtools/1.15.1 picard/2.26.2 bamtools/2.5.2 r/4.2.1 vcftools/0.1.14 bcftools/1.15.1

MSMCTOOLS=/work/FAC/FBM/DBC/nsalamin/software/nsalamin/clownfishes/agarciaj/msmc-tools

####################
# BEGINNNIN OF SCRIPT
####################

cd $workindir

# make directories for intermediate files-- will fail if these don't exist
mkdir -p ${OUTDIR}/vcf
mkdir -p ${OUTDIR}/masks
mkdir -p ${OUTDIR}/depthInfo

IND=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SampleList)
BAMFILE=${BAMDIR}/${IND}${BAMSUFFIX}

echo "Current script: 1A.generate_indVCF_indMASK.sh"
echo "Working on sample: $IND"
echo "Bamfile: $BAMFILE"

if [ -f "${BAMFILE}.bai" ]
        then
                echo "Bamfile index already exists, moving on!"
        else
                echo "Bamfile index does not exist, creating index"
                samtools index $BAMFILE ${BAMFILE}.bai
fi

if [[ $CHR == "all" ]]
then
        ### Calculate mean coverage (to be used as input for bamCaller.py):
        MEANCOV=$(samtools depth $BAMFILE | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.') # calculate mean coverage
        echo ${IND} $MEANCOV >> ${OUTDIR}/depthInfo/${IND}_DepthInfo.txt # save mean coverage in separate file
        echo "Mean coverage for this sample is: $MEANCOV"

        ### Generate a single-sample VCF and a mask-fill
        MASK_IND=${OUTDIR}/masks/ind_mask.${IND}.bed.gz # Individual mask file to be created
        VCF=${OUTDIR}/vcf/${IND}.vcf # VCF file to be created

    	echo "Generating VCF and mask for the whole genome of sample $IND..."
        bcftools mpileup --threads $SLURM_NTASKS -q 20 -Q 20 -C 50 -f $GENOME $BAMFILE | bcftools call --threads $SLURM_NTASKS -c -V indels | $MSMCTOOLS/bamCaller.py --minMapQ 20.0 --minConsQ 20.0 $MEANCOV $MASK_IND | bcftools view -Oz -o $VCF.gz
        bcftools index $VCF.gz
	echo "Filtered VCF and mask created for whole genome of sample ${IND}."

elif [[ $CHR == "split" ]]
then
    	Chromosomes=$(samtools idxstats $BAMFILE | cut -f 1 | grep "chr")
        for CHR in $Chromosomes
        do
		echo "Working on chromosome $CHR..."
                ### Calculate mean coverage (to be used as input for bamCaller.py):
                MEANCOV=$(samtools depth -r $CHR $BAMFILE | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.') # calculate mean coverage
                echo ${IND}.${CHR} $MEANCOV >> ${OUTDIR}/depthInfo/${IND}_DepthInfo.txt # save mean coverage in separate file
                echo "Mean coverage for this sample in chromosome ${CHR} is: $MEANCOV"

	        ### Generate a single-sample VCF and a mask-fill
        	MASK_IND=${OUTDIR}/masks/ind_mask.${IND}.${CHR}.bed.gz # Individual mask file to be created
        	VCF=${OUTDIR}/vcf/${IND}.${CHR}.vcf # VCF file to be created

                echo "Generating VCF and mask in chromosome $CHR for sample $IND..."
                bcftools mpileup --threads $SLURM_NTASKS -q 20 -Q 20 -C 50 -r $CHR -f $GENOME $BAMFILE | bcftools call --threads $SLURM_NTASKS -c -V indels | $MSMCTOOLS/bamCaller.py --minMapQ 20.0 --minConsQ 20.0 $MEANCOV $MASK_IND | bcftools view -Oz -o $VCF.gz
        	bcftools index $VCF.gz
        done
	echo "Filtered VCF and mask created for all chromosomes of sample ${IND}."
else
        echo "Working on chromosome $CHR..."
        ### Calculate mean coverage (to be used as input for bamCaller.py):
        MEANCOV=$(samtools depth -r $CHR $BAMFILE | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.') # calculate mean coverage
        echo ${IND}.${CHR} $MEANCOV >> ${OUTDIR}/depthInfo/${IND}_DepthInfo.txt # save mean coverage in separate file
        echo "Mean coverage for this sample in chromosome ${CHR} is: $MEANCOV"

        ### Generate a single-sample VCF and a mask-fill
        MASK_IND=${OUTDIR}/masks/ind_mask.${IND}.${CHR}.bed.gz # Individual mask file to be created
        VCF=${OUTDIR}/vcf/${IND}.${CHR}.vcf # VCF file to be created

        echo "Generating VCF and mask in chromosome $CHR for sample $IND..."
        bcftools mpileup --threads $SLURM_NTASKS -q 20 -Q 20 -C 50 -r $CHR -f $GENOME $BAMFILE | bcftools call --threads $SLURM_NTASKS -c -V indels | $MSMCTOOLS/bamCaller.py --minMapQ 20.0 --minConsQ 20.0 $MEANCOV $MASK_IND | bcftools view -Oz -o $VCF.gz
        bcftools index $VCF.gz

	echo "Filtered VCF and mask created for chromosome $CHR in sample ${IND}."
fi

        # -q = min. mapping qual; -Q = min. base qual; -C = coefficient for downgrading mapping qual for reads w/ excessive mismatches; -u = generate uncompressed VCF/BCF; -r = only specified region; -f = fasta.
        # bcftools: "call" will call SNPs/indels; "-V indels" skips indels; -c = consensus caller.

echo "Done with script."
