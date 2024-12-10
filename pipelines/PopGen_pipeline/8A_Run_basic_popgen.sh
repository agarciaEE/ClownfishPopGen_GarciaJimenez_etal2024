#!/bin/bash --login
# example run: sbatch -J AKA_popgen_stats 8A_Run_basic_popgen.sh AKA_sample.vcf.gz . MajorMinor PopGenStats/AKA_results AKA popList.txt 50000 10000 100 1000000 all Fst,relatedness,pi,het,TajimaD

#SBATCH --time=24:00:00
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=60G
#SBATCH --output=%x_%j_%a.out
#SBATCH --error=%x_%j_%a.err

# ARGUMENTS:
########################################
vcfFile=${1-AKA_sample.vcf.gz}  	# Vcf file
workindir=${2-.}                	# Working directory
inputfileDir=${3-MajorMinor}    	# Input directory containing the vcf file
outputDir=${4-PopGenStats}		# Output directory
prefix=${5-}	        		# prefix to add to output files
popListFile=${6-popList.txt}            # File containing sample ID in the first columns and pop ID in the 2nd
WinSize=${7-50000}			# Window size for pairwise Fst per site computation. Default: 50000
WinStep=${8-10000}			# Window step for pairwise Fst per site computation. Default: 10000
LDwinmin=${9-100}			# Lower end of the LD window in base pairs. Default: 100
LDwinmax=${10-1000000}			# Upper end of the LD window in base pairs. Default: 1000000 (1Mb)
POP=${11-"all"}				# Population. Leave empty to carry out all populations analyses
METRICS=${12-"all"}			# comma-separated list of metric to estimate (Fst,relatedness,pi,het,TajimaD,ld). Default: all.
##########################################

pop=$(if [ "$POP" = "all" ];then echo ""; else echo "$POP"; fi)
METRICS=$(if [ "$METRICS" = "all" ];then echo "Fst,relatedness,pi,het,TajimaD,ld"; else echo "$METRICS"; fi)

echo "Using following arguments:

VCF file:        		"$vcfFile"
Working directory:		"$workindir"
Input file directory:   	"$inputfileDir"
Output directory:		"$outputDir"
Population info file:   	"$popListFile"
Prefix:        			"$prefix"
Window size:   			"$WinSize"
Window step:  			"$WinStep"
LD min bp between two variants:	"$LDwinmin"
LD max bp between two varians:	"$LDwinmax"
Populations:			"$POP"
Estimating following metrics:	"$METRICS

echo "
If settings are wrong, please look at instructions and arguments to pass

"

####################
#### load module ###
####################

module load gcc/10.4.0 r/4.2.1 vcftools/0.1.14 bcftools/1.15.1 bzip2/1.0.8 samtools/1.15.1
module load gcc/10.4.0 miniconda3/4.10.3 r/4.2.1
conda activate PopGen

######################
#### Begin script ####
######################

cd $workindir

# Make output directory if it doesn't exist
mkdir -p "$outputDir"

bash SNPcall_pipeline/0_split_samples_by_pop.sh "$popListFile" "$prefix"

# Define the directory containing files
directory="$prefix"_popFiles

# Get a list of files in the directory
#files=("$directory"/*)
files=($(ls "$directory"/* | grep "$pop"))

weir_fst_POPS=()

# Iterate through each pairwise comparison
for ((i = 0; i < ${#files[@]}; i++)); do

  file1="${files[i]}"
  pop1=$(echo ${file1%samples.txt} | cut -d '/' -f 2 | sed "s/"$prefix"//g" | tr -d '_')

  weir_fst_POPS+=("--weir-fst-pop $file1 ")
  # if array files contain multiple population files then do population pairwise-fst
  if [ ${#files[@]} -gt 1 ]; then

    if [[ "$METRICS" == *"Fst"* ]]; then

      for ((j = i + 1; j < ${#files[@]}; j++)); do
        # Access the filenames for pairwise comparison
        file2="${files[j]}"

        # Check if the files are different
        if [ "$file1" != "$file2" ]; then

          pop2=$(echo ${file2%samples.txt} | cut -d '/' -f 2 | sed "s/"$prefix"//g" | tr -d '_')

          # Run vcftools --weir-fst for each pairwise combination
    	  echo "Calculating pairwise Fst between $pop1 and $pop2..."
	  vcftools --gzvcf "$inputfileDir"/"$vcfFile" --weir-fst-pop "$file1" --weir-fst-pop "$file2" --out "$outputDir"/"$prefix"_"${pop1}"-"${pop2}"

	  if [ ! -f "$outputDir"/"$prefix"_"${pop1}"-"${pop2}".weir.fst ]; then echo "Pairwise Fst between $pop1 and $pop2 failed."; else echo "Pairwise Fst between $pop1 and $pop2: success"; fi

    	  # Run vcftools --weir-fst for each pairwise combination per site defined by window size and window step
    	  echo "Calculating pairwise Fst per site between $pop1 and $pop2..."
    	  vcftools --gzvcf "$inputfileDir"/"$vcfFile" --weir-fst-pop "$file1" --weir-fst-pop "$file2" --fst-window-size "$WinSize" --fst-window-step "$WinStep" --out "$outputDir"/"$prefix"_"${pop1}"-"${pop2}"

          if [ ! -f "$outputDir"/"$prefix"_"${pop1}"-"${pop2}".windowed.weir.fst ]; then echo "Windowed pairwise Fst between $pop1 and $pop2 failed."; else echo "Windowed pairwise Fst between $pop1 and $pop2: success"; fi

     	fi
      done
    fi
  else

    # if not, do only population metrics
#    echo "Generating independent vcf for population $pop1..."
#    vcftools --gzvcf "$inputfileDir"/"$vcfFile" --keep "$file1" --recode --recode-INFO-all --out "$outputDir"/"${prefix}"_"${pop1}".subset

    if [[ "$METRICS" == *"pi"* ]]; then

      echo "Calculating nucleotide diversity in $pop1 population..."
      # Calculate nucleotide diversity
      vcftools --gzvcf "$inputfileDir"/"$vcfFile" --keep "$file1" --site-pi --out "$outputDir"/"$prefix"_"$pop1"
      #vcftools --vcf "$outputDir"/"${prefix}"_"${pop1}".subset.recode.vcf --site-pi --out "$outputDir"/"$prefix"_"$pop1"

      if [ ! -f "$outputDir"/"$prefix"_"$pop1".pi ]; then echo "Nucleotide diversity in $pop1 failed."; else echo "Population $pop1 nucleotide diversity: success"; fi

    fi

    if [[ "$METRICS" == *"relatedness"* ]]; then

      echo "Calculating relatedness among in dividuals from population $pop1..."
      # Calculate relatedness
      vcftools --gzvcf "$inputfileDir"/"$vcfFile" --keep "$file1"  --relatedness2 --out "$outputDir"/"$prefix"_"$pop1"

      if [ ! -f "$outputDir"/"$prefix"_"$pop1".relatedness2 ]; then echo "Population $pop1 individuals relatedness failed."; else echo "Population $pop1 relatedness: success"; fi

    fi

    if [[ "$METRICS" == *"het"* ]]; then

      echo "Calculating individual heterozygosity in $pop1 population..."
      # Calculate heterozygosity
       vcftools --gzvcf "$inputfileDir"/"$vcfFile" --keep "$file1" --het --out "$outputDir"/"$prefix"_"$pop1"
       #vcftools --vcf "$outputDir"/"${prefix}"_"${pop1}".subset.recode.vcf --het --out "$outputDir"/"$prefix"_"$pop1"

       if [ ! -f "$outputDir"/"$prefix"_"$pop1".het ]; then echo "Individual heterozygosity in $pop1 failed."; else echo "Individual heterozygosity in $pop1: success";  fi

    fi

    if [[ "$METRICS" == *"TajimaD"* ]]; then

       echo "Calculating Tajima's D in $pop1 population..."
       # Calculate Tajima's D
       vcftools --gzvcf "$inputfileDir"/"$vcfFile" --keep "$file1" --TajimaD "$WinStep" --out "$outputDir"/"$prefix"_"$pop1"
       #vcftools --vcf "$outputDir"/"${prefix}"_"${pop1}".subset.recode.vcf --TajimaD "$WinStep" --out "$outputDir"/"$prefix"_"$pop1"

      if [ ! -f "$outputDir"/"$prefix"_"$pop1".Tajima.D ]; then echo "Tajimas D in $pop1 failed."; else echo "Tajimas D in $pop1: success";  fi

    fi

    if [[ "$METRICS" == *"ld"* ]]; then

      echo "Calculating Linkage Disequilibrium (LD) in $pop1 population within windows between $LDwinmin and $LDwinmax..."
      # Calculate Linkage Disequilibrium
      vcftools --gzvcf "$inputfileDir"/"$vcfFile" --keep "$file1" --geno-r2 --ld-window-bp "$LDwinmax" --ld-window-bp-min "$LDwinmin" --out "$outputDir"/"$prefix"_"$pop1"
      #vcftools --vcf "$outputDir"/"${prefix}"_"${pop1}".subset.recode.vcf --geno-r2 --ld-window-bp "$LDwinmax" --ld-window-bp-min "$LDwinmin" --out "$outputDir"/"$prefix"_"$pop1"

      if [ ! -f "$outputDir"/"$prefix"_"$pop1".geno.ld ]; then echo "Linkage Disequilibrium in $pop1 failed.";
      else echo "Linkage Disequilibrium in $pop1: success"; gzip -c "$outputDir"/"$prefix"_"$pop1".geno.ld > "$outputDir"/"$prefix"_"$pop1".geno.ld.gz; fi

    fi

  fi

done

# if array files contain multiple populations files do overall metrics too.
if [ ${#files[@]} -gt 1 ]; then

  if [[ "$METRICS" == *"Fst"* ]]; then

    echo "Calculating Fst per site in all populations..."
    # Calculate pairwise Fst per site
    vcftools --gzvcf "$inputfileDir"/"$vcfFile" "${weir_fst_POPS[@]}" --fst-window-size "$WinSize" --fst-window-step "$WinStep" --out "$outputDir"/"$prefix"_allpops

    if [ ! -f "$outputDir"/"$prefix"_allpops.windowed.weir.fst ]; then echo "All populations windowed pairwise Fst failed."; else echo "Overall Fst: success"; fi

  fi

  if [[ "$METRICS" == *"relatedness"* ]]; then

    echo "Calculating relatedness among in dividuals from all populations..."
    # Calculate relatedness
    vcftools --gzvcf "$inputfileDir"/"$vcfFile" --relatedness2 --out "$outputDir"/"$prefix"

    if [ ! -f "$outputDir"/"$prefix".relatedness2 ]; then echo "All individuals relatedness failed."; else echo "Overall relatedness: success"; fi

  fi

  if [[ "$METRICS" == *"pi"* ]]; then

    echo "Calculating overall nucleotide diversity..."
    # Calculate nucleotide diversity
    vcftools --gzvcf "$inputfileDir"/"$vcfFile" --site-pi --out "$outputDir"/"$prefix"

    if [ ! -f "$outputDir"/"$prefix".pi ]; then echo "Overall nucleotide diversity estimation failed."; else echo "Overall nucleotide diveristy: success";fi

  fi

  if [[ "$METRICS" == *"het"* ]]; then

   echo "Calculating overall heterozygosity..."
    # Calculate heterozygosity
    vcftools --gzvcf "$inputfileDir"/"$vcfFile" --het --out "$outputDir"/"$prefix"

    if [ ! -f "$outputDir"/"$prefix".het ]; then echo "Overall individual heterozygosity failed."; else echo "Overall heterozygosity: success"; fi

  fi

  if [[ "$METRICS" == *"ld"* ]]; then

    echo "Calculating overall Linkage Disequilibrium (LD) within windows between $LDwinmin and $LDwinmax..."
    # Calculate Linkage Disequilibrium
    vcftools  --gzvcf "$inputfileDir"/"$vcfFile" --geno-r2 --ld-window-bp "$LDwinmax"  --ld-window-bp-min "$LDwinmin" --out "$outputDir"/"$prefix"

    if [ ! -f "$outputDir"/"$prefix".geno.ld ]; then echo "Estimation of overall Linkage Disequilibrium failed.";
    else echo "Overall LD: success"; gzip -c "$outputDir"/"$prefix".geno.ld > "$outputDir"/"$prefix".geno.ld.gz; fi

  fi
fi

#echo "Calculating linkage disequilibrium..."
# Calculate linkage disequilibrium decay
#bcftools +fill-tags "$inputfileDir"/"$vcfFile" -O z -o "$outputDir"/"$prefix"_filled.vcf.gz -- -t ^GT
#plink --vcf "$inputfileDir"/"$vcfFile" --double-id --allow-extra-chr --set-missing-var-ids @:# \
#	--maf 0.01 --geno 0.1 --mind 0.5 \
#	--r2 gz --ld-window 1000 --ld-window-kb 100 --ld-window-r2 0 \
#	--make-bed --threads "$SLURM_CPUS_PER_TASK" --out "$outputDir"/"$prefix"_LDdecay

#echo "Cleaning up intermediate files..."
# Cleanup intermediate files
#rm "$outputDir"/"$prefix"_filled.vcf.gz
#rm "$outputDir"/"$prefix"_filled.vcf.gz.tbi

echo "Analysis completed. Results saved in $outputDir with prefix $prefix."
