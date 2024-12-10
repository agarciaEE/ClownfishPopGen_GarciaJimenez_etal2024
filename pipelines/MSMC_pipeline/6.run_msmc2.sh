#!/usr/bin/bash  --login
# example run: sbatch -J AKA_MDG_msmc2_run 4_run_msmc2.sh AKA MDG AKA_popList.txt all . MSMC/inputPOP "1*2+25*1+1*2+1*3" 3 MSMC/POP_results test 100 false

#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G
#SBATCH --output=%x_%j_%a.out
#SBATCH --error=%x_%j_%a.err
#SBATCH --time=2-00:00:00

##############################
# Argument to pass:
PREFIX=${1-AKA} 		# Species prefix of the file. Default: AKA
POPS=${2-MDG} 			# Comma-separated list of populations to be estimated. Default: MDG
POPFILE=${3-AKA_popList.txt}	# File containing samplesID and popID by columns.
CHR=${4-all}			# Chromosomes to run. f "all", all will be used. Default: "all"
workinDir=${5-.}		# Working directory. Detault: current directory (.)
inputDir=${6-MSMC/inputPOP} 	# Directory where the multihet files by population are stored. Default: MSMC/inputPOPs
PATTERN=${7-"1*2+25*1+1*2+1*3"} # Segment patterning. Default: "1*2+25*1+1*2+1*3"
NR_IND=${8-all} 		# Number of individuals to use from the population. Picks random individuals of the file set. If "all", all infivifuals will be used. Default: "all"
outputDir=${9-MSMC/POP_results}	# Output directory
RUN_NAME=${10-test}		# msmc2 run name
MAXITER=${11-100}		# Number of maximum iterations
CROSS=${12-false}		# Whether running cross-population msmc2
############################################################

IFS=',' read -ra popList <<< "$POPS"
echo "Populations: ${popList[@]}"
npops=${#popList[@]}
echo $npops

if "$CROSS" ; then
	# if one population
	if [ $npops -lt 2 ]; then
	  echo "Not enough populations for a cross-population run."
	  exit 1
	# if more than two populations
        elif [ $npops -gt 2 ]; then
		echo "More than two populations for a cross-population run. Submitting pairwise cross-population runs."
                # compute pairwise comparisons and store as array
		for ((i = 0; i < ${#popList[@]}; i++)); do
			for ((j = i + 1; j < ${#popList[@]}; j++)); do
        			# Access the filenames for pairwise comparison
        			pop1="${popList[i]}"
				pop2="${popList[j]}"
				POPS="$pop1","$pop2"
				echo "Submitting job $SLURM_JOB_NAME_$pop1-$pop2..."
		                # run job for each pairwise comparison
				sbatch -J "$SLURM_JOB_NAME"_"$pop1"-"$pop2" $0 $PREFIX $POPS $POPFILE $CHR $workinDir $inputDir $PATTERN $NR_IND $outputDir $RUN_NAME $MAXITER $CROSS
			done
		done
	# if only two populations
	else
		echo "Running SLURM job on msmc2 cross-population..."

		echo "Using following arguments:

		Species prefix:         "$PREFIX"
		Populations:	        "$POPS"
		Population info file:	"$POPFILE"
		Chromosome(s):          "$CHR"
		Working directory:	"$workinDir"
		Input directory (relative to "$workinDir"):        "$inputDir"
		Segmentation pattern:   "$PATTERN"
		Number of individuals:  "$NR_IND"
		Output directory:	"$outputDir"
		MSMC run name:          "$RUN_NAME"
		Cross-population analysis:	"$CROSS"
		Maximum number of iterations:   "$MAXITER

		echo "
		If settings are wrong, please look at instructions and arguments to pass

		"

		####################
		# LOAD MODULES
		####################

		module load gcc miniconda3
		conda activate MSMC

		####################
		# BEGINNNIN OF SCRIPT
		####################

		cd $workinDir

		mkdir -p "$outputDir"/logs

		if [ "$CHR" == "all" ]; then
		   # Get all chromosome population-level multihet files for the msmc2
		   INPUTFILES=("$inputDir"/"$PREFIX".chr*.multihet.txt)
		   inputINDEX=("$inputDir"/"$PREFIX".chr*.multihet.index)
		   echo "List of chromosomes:"
		   for element in "${INPUTFILES[@]}"; do  echo "$element" | cut -d '.' -f 2; done
		else
		   # Get specified chromosome file
		   INPUTFILES=("$inputDir"/"$PREFIX"."$CHR".multihet.txt)
		   inputINDEX=("$inputDir"/"$PREFIX"."$CHR".multihet.index)
		fi

		if [ $(awk 'NR==1{print NF}' "$POPFILE") -ge 2 ]; then
	                # get samples of each population
			pop1=${popList[0]}
			pop2=${popList[1]}

			ind_pop1=$(cat "$POPFILE" | grep "$pop1" | awk '{print $1}')
                        ind_pop2=$(cat "$POPFILE" | grep "$pop2" | awk '{print $1}')

			# select NR_IND individuals of each population
			if [[ $NR_IND =~ ^[0-9]+$ ]]; then
                                # Create random subset of samples of each population
                                echo "Randomly selecting $NR_IND individuals from each population..."
                                ind_pop1=($(printf "%s\n" "${ind_pop1[@]}" | shuf -n $NR_IND))
                                ind_pop2=($(printf "%s\n" "${ind_pop2[@]}" | shuf -n $NR_IND))
			else
                                echo "Computing pairwise indexes with all individuals of both populations..."
			fi
			# get pairwise indexes
			fileINDEX=${inputINDEX[0]}
			INDEX=""
                        for ind1 in ${ind_pop1[@]}; do
                            ind1_index=$(cat "$fileINDEX" | grep "$ind1" | awk '{print $3}')
                            ind1_pos1=$(echo "$ind1_index" | cut -d',' -f1)
                            ind1_pos2=$(echo "$ind1_index" | cut -d',' -f2)
                            for ind2 in ${ind_pop2[@]}; do
                                ind2_index=$(cat "$fileINDEX" | grep "$ind2" | awk '{print $3}')
                                ind2_pos1=$(echo "$ind2_index" | cut -d',' -f1)
                                ind2_pos2=$(echo "$ind2_index" | cut -d',' -f2)

                                combination="$ind1_pos1-$ind2_pos1,$ind1_pos1-$ind2_pos2,$ind1_pos2-$ind2_pos1,$ind1_pos2-$ind2_pos2"
				INDEX+="$combination,"
                            done
                        done
                        INDEX=$(echo "$INDEX" | sed 's/,$//')
		else
			echo "Error: Info file provided does not contain a second column with population info."
			exit 1
		fi

		# run job with two populations
		outputFile="$outputDir"/"$PREFIX"_$(echo "$POPS" | cut -d',' -f1)-$(echo "$POPS" | cut -d',' -f2)."$CHR"."$RUN_NAME".msmc2

	        echo "Running msmc2 with $SLURM_CPUS_PER_TASK threads and $MAXITER maximum iterations on individuals $ind_pop1 from ${pop1[@]} and ${ind_pop2[@]} from $pop2 populations..."

	        msmc2_Linux -t "$SLURM_CPUS_PER_TASK" -p "$PATTERN" \
	                    -o "$outputFile" -s \
	                    -i "$MAXITER" -I "$INDEX" ${INPUTFILES[@]}

	        echo "Moving log and loop files to $outputDir/logs/."
	        mv "$outputFile".log "$outputFile".loop.txt "$outputDir"/logs/

        fi
else

	############################################################
	## Relounch the script as array with the number of jobs corresponding to the number of populations
	if [[ "$SLURM_ARRAY_TASK_ID" == "" ]] && [ -n "$popList" ]; then
	     # Relaunch this script as an array
	     exec sbatch -J "$SLURM_JOB_NAME" --array=1-$npops $0 $PREFIX $POPS $POPFILE $CHR $workinDir $inputDir $PATTERN $NR_IND $outputDir $RUN_NAME $MAXITER $CROSS
	fi
	############################################################
	POP=${popList[$(echo "$SLURM_ARRAY_TASK_ID - 1" | bc)]}

        echo "Running SLURM job on msmc2 on single population..."

	echo "Using following arguments:

	Species prefix:        	"$PREFIX"
	Population:		"$POP"
	Chromosome(s):		"$CHR"
	Working directory:	"$workinDir"
	Input directory (relative to "$workinDir"):        "$inputDir"
	Segmentation pattern:   "$PATTERN"
	Number of individuals:	"$NR_IND"
	Output directory:	"$outputDir"
	MSMC run name:		"$RUN_NAME"
	Maximum number of iterations:	"$MAXITER

	echo "
	If settings are wrong, please look at instructions and arguments to pass

	"

	####################
	# LOAD MODULES
	####################
	module load gcc miniconda3
	conda activate MSMC

	####################
	# BEGINNNIN OF SCRIPT
	####################

	cd $workinDir

	mkdir -p "$outputDir"/logs

	if [ "$CHR" == "all" ]; then
	   # Get all chromosome population-level multihet files for the msmc2
	   INPUTFILES=("$inputDir"/"$PREFIX"_"$POP".chr*.multihet.txt)
	   echo "List of chromosomes:"
	   for element in "${INPUTFILES[@]}"; do  echo "$element" | cut -d '.' -f 2; done
	else
	   # Get specified chromosome file
	   INPUTFILES=("$inputDir"/"$PREFIX"_"$POP"."$CHR".multihet.txt)
	fi

	# Running msmc2 for a single population or individual
	#####################################################

	# Get first line of genotypes in the multihet file
	file=${INPUTFILES[0]}

	echo "Using file $file to obtain the number of haplotypes..."

	var=$(cat $file | awk '{print $4}' | cut -d ',' -f1 | sed '1q;p')
	n=$(expr ${#var} - 2)
	n_inds=$(echo "${#var}/2" | bc)

	echo "There are $n_inds in the selected population."

	# if all, run all pairs
	if [ "$NR_IND" == "all" ] || [ "$NR_IND" -eq $n_inds ]; then

	 echo "Selecting all $n_inds individual(s)..."
	 # Get all individual genotype indexes
	 INDEX=$(for num in $(seq 0 ${n}); do echo -n "${num},"; done; echo "${n}+1" | bc)

	 # else, run pairs from random selection of the specified number of individuals
	 else
	 # if given maximum number of samples to run (instead of the whole file set)
	 # Create random index picking NR_IND out of the (n) max number of samples
	 echo "Randomly selecting $NR_IND out of $n_inds individual(s)..."

	 INDEX=$(bash MSMC_pipeline/0_pick_random_index.sh "$n" "$NR_IND")

	fi

	outputFile="$outputDir"/"$PREFIX"_"$POP"."$CHR"."$RUN_NAME".msmc2

	echo "Running msmc2 with $SLURM_CPUS_PER_TASK threads and $MAXITER maximum iterations on haplotypes $INDEX..."

	msmc2_Linux -t "$SLURM_CPUS_PER_TASK" -p "$PATTERN" \
		    -o "$outputFile" -s \
	            -i "$MAXITER" -I "$INDEX" ${INPUTFILES[@]}

	echo "Moving log and loop files to $outputDir/logs/."
	mv "$outputFile".log "$outputFile".loop.txt "$outputDir"/logs/

fi

echo "End of the script."
