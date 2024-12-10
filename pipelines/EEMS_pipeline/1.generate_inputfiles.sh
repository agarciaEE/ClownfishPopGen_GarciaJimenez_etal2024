#!/usr/bin/bash --login
# example run: sbatch -J prefix_EEMS_generate_inputfiles 1.generate_inputfiles.sh dataset workindir BEDdir BEDprefix outDir outDir chain mcmciter mcmcburn mcmcthin regfile distfile TRUE TRUE

#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# ARGUMENTS:
###############
dataset=${1-popgen_dataset.csv}                 # complete dataset
workindir=${2-.}                                # Working directory
BEDdir=${3-MajorMinor/Pruned}                   # Input BED file directory
BEDprefix=${4-AKA}				# BED file prefix
outDir=${5-MajorMinor/Results/Pruned/EEMS}      # Output file directory
ndemes=${6-200}					# Number of demes. Default 200
chain=${7-1}					# MCMC chain number. Default 1
mcmciter=${8-2000000}				# MCMC number of iterations. Default 2000000
mcmcburn=${9-1000000}				# MCMC brun-in. Default 1000000
mcmcthin=${10-999}				# MCMC thining size. Default 999
regfile=${11-EEMS/regfile}                      # Regions shape file
distfile=${12-EEMS/speciesdist/AKA.shp}         # species distribution shape file
translocate=${13-TRUE}                          # Wheter translocate longitud coordintates (+360). Default TRUE
overwrite=${14-TRUE}                            # Whether to overwrite input files. Default TRUE
indcol=${15-1}					# Column number containing samples ID in the dataset. Default 1
spscol=${16-2}					# Column number containing species names in the dataset. Default 2
abrcol=${17-3}					# Column number containing species abbreviations in the dataset. Default 3
popcol=${18-9}					# Column number containing population ID in the dataset. Default 9
loncol=${19-10}					# Column number containing longitud coordinates in the dataset. Default 10
latcol=${20-11}					# Column number containing latitud coordinates in the dataset. Default 11
###############

####################
#### load module ###
####################

module load gcc/10.4.0 r/4.2.1 miniconda3

conda activate FEEMS

# path to programs (check if this changes and modify script accordinly)
#######################################################################
plinkPATH=/work/FAC/FBM/DBC/nsalamin/software/nsalamin/clownfishes/agarciaj/plink
diffsPATH=/work/FAC/FBM/DBC/nsalamin/software/nsalamin/clownfishes/agarciaj/eems/bed2diffs/src-wout-openmp/bed2diffs_v1
makeouterPATH=/work/FAC/FBM/DBC/nsalamin/clownfish/agarciaj/EEMS_pipeline/eems_generate_outerfile.R

cd "$workindir"

echo "Using following parameters :

Datasset :                          		"$dataset"
Working directory :                    		"$workindir"
Input file directory (relative to $workindir) : "$BEDdir"
BED file prefix :				"$BEDprefix"
Output file directory :                    	"$outDir"
Number of demes :				"$ndemes"
MCMC chain ID :					"$chain"
Number of MCMC iterations :			"$mcmciter"
Number of MCMC iteration to burn :		"$mcmcburn"
MCMC thining size :				"$mcmcthin"
Translocate longitud coordinates : 		"$translocate"
Regions file :					"$regfile"
Species distribution file :			"$distfile"
Overwrite :					"$overwrite

echo "
If settings are wrong, please look at instructions and arguments to pass

"

mkdir -p "$outDir"

echo "Generating input files from dataset..."

Rscript EEMS_pipeline/eems_generate_input.R --infofile "$dataset" --beddir "$BEDdir" --bedprefix "$BEDprefix" --outputdir "$outDir" --overwrite "$overwrite" \
				--ndemes "$ndemes" --chain "$chain" --mcmciter "$mcmciter" --mcmcburn "$mcmcburn" --mcmcthin "$mcmcthin" \
				--indcol "$indcol" --spscol "$spscol" --abrcol "$abrcol" --popcol "$popcol" --loncol "$loncol" --latcol "$latcol" --translon "$translocate" \
				--regshfile "$regfile" --distshdir "$distfile" --plinkpath "$plinkPATH" --bed2diffspath "$diffsPATH" --makeouterpath "$makeouterPATH"


echo "Script done."
