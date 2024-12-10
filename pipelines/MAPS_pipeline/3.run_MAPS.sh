#!/bin/bash -login
# examples run: sbatch run_MAPS.sh params.ini

#SBATCH --job-name=MAPS
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mem=100G
#SBATCH --time=2-00:00:00

######################
### ARGUMENTS
######################
paramsFile=${1-params.ini}
######################

######################
### SOFTWARE
######################
maps=/work/FAC/FBM/DBC/nsalamin/software/nsalamin/clownfishes/agarciaj/MAPS/src/runeems2

######################
### LOAD DEPENDNENCIES
######################
module load miniconda3

##############################
### ACTIVATE CONDA ENVIRONMENT
##############################
conda activate MAPS

#### EXPORT BOOST LIBRARY PATH
export LD_LIBRARY_PATH=/work/FAC/FBM/DBC/nsalamin/software/nsalamin/clownfishes/agarciaj/Conda/envs/MAPS/lib


# Extract parameters from the file
ndemes=$(grep "nDemes" "$paramsFile" | awk '{print $3}')
full_path=$(grep "datapath" "$paramsFile" | awk '{print $3}')
lwrbd=$(grep "lowerBound" "$paramsFile" | awk '{print $3}')
uprbd=$(grep "upperBound" "$paramsFile" | awk '{print $3}')
outDir=$(grep "mcmcpath" "$paramsFile" | awk '{print $3}')

# Extract the path (everything before the last '/')
path="${full_path%/*}/"
# Extract the prefix (the last part after the last '/')
prefix="${full_path##*/}"

echo "Path: $path"
echo "Prefix: $prefix"
echo "Number of demes: $ndemes"
echo "Lower boundary: $lwrbd"
echo "Upper boundary: $uprbd"
echo "Output Directory: $outDir"

# Create the output directory
echo "Creating output directory..."
mkdir -p "${path}${outDir}"

# Renaming .edges and .demes files
edges_file="${full_path}_${ndemes}.edges"
demes_file="${full_path}_${ndemes}.demes"
ipmap_file="${full_path}_${ndemes}.ipmap"

if [[ -f "$edges_file" ]]; then
    echo "Copying $edges_file to ${full_path}.edges"
    cp "$edges_file" "${full_path}.edges"
else
    echo "Warning: $edges_file does not exist."
fi

if [[ -f "$demes_file" ]]; then
    echo "Copying $demes_file to ${full_path}.demes"
    cp "$demes_file" "${full_path}.demes"
else
    echo "Warning: $demes_file does not exist."
fi

if [[ -f "$ipmap_file" ]]; then
    echo "Copying $ipmap_file to ${full_path}.ipmap"
    cp "$ipmap_file" "${full_path}.ipmap"
else
    echo "Warning: $ipmap_file does not exist."
fi

# Duplicate coords in corrd file if not done
coord_file="${full_path}.coord"

nind=$(cat ${prefix}_MAPSFiles/${prefix}_ndemes${ndemes}_params-chain1.ini | grep nIndiv | awk '{print $3}')
ncoord=$(cat  ${coord_file} | wc -l)

if [[ $ncoord -eq $((nind / 2)) ]]; then
	echo "number of coords (${ncoord}) equal to half of number of individuals (${nind}). Duplicating coords file to match number of individuals."
	cat ${coord_file} > temp_file
	cat temp_file >> ${coord_file}
	rm temp_file
fi

if [[ $(echo ${outDir} | grep pyrho | wc -l) -eq 1 ]]; then

RecombMap=$(echo "$(basename "$outDir")" | rev | cut -d'-' -f1-2 | rev)

echo "Using ${RecombMap}..."

sims_file="${path}ibd_pipeline/RecombMapsTests/${RecombMap}/${prefix}.maps.${lwrbd}_${uprbd}.sims"

else

# Renaming .sims file
sims_file="${path}ibd_pipeline/${prefix}.maps.${lwrbd}_${uprbd}.sims"

fi

if [[ -f "$sims_file" ]]; then
    echo "Copying $sims_file to ${path}${prefix}.sims"
    cp "$sims_file" "${path}${prefix}.sims"
else
    echo "Warning: $sims_file does not exist."
fi

echo "Running MAPS..."
$maps --params "$paramsFile"

#### DONE
echo "Job finished."


