#!/bin/bash

# Arguments
prefix=$1
workindir=${2-.}
inputfileDir=${3-Admixture}

cd "$workindir"

if [ -f "$inputfileDir"/"$prefix"_CV.error ]; then
	# remove outptut filename if exists
	rm "$inputfileDir"/"$prefix"_CV.error
fi

for i in $(ls "$inputfileDir"/"$prefix"*.out)
do
awk '/CV/ {print $3,$4}' "$i" | cut -c 4,7-20 >> "$inputfileDir"/"$prefix"_CV.error
done

echo "CV errors of $prefix Admixture files have been saved in "$prefix"_CV.error."

