#!/usr/bin/env bash
#example run: bash 0_extract_samples_names.sh . "" "" "[A-Z]*_[0-9]*" samplesnames.txt

INPUTDIR=${1-.} # input directory where to extract samples names
PREFIX=${2-""} # prefix to remove from file names (optional)
SUFFIX=${3-""} # suffix to remove from file names (optional)
PATTERN=${4-"[A-Z]*_[0-9]*"} # pattern to extract from file names
OUTPUTFILE=${5-samplesnames.txt} # name of generated output file

list=$(ls -p $INPUTDIR | grep -v /)

declare -a SAMPLES
for file in $list
do
	# remove prefix and suffix
	string="${file#$PREFIX}"
	string="${string%$SUFFIX}"
	# extract pattern
	string="$(echo $string | sed 's/\('$PATTERN'\).*/\1/')"
	SAMPLES+=("$string")
done

SAMPLES=$(printf "%s\n" "${SAMPLES[@]}" | sort -u)

printf "%s\n" "${SAMPLES[@]}" | sort -u  > $OUTPUTFILE
