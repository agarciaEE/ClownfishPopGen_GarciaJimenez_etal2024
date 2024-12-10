#!/bin/bash

## Arguments
input_file=$1	 # input file
output_prefix=$2 # shared prefix for output files

# Check if the input file exists
if [ ! -e "$input_file" ]; then
  echo "Input file '$input_file' does not exist."
  exit 1
fi

# Create a directory to store the output files
output_dir=$output_prefix"_popFiles"
# remove dir if exists
rm -rf "$output_dir"
mkdir -p "$output_dir"

# Read the input file and split based on the second column
while read -r line; do
  # Extract the first and second columns using 'awk'
  sampleID=$(echo "$line" | awk '{print $1}')
  pop=$(echo "$line" | awk '{print $2}')

  # Create a separate text file for each unique value in the second column
  output_file="${output_dir}/${output_prefix}_${pop}_samples.txt"
  echo "${sampleID}" >> "${output_file}"
done < "$input_file"

echo "Files have been split based on each population and saved in the '$output_dir' directory."

