#!/usr/bin/env python3
# example run: python convert_recomb_map_2PLINK.py --intput recomMap.map --output recomMap_plink.map --split_by_chrom False

import argparse
import pandas as pd
import numpy as np
import os

parser = argparse.ArgumentParser()
parser.add_argument("--input", type=str, help="input file.")
parser.add_argument("--output", type=str, help="Name of the output file.")
parser.add_argument("--split_by_chrom", action="store_true", default=True, help="Split output by chromosome.")
args = parser.parse_args()

input_file = args.input
output_file = args.output
split_by_chrom = args.split_by_chrom

# Read the data
data = []
with open(input_file, 'r') as file:
    for line in file:
        chrom, start, end, distance = line.strip().split()
        data.append([chrom, int(start), int(end), float(distance)])

df = pd.DataFrame(data, columns=['chrom', 'start', 'end', 'distance'])

# Renumber chromosomes

## get unique chromosomes
unique_chroms = df['chrom'].unique()
## create dictionary
chrom_dict = {c: i + 1 for i, c in enumerate(unique_chroms)}
## map dictionary into df to create a new column
df['chrom_num'] = df['chrom'].map(chrom_dict)

# Group neighboring positions with the same genetic distance

## Create a boolean mask for group changes
group_change = (df['distance'].diff() != 0) | (df['start'] != df['end'].shift())

## Use cumsum to create group IDs
group_ids = group_change.cumsum()

## Group by these IDs and aggregate
grouped = df.groupby(group_ids).agg({
    'chrom': 'first',
    'chrom_num': 'first',
    'start': 'first',
    'end': 'last',
    'distance': 'first'
})

# Create the final dataframe
result = pd.DataFrame(grouped)

# Create the group identifier
result['group_id'] = result['chrom_num'].astype(str) + '_' + result['end'].astype(str)

# Convert genetic distance (pyrho 'r' = per-base recombination rate) to centimorgans
result['distance_cM'] = result['distance'] * 1e6
result['cumulative_sum'] = result['distance_cM'].cumsum()

# Select and reorder columns
final_df = result[['chrom_num', 'group_id', 'cumulative_sum', 'end']]

# Rename columns
final_df.columns = ['chromosome', 'group_id', 'genetic_distance_cM', 'bp_coordinate']

#print(final_df)

if split_by_chrom:
    # Split by chromosome and write each independently
    for chrom, group in final_df.groupby('chromosome'):
 	# Subtract the minimum value of genetic_distance_cM from the column so each chromosome starts from 0 cM
        group["genetic_distance_cM"] = group["genetic_distance_cM"] - group["genetic_distance_cM"].min()
        
        output_chrom_file = f"{os.path.splitext(output_file)[0]}_chr{chrom}{os.path.splitext(output_file)[1]}"
        group.to_csv(output_chrom_file, sep='\t', index=False, header=False)
    print(f"Files have been split by chromosome and saved with prefix: {os.path.splitext(output_file)[0]}_chr")
else:
    # Save to a single file
    final_df.to_csv(output_file, sep='\t', index=False, header=False)
    print(f"File has been saved as: {output_file}")

