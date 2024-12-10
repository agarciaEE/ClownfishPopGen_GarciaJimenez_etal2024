#!/usr/bin/env python3

import pandas as pd
import argparse
import plot_utils

parser = argparse.ArgumentParser()
parser.add_argument("--input", type=str, help="MSMC2 result file.")
parser.add_argument("--mu", type=float, default=4e-08, help="mutation rate.")
parser.add_argument("--gen", type=int, default=5, help="generation time.")
parser.add_argument("--lambda_index", type=int, default=0, help="Index of lambda column.")
parser.add_argument("--resolution", type=int, default=10, help="Sets the number of time points in each time-segment, "
                                                               "default=10.")
parser.add_argument("--cdf", type=bool, default=False, help="Cumulative density distribution.")
parser.add_argument("--outfile", type=str, default="msmc_result", help="Name of the output plot.")
args = parser.parse_args()

mu = args.mu
gen = args.gen
lambda_index = args.lambda_index
resolution = args.resolution
cdf = args.cdf
output_file = args.outfile
filename = args.input

x, y = plot_utils.tmrcaDistribution(filename, resolution=resolution, lambda_index=lambda_index, mu=mu, gen=gen, cdf=cdf)

# Create a DataFrame from lists
data = {'time': x, 'probability': y}
df = pd.DataFrame(data)

# Write DataFrame to a TSV file
df.to_csv(f"{output_file}_tmrca_distribtion.tsv", sep='\t', index=False)
