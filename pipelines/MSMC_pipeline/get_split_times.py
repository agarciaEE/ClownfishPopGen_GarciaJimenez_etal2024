#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import argparse

from pandas import NA


def getCCRintersect(df, val, raiseError=True):
    xVec = gen * ((df.left_time_boundary + df.right_time_boundary) / 2) / mu
    yVec = 2.0 * df.lambda_01 / (df.lambda_00 + df.lambda_11)
    yVec = yVec.to_numpy()
    xVec = xVec.to_numpy()
    i = 0
    while yVec[i] < val:
        i += 1
    if not 0 < i <= len(yVec):
        if raiseError:
            raise ValueError("CCR intersection index out of bounds: {}".format(i))
        else:
            return np.nan
    else:
        intersectDistance = (val - yVec[i - 1]) / (yVec[i] - yVec[i - 1])
        return xVec[i - 1] + intersectDistance * (xVec[i] - xVec[i - 1])


def plot_msmc_files(msmc_files, labels, drop_recent=5, drop_past=5, output_file="MSMC_plot.pdf"):
    plt.figure(figsize=(8, 10))

    labLevels = set("-".join(labels).split("-"))
    color_palette = sns.color_palette("colorblind", n_colors=len(labLevels))
    labCols = {label: color_palette[i] for i, label in enumerate(labLevels)}

    color_palette = sns.color_palette("hls", n_colors=len(msmc_files))

    pop1_name = []
    pop2_name = []
    rCCR = []
    for i, msmc_file in enumerate(msmc_files):
        msmc_out = pd.read_csv(msmc_file, sep='\t', header=0)

        msmc_out = msmc_out.iloc[drop_recent:]

        msmc_out = msmc_out.iloc[:-drop_past]

        rCCR.append(getCCRintersect(msmc_out, 0.5, raiseError=False))

        t_years = gen * ((msmc_out.left_time_boundary + msmc_out.right_time_boundary) / 2) / mu
        label0 = labels[i].split("-")[0]
        label1 = labels[i].split("-")[1]

        plt.subplot(211)
        if label0 not in pop1_name + pop2_name:
            plt.semilogx(t_years, (1 / msmc_out.lambda_00) / (2 * mu), drawstyle='steps', color=labCols[label0],
                         label=label0)
        if label1 not in pop1_name + pop2_name:
            plt.semilogx(t_years, (1 / msmc_out.lambda_11) / (2 * mu), drawstyle='steps', color=labCols[label1],
                         label=label1)

        plt.subplot(212)
        relativeCCR = 2.0 * msmc_out.lambda_01 / (msmc_out.lambda_00 + msmc_out.lambda_11)
        plt.semilogx(t_years, relativeCCR, color=color_palette[i], drawstyle='steps')
        if rCCR[i] is not np.nan:
            plt.axvline(x=rCCR[i], color=color_palette[i], linestyle="--", label=f'{labels[i]} rCCR=0.5')

        pop1_name.append(label0)
        pop2_name.append(label1)

    plt.subplot(211)
    plt.xlabel("years ago")
    plt.ylabel("population Sizes")
    plt.legend()

    plt.subplot(212)
    plt.xlabel("years ago")
    plt.ylabel("Relative CCR")
    plt.legend()

    plt.savefig(f"{output_file}.png")

    # Create a DataFrame from lists
    data = {'pop1': pop1_name, 'pop2': pop2_name, 'splitTime': rCCR}
    df = pd.DataFrame(data)

    # Write DataFrame to a TSV file
    df.to_csv(f"{output_file}_splitTimes.tsv", sep='\t', index=False)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("msmc_files", nargs='+', type=str, help="List of MSMC2 result files.")
    parser.add_argument("--mu", type=float, default=4e-08, help="mutation rate.")
    parser.add_argument("--gen", type=int, default=5, help="generation time.")
    parser.add_argument("--out", type=str, default="msmc_result", help="Name of the output plot.")
    parser.add_argument("--labels", type=str, default=None, help="Comma-separated labels for each file.")
    parser.add_argument("--drop_recent", type=int, default=5, help="Drop n first time entries.")
    parser.add_argument("--drop_past", type=int, default=5, help="Drop n last time entries.")
    args = parser.parse_args()

    mu = args.mu
    gen = args.gen

    if len(args.msmc_files) == 1:
        msmc_files = args.msmc_files[0].split("\n")
    else:
        msmc_files = args.msmc_files

    if args.labels is None:
        labels = [file.split("/")[-1].split(".")[0].split("_")[-1] for file in msmc_files]
    else:
        labels = args.labels.split(",")

    plot_msmc_files(msmc_files, labels, drop_recent=args.drop_recent, drop_past=args.drop_past, output_file=args.out)
