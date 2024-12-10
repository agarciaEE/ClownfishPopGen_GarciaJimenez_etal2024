---
title: "EEMS Input Generation Pipeline"
author: "Alberto Garcia Jimenez"
date: "November 2024"
output: html_document
---

# **Introduction**

This pipeline facilitates the creation of input files required for running EEMS (Estimated Effective Migration Surfaces) analyses. It includes scripts for generating input files, control files, and other necessary components for EEMS.

---

# **Pipeline Overview**

The pipeline performs the following tasks:

Processes a population genetics dataset (.csv) to create BED files and other inputs for EEMS.

Generates control and outer files required for running EEMS simulations.

Supports parallel processing through SLURM job submission (sbatch).

---

# **Workflow**

## 1. Main Script: 1.generate_inputfiles.sh

### Usage

Submit the script using the following command:

```bash
sbatch -J prefix_EEMS_generate_inputfiles 1.generate_inputfiles.sh [dataset] [workindir] [BEDdir] [BEDprefix] [outDir] [chain] [mcmciter] [mcmcburn] [mcmcthin] [regfile] [distfile] [translocate] [overwrite]
```

| **Argument**     | **Description**                                     | **Default**                          |
|-------------------|-----------------------------------------------------|--------------------------------------|
| `dataset`        | Path to the dataset file                            | `popgen_dataset.csv`                 |
| `workindir`      | Working directory                                   | `.`                                  |
| `BEDdir`         | Directory of input BED files                        | `MajorMinor/Pruned`                  |
| `BEDprefix`      | Prefix of BED files                                 | `AKA`                                |
| `outDir`         | Output directory for results                        | `MajorMinor/Results/Pruned/EEMS`     |
| `ndemes`         | Number of demes for EEMS                            | `200`                                |
| `chain`          | MCMC chain ID                                       | `1`                                  |
| `mcmciter`       | Number of MCMC iterations                           | `2000000`                            |
| `mcmcburn`       | MCMC burn-in iterations                             | `1000000`                            |
| `mcmcthin`       | MCMC thinning size                                  | `999`                                |
| `regfile`        | Regions shape file                                  | `EEMS/regfile`                       |
| `distfile`       | Species distribution shape file                     | `EEMS/speciesdist/AKA.shp`           |
| `translocate`    | Whether to translocate longitude coordinates (+360) | `TRUE`                               |
| `overwrite`      | Whether to overwrite existing files                 | `TRUE`                               |


## 2. Supporting R Scripts

### eems_generate_input.R

Processes the dataset to create EEMS-compatible input files (e.g., .diffs and .outer files).

Parameters:

Input dataset, BED directory, and prefix.

Output directory, number of demes, and MCMC settings.

### eems_create_control_file.R

Creates the EEMS control file based on the provided parameters.

### eems_generate_outerfile.R

Generates the "outer" file for EEMS, defining the geographic boundaries.

# **Example Workflow**

1. Prepare Input Files: Organize your dataset (popgen_dataset.csv) and BED files into the appropriate directories.

2. Submit the Script: Use sbatch to run the 1.generate_inputfiles.sh script:

```bash
sbatch -J my_EEMS_job 1.generate_inputfiles.sh dataset.csv . MajorMinor/Pruned AKA Results/EEMS 1 2000000 1000000 999 EEMS/regfile EEMS/speciesdist/AKA.shp TRUE TRUE
```

3. Generated Outputs:

Input files for EEMS (.diffs, .outer, .control).

Logs of the process.

# **Software Requirements**

## Modules:

gcc
r
miniconda3

## Conda Environment: Activate environment before running:

```bash
conda activate EEMS
```

## Dependencies:

PLINK
bed2diffs [https://github.com/dipetkov/eems]
R scripts for input generation

# **References**

EEMS: https://github.com/dipetkov/eems

PLINK: https://www.cog-genomics.org/plink/

For questions or issues, contact Alberto Garcia Jimenez (agarcia26286@gmail.com).
