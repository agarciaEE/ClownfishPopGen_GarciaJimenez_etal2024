---
title: "Migration And Population Surfaces (MAPS) Pipeline"
author: "Alberto Garcia Jimenez"
date: "November 2024"
output: html_document
---

## **Introduction**

This document provides an overview of the pipeline for **Migration And Population Surfaces (MAPS)**, based on [Al-Asadi's MAPS GitHub repository](https://github.com/halasadi/MAPS). The pipeline processes BED files to build IBD data and generate sims input file needed for MAPS. Includes all steps from configuration file creation to MAPS execution.

## **Folder Structure**

### Root Directory

| File/Folder                  | Description                                                                                 |
|------------------------------|---------------------------------------------------------------------------------------------|
| `1.create_all_controlFiles.sh` | Creates control files for MAPS runs containing required parameters. Lower and upper IBD length boundaries as input. |
| `2.make_sims.sh`             | SLURM script to generate `.sims` files using `beagle2sims.R` for each species specifying IBD segment length. |
| `3.run_MAPS.sh`              | Executes MAPS using the generated `.sims` and control files.                                |
| `MAPS_create_control_file.R` | R script to generate MAPS control files based on user-defined parameters.                   |
| `beagle2sims.R`              | R script to transform BEAGLE IBD results into `.sims` files for MAPS.                       |

### ibd_data_pipeline/

#### Scripts

| File/Folder                   | Description                                                   |
|-------------------------------|---------------------------------------------------------------|
| `0.make_recombmap_configFiles.sh` | Creates configuration files for testing recombination maps.            |
| `1.make_configFiles.sh`      | Generates configuration files for each species for the IBD pipeline.           |
| `2.run_ibdpipeline.sh`       | Runs the IBD data pipeline.                                    |


#### Utility Scripts (Copied or adapted from https://github.com/halasadi/MAPS).

| File               | Description                                                   |
|--------------------|---------------------------------------------------------------|
| `interpolate.R`    | Interpolates recombination maps for finer resolution.          |
| `remove-gaps-fibd.py` | Removes gaps from IBD files.                                   |
| `remove_sparse.py` | Filters sparse regions in recombination maps.                  |
| `rename.py`        | Renames sample IDs in IBD files for consistency.               |
| `rewrite_bim.sh`   | Modifies `.bim` files for compatibility.                       |
| `sparse_regions.py`| Identifies and excludes sparse genomic regions.                |

## **Pipeline Overview**

### Step 1: Create Control Files
Run `1.create_all_controlFiles.sh` to generate control files. These are already stored in `ibd_data_pipeline/config_files`.

### Step 2: Generate Simulation Files
Submit the array job with `2.make_sims.sh` to process IBD data and create `.sims` files:
```bash
sbatch 2.make_sims.sh lower_bound upper_bound
```

### Step 3: Run MAPS
Use `3.run_MAPS.sh` to execute MAPS:
```bash
bash 3.run_MAPS.sh params.ini
```

## **Dependencies**

- **R** (with required packages):
  - `optparse`
- **Python** (for preprocessing scripts):
  - `pandas`, `numpy`
- **SLURM** (for batch job submissions)
- **Software**:
  - [BEAGLE](https://faculty.washington.edu/browning/beagle/b4_1.html)
  - [PLINK](https://www.cog-genomics.org/plink/)

## **References**

- MAPS: [Al Asadi's MAPS GitHub Repository](https://github.com/halasadi/MAPS)  
- BEAGLE: [https://faculty.washington.edu/browning/beagle/b4_1.html](https://faculty.washington.edu/browning/beagle/b4_1.html)  
- PLINK: [https://www.cog-genomics.org/plink/](https://www.cog-genomics.org/plink/)

## **Acknowledgments**

This pipeline builds upon scripts and methodologies from the MAPS pipeline by Al Asadi. 

For questions or issues, contact Alberto Garcia Jimenez (agarcia26286@gmail.com).
