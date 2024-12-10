---
title: "Extended Data of Habitat Specialization Impacts Clownfish Demographic Resilience to Pleistocene Sea-Level Fluctuations"
author: Alberto Garcia Jimenez
output: html_document 
date: "2024-12-05"

---

# Project Repository README

This repository contains data, analyses and scripts supporting Garcia Jimenez et al 2025 "Habitat Specialization Impacts Clownfish Demographic Resilience to Pleistocene Sea-Level Fluctuations". 
Below is an overview of the repository structure and its contents.

---

## Raw Data

The repository includes several raw datasets essential for the analyses:

- **Samples information:** Dataset containing all samples information.
- **<species>\_popList:** species specific .txt files with sample ID and population ID information.
- **Genomic data:** Genomic data can be retrieved from the NCBI SRA given the SRR contained in samples_info.csv

---

## Analyses

### **Summary Genomics**

This folder contains genomic summaries for samples and species, along with their corresponding SRR identifiers.

### **MSMC2 (Multiple Sequentially Markovian Coalescent)**
MSMC2 analyses were conducted to infer population size changes over time and coalescent processes.

#### **Components:**
- Individual Bootstraps: Estimates of population variability derived from bootstrap replicates based on individual genome selections.
- Multi-Heterozygosity Bootstraps: Assessments of consistency across genomic regions to validate demographic reconstructions.
- Combined Results: Pairwise cross-coalescence rate (CCR) estimates used to calculate population split times.
- Raw MSMC Datasets: Complete results from individual and multi-heterozygosity bootstrap analyses for all species and populations.
- Filtered MSMC Datasets: Cleaned and curated results from bootstrap analyses for all species and populations.
- Geodesic Distance and Split Times: Dataset summarizing estimated geodesic distances and split times between populations.

### **PopGen**

The `popgen` folder contains produced population genomic data for various clownfish species with pairwise genetic differentiation measures (FST), nucleotide diversity (\(\pi\)) and genetic divergence (\(d_{xy}\)). Below is a detailed description of the contents:

- **`FST_produced_dataset.csv`**  
  Consolidated average FST values for all species and population pairs produced with vcftools.

- **`Pi_produced_dataset.csv`**  
  Summarized nucleotide diversity (\(\pi\)) for all species and populations produced with vcftools.

- **`dxy_produced_dataset.csv`**  
  Combined averaged \(d_{xy}\) values across all species and population pairs produced with `popgenWindows.py` from https://github.com/simonhmartin/genomics_general.

These data sets along with `pca` and `admixture` results can be reproduced following PopGen pipeline in the corresponding `pipelines` folder.

### **EEMS**
EEMS analyses were performed for multiple species to estimate migration rates and population isolation patterns.
This folder structure contains the input data and results of the EEMS (Estimated Effective Migration Surfaces) analyses for different species, demes, and chains. 
Each subfolder corresponds to a specific species and contains the output data for different analysis runs, distinguished by the number of demes and the chain number.

#### **Species Analyzed:**
- AKA
- AKY
- CLK
- CRP
- MEL
- POL
- PRD
- SAN

#### **Configurations:**
- Different deme sizes: 50, 200, and 500.
- Three independent chains for each configuration.

#### **File Descriptions**

##### <species>_all_log_posteriors.png
- A summary image for each species showing log posterior distributions.

##### <species>_eemsFiles/
- Subdirectory containing detailed output files for the specific species.

###### <species>-EEMS-ndemes<deme_count>-chain<chain_number>/
- Each subfolder represents a unique run defined by the number of demes and the chain number. It includes the following files:
  - **demes.txt**: Information on the demes.
  - **edges.txt**: Details about the edges in the model.
  - **eemsrun.txt**: General run information.
  - **ipmap.txt**: IP map data.
  - Various `.txt` files related to model parameters and MCMC chains (e.g., `lastdfpars.txt`, `mcmcmhyper.txt`, etc.).

###### <species>.coord
- Coordinates of the sampled locations.

###### <species>.diffs
- File generated by the bed2diffs program from the EEMS. Contains a pairwise genetic dissimilarity matrix. This matrix quantifies the genetic differences between pairs of demes (spatially defined populations) based on input genotype data, such as a .bed file. The values in the .diffs file represent the genetic distance or dissimilarity, which is used by EEMS to estimate and visualize patterns of gene flow and barriers to migration across a landscape.

###### <species>.order
- Order of data points.

###### <species>.outer
- Outer coordinates indicating the geographical boundary of the data.

###### <species>_50demes-chain<chain_number>.ini
- Configuration files for runs with 50 demes.

###### <species>_200demes-chain<chain_number>.ini
- Configuration files for runs with 200 demes.

###### <species>_500demes-chain<chain_number>.ini
- Configuration files for runs with 500 demes.

###### <species>_area.tiff
- Area visualization map.

###### <species>_excluded_indv.txt
- List of excluded individuals (if any).

###### <species>_wcoords.log
- Plink log file.

###### <species>_wcoords.nosex
- Samples name and sex code.

### **MAPS**

This folder contains input and results of MAPS (Migration And Population Size Surfaces), designed to infer spatial and temporal heterogeneity in population sizes and migration rates across landscapes.

#### **Species Analyzed:**
- AKA
- AKY
- CLK
- CRP
- MEL
- POL
- PRD
- SAN

#### **Configurations:**
- Different deme sizes: 50, 200, and 500.
- Three independent chains for each configuration.

#### **File Descriptions**

##### <species>/
- Subdirectory containing detailed intput files for the specific species.
    - **<species>.coord**: Coordinates of sampled locations.
    - **<species>.demes**: Information about the demes in the model.
    - **<species>.edges**: Details about the edges in the spatial model.
    - **<species>.ipmap**: IP map data.
    - **<species>.maps.0.1_Inf.sims**, **<species>.maps.2_6.sims**, **<species>.maps.6_Inf.sims**: Main input MAPS files for different ibd segment lengths.
    - **<species>.outer**: Outer boundary data.
    - **<species>_200.demes**, **<species>_50.demes**, **<species>_500.demes**: Configuration files specific for different deme sizes.
    - **<species>_ndemes200_params-chain1.ini**, **<species>_ndemes500_params-chain2.ini**, etc.: Configuration files for different chain runs and deme sizes.

##### Additional Files
- **best_deme_bysp.csv**: A CSV file summarizing the best number of demes per species based on log posteriors.

### **Other Analyses**
- **Species Distribution Modeling (SDM):** Predictions of current species ranges performed with ENMTML R package (de Andrade et al. 2020).

---

## Scripts

Scripts used for data processing, analysis, and visualization are included. Detailed descriptions of each script are available in the `scripts/` directory.

### Main Script Files
- **eems_plots.R**: Script for generating plots related to the **EEMS** model output.
- **fst_dxy_pi_plots.R**: Script for plotting genetic differentiation metrics such as **FST**, **DXY**, and **π**.
- **maps_plots.R**: Script for generating plots from **MAPS** data.
- **process_msmc_reconstructions.R**: Script for processing MSMC reconstructions data.
- **msmc_analysis.R**: Script for the analysis of MSMC data and generating plots.
- **pca_admx_plots.R**: Script for performing PCA and ADMIX analysis plots.
- **test_MAPS_recombmaps.R**: Script for visualizing MAPS test on various recombination maps.
- **convert_recomb_map_2PLINK.py**: Python script for converting pyrho recombination maps to PLINK-compatible formats for downstream analyses.
- **get_tmrca.py**: Python script to compute the Time to Most Recent Common Ancestor (TMRCA) from msmc2 cross-coalescence results
- **plot_utils.py**: A collection of Python utility functions from msmc-tools https://github.com/stschiff/msmc-tools
- **custom_functions.R**: An R script containing reusable functions for statistical analysis, data manipulation, and plotting, tailored for this project.

---

## Pipelines

### SNPCall pipeline

The SNPCall pipeline facilitates the identification of Single Nucleotide Polymorphisms (SNPs) from raw sequencing data. It includes preprocessing, mapping reads to a reference genome, variant calling, and filtering steps to generate high-quality SNP datasets. See its own README for more information.

### PopGen pipeline

The PopGen pipeline encompasses analyses of population genetics metrics, including nucleotide diversity (π), genetic differentiation (FST, DXY), and population structure (PCA and ADMIXTURE). It integrates genomic datasets and generates plots to visualize population-level patterns. See its own README for more information.

### EEMS pipeline

The EEMS pipeline supports the application of the EEMS (Estimated Effective Migration Surfaces) model, which maps genetic differentiation across spatial landscapes. This pipeline automates data formatting, model execution, and the creation of migration surface visualizations. See its own README for more information.

### MAPS pipeline

The MAPS pipeline (Migration and Population Surfaces) is based on the framework described in Hussain Al-Asadi et al., 2017. It integrates genetic data to infer and visualize migration rates and population density across spatial landscapes. This pipeline processes input data, estimates migration and population parameters, and generates spatial visualizations to interpret historical and contemporary population dynamics. See its own README for more information.

### MSMC pipeline

The MSMC pipeline processes data for Multiple Sequentially Markovian Coalescent (MSMC) analyses. It includes steps to estimate historical effective population size, preprocess raw output, and perform cross-population coalescence to estimate populations split time. See its own README for more information.

---

## Usage

Follow the instructions below to use this repository effectively:

1. **Navigate to the Analysis Folder**: Start by moving to the appropriate folder corresponding to the analysis you wish to perform. 

2. **Review Analysis Instructions**: Check the subfolders for specific instructions or documentation related to each analysis or refer to the 'Materials and Methods' section on the manuscript for an overview of the analysis protocols and methodology.

3. **Run Scripts for Analysis**:
   - Use the scripts located in the `pipelines` to produce results and `scripts/` folder to replicate and visualize analyses. These scripts are designed to process data and generate results as described in the project.
   - Ensure that you set the correct working directory in the scripts by adjusting the path to match the location of your project files. You can do this using `setwd("/path/to/your/project")` in R or `os.chdir("/path/to/your/project")` in Python.
   - Make sure to install and load any required R or Python packages and external tools as outlined in the **Required Software** section.

4. **Adapt and Customize**:
   - Customize scripts as needed for different datasets or analysis settings. Be mindful of input data formats and required parameters when modifying the code.
   - Refer to the comments within the scripts for additional guidance on how to run specific sections or modify parameters.

5. **Output and Results**:
   - Outputs will be generated in the relevant subfolders and can include figures, log files, or processed data files. Ensure that you have sufficient disk space and the necessary permissions for file output.

By following these steps, you can reproduce analyses and generate results consistent with the project's workflow.

---

## Dependencies

#### **Required Software:**
- R and Python are required to run most scripts.
- Additional packages such as `ggplot2`, `dplyr`, `pandas`, `numpy`, etc., may be required based on the script.
- **ADMIXTURE**: A software tool for estimating individual ancestry proportions.
- **EEMS (Effective Migration Surfaces)**: For spatial modeling of gene flow (available at [EEMS GitHub](https://github.com/dipetkov/eems)).
- **MAPS (Migration And Population size Surfaces)**: For generating spatial maps of genetic differentiation (available at [MAPS GitHub](https://github.com/halasadi/MAPS)).
- **MSMC2 (Multiple Sequential Markovian Coalescent)**: For estimating demographic history from genomic data (available at [MSMC2 GitHub](https://github.com/stschiff/msmc2)).
- **bcftools**: For processing VCF files, variant calling, and filtering.
- **bamtools**: For manipulating and analyzing BAM files.
- **vcftools**: For analyzing VCF files and performing various genomic data operations.
- **plink**: For genome-wide association studies and data manipulation.

#### **Additional Requirements**
- Ensure that all tools and packages are properly installed and accessible in your system's PATH.
- Depending on the analysis, you may need access to a high-performance computing environment or cluster to handle large datasets efficiently.

For installation guides, refer to the documentation of each tool and package. If you encounter any issues, consult the respective support forums or the documentation for troubleshooting.

---

## Contributors to the repository

- **Alberto García Jiménez**: Produced the reproducible repository and wrote the README files.
- **Marion Talbi** and **Milan Malinsky**: Produced the recombination maps shared in the repository.
- Other contributors: Fieldwork expeditions and sample collection by Théo Gaboriau, Lucy M. Fitzgerald, Sara Heim, Anna Marcionetti, Sarah Schmid, Joris Bertrand, Bruno Frédérich, Fabio Cortesi, Marc Kochzius, Ploypallin Rangseethampanya, Phurinat Ruttanachuchote, Wiphawan Aunkhongthong, Sittiporn Pengsakun, Makamas Sutthacheep, Thamasak Yeemin.

## Authors

- **Alberto García Jiménez**: Conceptualized and designed the study, conducted research, performed analyses, interpreted results, and wrote the manuscript.
- **Marion Talbi** and **Milan Malinsky**: Generated the recombination maps.
- **Théo Gaboriau** and **Nicolas Salamin**: Contributed to design the study, data interpretation and manuscript writing.
- Additional authors have contributed to fieldwork and data collection and the revision of the manuscript.

---

## License

This project is licensed under the GNU General Public License (GPL) v3. This means you are free to use, modify, and distribute the code, provided that any derivative works are also licensed under the GPL v3. For more details, please refer to the GPL v3 License.

---
