---
title: "Population Genomics Pipeline"
author: "Alberto Garcia Jimenez"
date: "November 2024"
output: html_document
---

# **Introduction**

This pipeline facilitates population genomics analyses, including linkage disequilibrium (LD) pruning, principal component analysis (PCA), admixture analysis, and basic population genetic statistics. The pipeline is designed to work with genotype data in **PLINK** format and assumes the necessary input files are prepared beforehand.

---

## **Pipeline Overview**

The pipeline consists of five main steps, each represented by a script:

1. **Pruning for Linkage Disequilibrium (LD):** Filters SNPs to reduce redundancy caused by LD.
2. **Principal Component Analysis (PCA):** Performs PCA to explore population structure.
3. **Admixture Analysis:** Estimates individual ancestry proportions using ADMIXTURE.
4. **Cross-Validation for Admixture:** Calculates cross-validation errors to determine the optimal number of ancestral populations (K).
5. **Basic Population Genetics Statistics:** Computes summary statistics such as heterozygosity, allele frequencies, and F-statistics.

---

## **Pipeline Components**

| **Script**          | **Description**                                                          |
|----------------------|--------------------------------------------------------------------------|
| `6A_PruneLD.sh`     | Prunes SNPs for LD using PLINK.                                           |
| `6B_PlinkPCA.sh`    | Runs PCA on pruned SNPs using PLINK.                                      |
| `7A_Admixture.sh`   | Executes ADMIXTURE analysis to estimate ancestry proportions.             |
| `7B_getCVerrors.sh` | Computes cross-validation errors for ADMIXTURE runs.                     |
| `8A_Run_basic_popgen.sh` | Performs basic population genetic analyses using PLINK and other tools. |

---

## **Installation and Setup**

### **Software Requirements**

The pipeline requires the following software to be installed and accessible in the system's `PATH`:

- **PLINK**: For genotype manipulation and analysis.
- **ADMIXTURE**: For ancestry inference.

## **Analyses**

Below are example usage commands for each script. Adjust the arguments as needed for your dataset and analysis.

---

### **1. LD Pruning**

```bash
sbatch 6A_PruneLD.sh file.vcf.gz . inputFolder outputFolder prefix WindowSize SlidingWindow r2threshold
```

- **Input File:** `input.vcf.gz`
- **Output Files:** LD-pruned genotype files (`output_prefix.*`)

---

### **2. PCA Analysis**

```bash
sbatch 6B_PlinkPCA.sh file.vcf.gz . inputFolder outputFolder prefix file.prune.in popList.txt
```

- **Input Files:** `file.vcf.gz`, `file.prune.in` and `poplist.txt`
- **Output Files:** PCA eigenvalues and eigenvectors (`output_prefix.eigenval`, `output_prefix.eigenvec`)

---

### **3. Admixture Analysis**

```bash
sbatch 7A_Admixture.sh input.bed . inputFolder outputFolder prefix "1-K" CV nBootstraps
```

- **Input Files:** PLINK binary files (`input.bed`)
- **Output Files:** ADMIXTURE results (`output_dir/*`)

---

### **4. Cross-Validation Errors**

```bash
bash 7B_getCVerrors.sh prefix workingDir input_dir
```

- **Input Files:** Same as Admixture Analysis.
- **Output Files:** Cross-validation error file (`input_dir/cv_errors.txt`)

---

### **5. Basic Population Genetics**

```bash
sbatch 8A_Run_basic_popgen.sh input.vcf.gz . inputDir outputDir prefix popList.txt 50000 10000 100 1000000 all Fst,pi
```

- **Input Files:** `input.vcf.gz`, `poplist.txt`
- **Output Files:** Summary statistics (`output_dir/*`)

Dxy analysis, a measure of genetic divergence between populations, can be performed using the tools provided in the genomics_general repository by Simon H. Martin. This toolkit is specifically designed for genome-wide population genetic analyses, including calculations of Dxy and other summary statistics.

For detailed instructions and examples, please refer to the repository:
genomics_general by Simon H. Martin

#### Key Steps to Perform Dxy Analysis:

**Prepare Input Files:**

Ensure you have genotype data in the appropriate format. The toolkit typically expects .geno or .vcf files containing SNP information for your populations.

**Install Required Tools:**

Clone the repository and install any dependencies, such as Python libraries (e.g., NumPy, SciPy). Use the following commands to set up the repository:

```bash
git clone https://github.com/simonhmartin/genomics_general.git
cd genomics_general
```

**Run Dxy Analysis:**

Use the provided Python scripts to calculate Dxy across genomic windows. For example:

```bash
python popgenWindows.py -g input.geno -o output.txt -w 50000 -s 10000 -p Pop1 -p Pop2 --stats Dxy
```

| **Argument** | **Description**                                                | **Example**        |
|--------------|----------------------------------------------------------------|--------------------|
| `-g`         | Input genotype file.                                           | `input_genotype.vcf` |
| `-o`         | Output file for Dxy results.                                  | `dxy_results.csv`  |
| `-w`         | Window size for analysis (e.g., 50,000 bp).                    | `50000`            |
| `-s`         | Step size for sliding windows.                                | `1000`             |
| `-p`         | Population labels (must match those in the input file).        | `pop1,pop2`        |
| `--stats`    | Specify the statistic to calculate (e.g., Dxy).                | `Dxy`              |

## **Visualize Results**

You can visualize the results of each analysis using the R scripts provided in the scripts folder. These scripts are designed to help you generate plots and figures of the manuscript.

---

## **Input File Format**

The pipeline expects genotype data in **VCF** and **PLINK** format, consisting of the following files:

- `.vcf.gz`: Compressed VCF file containing the variant data, including information on SNPs, genotypes, and other genomic annotations.
- `.bed`: Binary genotype file.
- `.bim`: SNP information file.
- `.fam`: Sample information file.

---

## **Output Files**

The pipeline generates a variety of output files, including:

- **PCA Results:** Eigenvalues and eigenvectors summarizing population structure.
- **Admixture Results:** Ancestry proportions and cross-validation errors.
- **Summary Statistics:** Population genetics metrics, such as nucleotide diversity (π) in .sites.pi files, FST in .weir.fst files, and Dxy in .csv files, providing insights into genetic variation and differentiation among populations.

---

## **Contact and Support**

For questions or issues, please contact:

- **Name:** Alberto Garcia Jimenez
- **Email:** agarcia26286@gmail.com

---

