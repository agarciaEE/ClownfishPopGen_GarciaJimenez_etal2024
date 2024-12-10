---
title: "Demographic Inference Pipeline Using MSMC2"
author: "Alberto Garcia Jimenez"
date: "November 2024"
output: html_document
---

# **Introduction**

This pipeline leverages **MSMC2** to infer past effective population sizes and split times between populations from phased genomic data. The steps cover data preparation, phasing, multi-het file generation, bootstrap analysis, and demographic inference.

The pipeline assumes input data in the form of VCF files and produces output ready for plotting demographic histories.

This pipeline is adapted and extended from the resources and tools available at the following repositories:  
- [MSMC Tools](https://github.com/stschiff/msmc-tools)  
- [MSMC2](https://github.com/stschiff/msmc2)  

These repositories provide the foundational scripts and guidelines for implementing MSMC2 analyses. For more information on MSMC2, refer to https://github.com/stschiff/msmc-tools and https://github.com/stschiff/msmc2 repositories or the **MSMC2 tutorial** included in this repository: `MSMC_Tutorial_v1.2.pdf`.

---

# **Pipeline Overview**

Each script corresponds to a specific step in the MSMC2 pipeline.  
The steps include the following:

1. **Concatenating Phased VCFs**: Combines phased VCF files across chromosomes into a single file.
2. **Phasing and Mask Generation**: Prepares individual-specific masks and phased VCFs.
3. **Global Mask Generation**: Generates a mask file to account for callable regions across the genome.
4. **Multi-het File Generation**: Prepares input files for MSMC2 from phased VCFs.
5. **Bootstrap Analysis**: Generates bootstrap replicates for statistical robustness.
6. **Running MSMC2**: Performs demographic inference using MSMC2.
7. **Cross-Population Analysis**: Combines outputs from different populations for split time estimation.
8. **Split Time Analysis**: Extracts and interprets split times using a custom script.

This pipeline is highly modular and can be adapted to various genomic datasets.

---

# **Installation and Setup**

## Installing MSMC2

Ensure **MSMC2** is installed and accessible. The recommended installation method is via `conda` or direct compilation from the [MSMC2 GitHub repository](https://github.com/stschiff/msmc2).

For example, using `conda`:

```bash
conda create -n MSMC2
conda activate MSMC2
conda install -c bioconda msmc2
conda deactivate
```

## Dependencies

Ensure the following tools are also installed:
- **vcftools** for handling VCF files.
- **bcftools** for processing VCFs.
- **python3** for running scripts like `get_split_times.py`.

---

# **Workflow Steps**

## 0. Accessory Scripts

### Concatenate Phased VCFs

Combine phased VCF files across chromosomes for each individual:

```bash
sbatch 0_concat_phasedVCF.sh samplesList.txt . MajorMinor/Phased 8
```

### Random Index Selection

Randomly select indices for bootstrap replicates:

```bash
bash 0_pick_random_index.sh 20 3
```

## 1. Generate Individual VCFs and Masks

Prepare phased VCFs and individual-specific masks:

```bash
sbatch 1.generate_indVCF_indMASK.sh samplesnames.txt . AclarkiiRef/A.clarkii_FinalAssembly.fasta MSMC FilteredMapping split .BWA.Aclarkii.Sort.Filt_mergedReads.bam 
```

## 2. Phasing VCFs

Phase VCF files to generate input for MSMC2:

```bash
sbatch 2.PhasingVCF.sh samplesnames.txt . MSMC/vcf FilteredMapping MSMC/phasedVCF whatshap_phased SCAFFOLDS.txt AclarkiiReference/A.clarkii_FinalAssembly.fasta
```

## 3. Generate Global Mask

Create a global mask for callable genome regions:

```bash
sbatch 3.generate_globalmask.sh AclarkiiReference/Amphiprion_clarkii.ref.fa Aclarkii_reference
```

## 4. Generate Multi-Hetsep Files

Prepare multi-hetsep files required for MSMC2:

```bash
sbatch 4.generate_multihet.sh samplesnames.txt MSMC/masks MSMC/ref_masks MSMC/vcf_phased SCAFFOLDS.txt MSMC/inputIND prefix
```

## 5. Generate Bootstrap Multi-Hetsep Files

Create bootstrap replicates for MSMC2:

```bash
sbatch 5.generate_multihetsep_bootstrap.sh prefix . MSMC/inputPOP/ MSMC/inputPOP/bootstraps/ popList.txt 50 2000000 10
```

## 6. Run MSMC2

Run MSMC2 to infer demographic history:

```bash
sbatch 6.run_msmc2.sh prefix POP popList.txt all . MSMC/inputPOP "1*2+25*1+1*2+1*3" 3 MSMC/POP_results runname 100 false
```

## 7. Combine Outputs for Cross-Population Analysis

Merge MSMC2 outputs from multiple populations:

```bash
sbatch 7.combine_CROSS_POP.sh Crossfiles.txt . MSMC/CROSS_results MSMC/POP_results MSMC/Combined_results
```

## 8. Extract and Analyze Split Times

Use a Python script to extract split times:

```bash
sbatch 8_get_split_times.sh PREFIX . MSMC/Combined_results MSMC/Combined_results PREFIX_MSMC 5 2 4e-8 5
```

---

# **Reference**

For detailed instructions, refer to the tutorial included in this repository: `MSMC_Tutorial_v1.2.pdf`.

If you use this pipeline in your research, please cite:

- **Schiffels & Wang (2020)** for MSMC2.
- This repository and associated scripts.

For questions or issues, contact Alberto Garcia Jimenez (agarcia26286@gmail.com).

