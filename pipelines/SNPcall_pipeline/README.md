---
title: "SNP Calling Pipeline Using ATLAS"
authors: "Alberto Garcia Jimenez and Anna Marcionetti"
date: "November 2024"
output: html_document
---

# **Introduction**

This pipeline processes raw sequencing data to identify SNPs using **ATLAS** v0.9.9 [Link et al. 2017]. It includes steps for quality control, mapping, genotype likelihood estimation, and variant calling. 
The default reference genome is the **_Amphiprion clarkii_ reference genome** (GenBank; ID: JALBFV000000000) [Moore et al., 2023].

---

# **Pipeline Overview**

Each script corresponds to a specific step in the pipeline.  
The steps encompass the entire process, starting from raw sequencing data to the generation of high-quality SNP datasets. This includes:

1. **Data Preparation**: Splitting samples by population, extracting sample names, and counting errors.
2. **Quality Control and Trimming**: Removing adapters, filtering low-quality reads, and ensuring read quality for downstream analyses.
3. **Mapping and Alignment**: Indexing the reference genome, mapping reads, and generating sorted BAM files.
4. **Genotype Likelihood Estimation**: Producing genotype likelihood files (GLF) for either the whole genome or individual chromosomes.
5. **Variant Calling**: Identifying major/minor alleles and generating VCF files.
6. **Post-Processing and Filtering**: Compressing, annotating, and filtering VCF files based on quality and allele frequency.

This modular pipeline allows for flexibility in customizing steps based on the research needs.

---

# **Installation and Setup**

## Installing ATLAS

The SNP calling is performed with **ATLAS**, which can be installed via `conda`:

```bash
module load gcc miniconda3
conda init bash  # Ignore the sudo password prompt (Ctrl+C to exit if it hangs)
conda create -n ATLAS
conda activate ATLAS
conda install -c conda-forge -c bioconda atlas
conda deactivate
```

Ensure the reference genome file is stored in an accessible location, such as:

```bash
/user/janedoe/Reference_Genomes/Amphiprion_clarkii/
```
# **Workflow Steps**

## 1. Trimming and Quality Control
Trims adapters, filters low-quality reads, and checks read quality.

Example command:

```bash
sbatch 1B_Trimming_reads.sh Samples.txt . RawData adapters/adapters.fa 20 50 10
```

## 2. Mapping and Filtering
Aligns reads to the reference genome and processes mapping outputs.

### Index the Reference Genome:

```bash
sbatch 2A_Run_BWA_RefIndexing.sh ClarkiiReference/ClarkiiGenome
```

### Map Reads and Generate BAM Files:

```bash
sbatch 2B_Mapping_generating_sorted_BAM.sh Samples.txt . TrimmedReads ClarkiiReference/ClarkiiGenome RemoveSAM
```

### Filter and Process Mapping Results:

```bash
sbatch 2C_MappingProcessing.sh Samples.txt . MappingOutput 30 /path/to/atlas
```

## 3. Genotype Likelihood Estimation

Generates Genotype Likelihood Files (GLF).

### Whole Genome:

```bash
sbatch 3A_GenotypeLikelihoods.sh Samples.txt . FilteredMapping 100000 /path/to/atlas
```
or

### By Chromosomes:

```bash
sbatch 3B_GenotypeLikelihoods_ByChromosomes.sh Samples.txt . FilteredMapping 100000 atlas Chromosomes.txt
```

## 4. Variant Calling (VCF Generation)

Identifies SNPs and generates VCF files.

### Whole Genome:

```bash
sbatch 4A_Run_MajorMinor.sh Samples.txt . GLF GLFFiles ATLAS_Clarkii_majorMinor BWA.Clarkii.Sort.Filt_mergedReads.glf.gz atlas
```

or

### By Chromosomes:

```bash
sbatch 4B_Run_MajorMinor_ByChromosomes.sh Samples.txt . GLF GLFFiles ATLAS_Clarkii_majorMinor BWA.Clarkii.Sort.Filt_mergedReads.glf.gz atlas Chromosomes.txt
```

## 5. Post-Processing

Includes VCF compression, reformatting, annotation, and filtering.

Run 5A_vcf_Statistics.sh to compute VCF statistics and 5B_produce_vcfStatistics_report.sh to produce a report.

Filter variants based on quality (5C_FilterQUAL_vcf.sh), minor allele frequency (5D_FilterMAF_vcf.sh), and other parameters.

```bash
sbatch -J QUALFilterVCF 5C_FilterQUAL_vcf.sh file.vcf.gz . . MajorMinor/ parameters.txt 
sbatch -J MAFFilterVCF 5D_FilterMAF_vcf.sh file.vcf.gz . MajorMinor/Filter MajorMinor/Filter parameters.txt 
sbatch -J VARFilterVCF 5E_FilterVAR_vcf.sh file.vcf.gz . MajorMinor/Filter MajorMinor/Filter variants.txt keep remove_ind.txt
```

# **Example Workflow Execution**

An example workflow combining all steps:

```bash
# Step 1: Quality Control
sbatch 1B_Trimming_reads.sh Samples.txt . RawData adapters/adapters.fa 20 50 10

# Step 2: Mapping
sbatch 2A_Run_BWA_RefIndexing.sh ClarkiiReference/ClarkiiGenome
sbatch 2B_Mapping_generating_sorted_BAM.sh Samples.txt . TrimmedReads ClarkiiReference/ClarkiiGenome RemoveSAM
sbatch 2C_MappingProcessing.sh Samples.txt . MappingOutput 30 /path/to/atlas

# Step 3: Genotype Likelihoods
sbatch 3A_GenotypeLikelihoods.sh Samples.txt . FilteredMapping 100000 /path/to/atlas

# Step 4: Variant Calling
sbatch 4A_Run_MajorMinor.sh Samples.txt . GLF GLFFiles ATLAS_Clarkii_majorMinor BWA.Clarkii.Sort.Filt_mergedReads.glf.gz atlas

# Step 5: Filter Variants
sbatch -J QUALFilterVCF 5C_FilterQUAL_vcf.sh Samples.vcf.gz . . MajorMinor/ parameters.txt  
sbatch -J MAFFilterVCF 5D_FilterMAF_vcf.sh Samples.vcf.gz . MajorMinor/Filter MajorMinor/Filter parameters.txt  
sbatch -J VARFilterVCF 5E_FilterVAR_vcf.sh Samples.vcf.gz . MajorMinor/Filter MajorMinor/Filter variants.txt keep remove_ind.txt
```

# **Additional Notes**

Ensure all paths (e.g., reference genome, ATLAS binaries) are correct for your system.

Edit Samples.txt and other input files as needed.

Use a high-performance computing (HPC) environment with SLURM for batch submission (sbatch commands).

# **Citation**

If you use this pipeline in your research, please cite:

Moore et al., 2023 for the reference genome, Link et al. 2017 for the use of ATLAS, and this repository and any relevant software tools.

For questions or issues, contact me at agarcia26286@gmail.com.






