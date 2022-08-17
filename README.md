# Meta-scRNA-seq
Author: Michael Wang (fw262@cornell.edu)

![temp](https://user-images.githubusercontent.com/56937181/185178404-809f2327-d4a5-439d-8b90-5db2511f6f64.png)

## Outline
We developed meta-scRNA-seq, a pipeline for unbiased detection of non-host transcriptomic information from scRNA-seq data. To achieve this, meta-scRNA-seq aligns scRNA-seq data against the host-genome reference using standard approaches, collected single-cell tagged unmapped reads, labeled them based on sequence similarity against a large metagenomic database, and demultiplexed the reads to generate a cell-by-metagenome count matrix in parallel with the standard cell-by-gene (host) matrix.

## Required Software
This workflow requires the following packages listed below. Please ensure that tool can be called from the command line (i.e. the paths to each tool is in your path variable).

### 1. [Snakemake](https://snakemake.readthedocs.io/en/stable/)
### 2. [STAR Aligner](https://github.com/alexdobin/STAR/releases)
```
conda install -c bioconda star
```
### 3. [R, version 3.6 or greater](https://www.r-project.org/)
Please also ensure that you have downloaded the following R packages. They will be used throughout the pipeline.
- [Seurat, version 3](https://satijalab.org/seurat/install.html)
- [dplyr](https://www.r-project.org/nosvn/pandoc/dplyr.html)
- [plyr](https://www.rdocumentation.org/packages/plyr/versions/1.8.7)
- [ggplot2](https://ggplot2.tidyverse.org/)
- [argparse](https://cran.r-project.org/web/packages/argparse/index.html)
- [DropletUtils](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html)
### 4. [Samtools](http://www.htslib.org/)
```
conda install -c bioconda samtools
```
### 5. [Kraken2](https://ccb.jhu.edu/software/kraken2/)
Please make sure this tool is available in your working environment. Please also download the reference database.

## Procedure

### 1. Clone this repository.
Run the following command in your command line.
```
git clone https://github.com/fw262/Meta-scRNA-seq.git
```

### 2. Download required software listed above.

Please ensure to include all required software before starting.

### 3. Store or link paired end sequencing files.

Please move raw fastq files for each experiment into one data directory. Please ensure the sequence files end in "{sample}\_R1_001.fastq.gz" and "{sample}\_R1_001.fastq.gz" in your data directory.

### 4. Create the STAR reference of the host genome.

### 5. Edit the config.yaml file for your experiment.

Please change the variable names in the config.yaml as required for your analysis. This includes the following changes:
- **Samples**: Samples prefix (before the \_R1_001.fastq.gz)
- **STAR_IND**: Path to your STAR generated index folder.
- **DATADIR**: Path to where the sequencing samples ({sample}\_R1_001.fastq.gz) are stored.
- **PIPELINE_MAJOR**: Directory where the outputs (expression matrices, plots) are stored.
- **GLOBAL**: Define global variables for pipeline including number of mismatches allowed in STAR, cell barcode base pair range in read 1, and UMI base pair range in read 1.
- **STAREXEC**: Path to STAR.
- **KRAKEN**: Path to Kraken2.
- **KRAKEN_DB**: Path to Kraken2 database.

- **CORES**: Number of cores used in each step of the pipeline. To run multiple samples in parallel, please specify total number of cores in the snakemake command (i.e. "snakemake -j {total cores}").

### 6. Run snakemake with the command "snakemake".

Please ensure the Snakefile and config.yaml files as well as the scripts folder are in the directory where you intend to run the pipeline.

## Output
- Merged transcriptome + metagenomice expression matrices are stored in "**[PIPELINE_MAJOR]/[Samples]_solo/Solo.out/merged**" folder.
- The "**[PIPELINE_MAJOR]/[Samples]_solo/plots**" folder contains several useful plots including UMAP projection of the data, level of unmapped reads for each cell cluster, as well as cell-cluster specific expression of all metagenomic features, differentially expressed genes, and differentially expressed metagenomic features.
