# HiCov

HiCov is a Bash-based bioinformatics pipeline for processing paired-end Illumina sequencing data from human coronavirus (**hCoV**) samples.

The pipeline performs read trimming and quality control, reference mapping, primer trimming, mapping statistics, consensus sequence generation, and variant calling. Additional R scripts are included for downstream coverage analysis, visualization, and variant-data preparation.

HiCov was developed through collaboration between **IMI MF, University of Ljubljana** and the **Hiscox Lab, University of Liverpool**.

**Author:** Alen Suljič (alen.suljic@mf.uni-lj.si)

---

## Prerequisites

Before running HiCov, install:

* **Singularity**
  Installation instructions: https://singularity-tutorial.github.io/01-installation/

The repository provides:

* `hicov.sh` — main HiCov analysis pipeline
* `20240710_hicov.def` — Singularity definition file used to build the analysis container
* `reference/` — reference genomes, genome annotations, and primer coordinates required by the pipeline
* R scripts for downstream coverage and variant analysis

The `reference/` directory contains corresponding:

```text
<reference>.fasta
<reference>.gff3
<reference>.primer.bed
```

files used for reference mapping, variant annotation, and primer trimming.

The required reference files are already included in the repository and do not need to be downloaded or prepared separately.

The Singularity container provides the bioinformatics software used by the main pipeline, including:

* `fastp`
* `bwa`
* `samtools`
* `ivar`
* `seqtk`

The downstream R scripts require:

* R
* `tidyverse`

---

## Repository overview

HiCov consists of a primary Bash pipeline for sequencing-data processing and two R scripts for downstream analysis and visualization.

```text
Raw paired-end FASTQ files
          │
          ▼
       hicov.sh
          │
          ├────────► QC reports
          │
          ├────────► Final BAM files
          │
          ├────────► Consensus sequences
          │
          ├────────► Mapping statistics
          │
          ├────────► coverage.csv
          │                │
          │                ▼
          │        Coverage analysis R script
          │                │
          │                ▼
          │        Genome coverage figures
          │
          └────────► sleek_variants.tsv
                           │
                           ▼
                 variants_data_prep.R
                           │
                           ▼
                 variants_enhanced.csv
```

### Scripts

| Script                     | Purpose                                                                                                                          |
| -------------------------- | -------------------------------------------------------------------------------------------------------------------------------- |
| `hicov.sh`                 | Main sequencing pipeline: read trimming, mapping, primer trimming, mapping statistics, consensus generation, and variant calling |
| `coverage_visualisation.R`   | Calculates coverage metrics and generates genome-wide coverage visualizations for hCoV-229E, NL63, OC43, and HKU1                |
| `variants_data_prep.R`     | Cleans, classifies, and enriches the consolidated iVar variant table for downstream analysis                                     |

---

# Usage

## 1. Build the Singularity container

Build the `hicov.sif` container from the supplied Singularity definition file:

```bash
sudo singularity build hicov.sif 20240710_hicov.def
```

If you do not have `sudo` privileges, build the container using `--fakeroot`:

```bash
singularity build --fakeroot hicov.sif 20240710_hicov.def
```

After a successful build, the repository should contain approximately:

```text
HiCov/
├── hicov.sh
├── hicov.sif
├── 20240710_hicov.def
├── reference/
│   ├── <reference>.fasta
│   ├── <reference>.gff3
│   └── <reference>.primer.bed
├── variants_data_prep.R
└── ...
```

---

## 2. Prepare the input FASTQ files

Input data must consist of paired-end, gzip-compressed FASTQ files following the naming convention:

```text
<sample>_R1.fastq.gz
<sample>_R2.fastq.gz
```

For example:

```text
sample01_R1.fastq.gz
sample01_R2.fastq.gz
sample02_R1.fastq.gz
sample02_R2.fastq.gz
```

Sample names are inferred from the part of each filename preceding the first underscore.

For example:

```text
sample01_R1.fastq.gz
```

is interpreted as:

```text
sample01
```

> **Important:** Sample identifiers should therefore not contain underscores unless the sample-detection logic in `hicov.sh` is modified.

---

## 3. Run HiCov

Run the pipeline inside the Singularity container:

```bash
singularity exec hicov.sif bash hicov.sh <reference_name> <fastq_directory>
```

For example:

```bash
singularity exec hicov.sif bash hicov.sh hcov_229e /home/imi.mf/asuljic/experiments/seasonal_coronas/run_id201384/data/e229
```

The first positional argument specifies the basename of one of the reference sets included in the repository.

For example:

```text
hcov_229e
```

corresponds to:

```text
reference/hcov_229e.fasta
reference/hcov_229e.gff3
reference/hcov_229e.primer.bed
```

The second positional argument specifies the directory containing the paired-end FASTQ files.

---

# Main script overview

The main `hicov.sh` script automates the complete processing of paired-end Illumina sequencing data from raw FASTQ files to consensus sequences, variant calls, and sequencing statistics.

For each sample, the pipeline performs the following main steps:

1. **Sample detection**
   Identifies samples from paired-end FASTQ filenames.

2. **Read trimming and QC**
   Uses `fastp` to remove low-quality bases, adapters, poly-G/poly-X sequences, and short reads. HTML and JSON QC reports are generated for each sample.

3. **Reference indexing**
   Indexes the selected reference genome using `bwa` and `samtools`.

4. **Read mapping**
   Aligns trimmed paired-end reads to the selected reference genome using `bwa mem`.

5. **Primer trimming and BAM processing**
   Removes primer-derived sequences using `ivar trim`, fixes mate information, sorts the alignments, removes duplicates, and produces final indexed BAM files.

6. **Mapping and coverage statistics**
   Calculates alignment statistics, genome coverage, and per-position sequencing depth using `samtools`.

7. **Consensus generation**
   Generates a consensus sequence for each sample using `samtools mpileup` and `ivar consensus`. Ambiguous nucleotide calls and gap characters are converted to `N`.

8. **Variant calling and annotation**
   Calls variants using `ivar variants` and annotates them using the selected reference FASTA and GFF3 files.

9. **Result consolidation**
   Combines individual sample outputs into summary files containing mapping statistics, consensus sequences, variants, and coverage data.

The workflow can be summarized as:

```text
Paired-end FASTQ
      │
      ▼
Read trimming + QC
      │
      ▼
Reference mapping
      │
      ▼
Primer trimming + BAM processing
      │
      ├──────────────► Mapping / coverage statistics
      │
      ├──────────────► Consensus sequences
      │
      └──────────────► Variant calling
                              │
                              ▼
                       Consolidated results
```

---

## Pipeline parameters

The main configurable parameters are defined near the beginning of `hicov.sh`:

```bash
# thread count
thr=24

# minimum Phred quality score
qqp=20

# minimum read length after trimming
lr=30

# minimum read depth for consensus
dc=10
```

| Parameter | Default | Description                                       |
| --------- | ------: | ------------------------------------------------- |
| `thr`     |    `24` | Number of processing threads                      |
| `qqp`     |    `20` | Minimum Phred quality threshold                   |
| `lr`      |    `30` | Minimum read length retained after trimming       |
| `dc`      |    `10` | Minimum read depth required for consensus calling |

These values can be modified directly in `hicov.sh` before execution.

---

# Output files

The main consolidated output files generated by `hicov.sh` are written to the `results/` directory:

```text
results/
├── mapstats.tsv
├── consensus_sequences.fasta
├── sleek_variants.tsv
└── coverage.csv
```

## `mapstats.tsv`

Contains consolidated mapping and sequencing statistics for all samples.

The table includes:

```text
sample
rname
startpos
endpos
numreads
covbases
coverage
meandepth
meanbaseq
meanmapq
r1_nreads
r2_nreads
```

The `r1_nreads` and `r2_nreads` columns contain the number of reads remaining in the trimmed R1 and R2 FASTQ files.

---

## `consensus_sequences.fasta`

Multi-FASTA file containing the final consensus sequence generated for each sample.

Consensus calling uses an allele-frequency threshold of:

```text
0.5
```

and, by default, a minimum sequencing depth of:

```text
10
```

Ambiguous IUPAC nucleotide codes and gap characters are converted to `N`.

---

## `sleek_variants.tsv`

Combined iVar variant table containing variants detected across all samples.

Individual sample variant files generated by `ivar variants` are consolidated into this table for downstream analysis.

---

## `coverage.csv`

Combined per-position sequencing-depth data generated from `samtools depth`.

This file is used by the downstream coverage-analysis R script.

---

# Complete output structure

After a successful analysis, the working directory will contain approximately:

```text
working_directory/
├── hicov.sh
├── hicov.sif
├── 20240710_hicov.def
│
├── reference/
│   ├── <reference>.fasta
│   ├── <reference>.gff3
│   ├── <reference>.primer.bed
│   └── reference index files
│
├── trimmed/
│   ├── <sample>_trim_R1.fastq.gz
│   └── <sample>_trim_R2.fastq.gz
│
├── qc/
│   ├── <sample>.fastp.html
│   └── <sample>.fastp.json
│
├── mappings/
│   ├── <sample>_final.bam
│   └── <sample>_final.bam.bai
│
├── stats/
│   ├── <sample>_stats.log
│   └── <sample>.covdepth
│
├── consensus/
│   └── <sample>.fasta
│
├── variants/
│   └── <sample>.tsv
│
├── results/
│   ├── mapstats.tsv
│   ├── consensus_sequences.fasta
│   ├── sleek_variants.tsv
│   └── coverage.csv
│
└── logs/
    ├── experiment.log
    └── samples
```

---

# Downstream analysis

The repository contains two R scripts for downstream processing of the files generated by `hicov.sh`.

---

## Coverage analysis

The coverage analysis R script processes `coverage.csv` files generated by HiCov and summarizes sequencing coverage across the viral reference genomes.

The script currently includes virus-specific analyses for:

* **hCoV-229E**
* **hCoV-NL63**
* **hCoV-OC43**
* **hCoV-HKU1**

The script requires:

```r
library(tidyverse)
```

For each coronavirus, the script:

1. loads the corresponding `coverage.csv` file
2. renames the input columns to sample, genomic position, and sequencing depth
3. calculates mean and median sequencing depth for individual genomic positions
4. calculates overall sequencing-depth statistics
5. determines whether genomic positions have sufficient coverage for consensus generation
6. calculates consensus completeness for each sample
7. assigns genomic positions to virus-specific annotated genes and genomic regions
8. calculates gene-level coverage statistics
9. generates genome-wide coverage plots
10. exports publication-quality TIFF figures

The analysis can be summarized as:

```text
results/coverage.csv
        │
        ▼
Load coverage data
        │
        ▼
Calculate per-position coverage
        │
        ├────► Mean / median sequencing depth
        │
        ├────► Consensus completeness
        │
        └────► Per-sample coverage statistics
        │
        ▼
Assign genomic regions
        │
        ▼
Calculate region-level coverage
        │
        ▼
Generate genome-wide coverage plot
        │
        ▼
coverage_<virus>.tiff
```

The script uses virus-specific genome coordinates to annotate regions including:

```text
5'UTR
ORF1ab
S
E
M
N
3'UTR
```

along with additional virus-specific genes and ORFs.

For example, the genome annotations differ between hCoV-229E, NL63, OC43, and HKU1.

A genomic position is considered sufficiently covered for consensus-completeness calculations when its sequencing depth is greater than 9 reads, corresponding to a minimum depth of 10.

### Input

The script reads the `coverage.csv` output generated by HiCov.

An example directory organization is:

```text
coverage_data/
├── e229/
│   └── coverage.csv
├── nl63/
│   └── coverage.csv
├── oc43/
│   └── coverage.csv
└── HKU1/
    └── coverage.csv
```

Before running the script, replace:

```r
setwd("/path/to/coverage/data")
```

with the appropriate directory containing the HiCov coverage results.

### Output

The script generates:

```text
coverage_229e.tiff
coverage_nl63.tiff
coverage_oc43.tiff
coverage_hku1.tiff
```

The output location is controlled by the `path` argument supplied to `ggsave()` in the script and should be adjusted for the local analysis environment.

---

## Variant data preparation

The `variants_data_prep.R` script performs downstream processing of the consolidated iVar variant table generated by `hicov.sh`.

The script requires:

```r
library(tidyverse)
```

and reads:

```text
sleek_variants.tsv
```

as input.

The script performs the following main steps:

### 1. Load variants

The consolidated iVar variant table is loaded using:

```r
d <- read_tsv("sleek_variants.tsv")
```

### 2. Filter variants

By default, only variants satisfying:

```r
PASS == TRUE
```

are retained.

The filtering line can be removed from the R script if variants should be analysed regardless of the Fisher's exact test result.

### 3. Classify variants by allele frequency

Variants are divided into two frequency classes:

```text
ALT_FREQ >= 0.5  → major
ALT_FREQ < 0.5   → minor
```

### 4. Extract gene annotations

Gene information is extracted from the iVar `GFF_FEATURE` annotation.

### 5. Generate nucleotide mutation identifiers

Reference nucleotide, genomic position, and alternate allele information are combined to generate a compact nucleotide mutation identifier.

### 6. Classify mutation type

Variants are categorized as:

* `Deletion`
* `Insertion`
* `SNP`
* `Non-coding`

### 7. Classify coding substitutions

Coding SNPs are further classified as:

* `Synonymous`
* `Non-synonymous`

### 8. Remove duplicate variants

Duplicate sample/mutation combinations are removed from the dataset.

### 9. Export the enhanced variant table

The resulting dataset is exported as:

```text
variants_enhanced.csv
```

The workflow can be summarized as:

```text
results/sleek_variants.tsv
           │
           ▼
      Load iVar results
           │
           ▼
    Filter PASS variants
           │
           ▼
  Classify allele frequency
      │             │
    major          minor
           │
           ▼
    Extract gene annotation
           │
           ▼
    Identify mutation type
     │       │       │
    SNP   Insertion Deletion
     │
     ├── Synonymous
     └── Non-synonymous
           │
           ▼
     Remove duplicates
           │
           ▼
 variants_enhanced.csv
```

The script contains placeholder input and output paths that should be adjusted before running the analysis.

---

# Logging

The main pipeline writes standard output and standard error to:

```text
experiment.log
```

using:

```bash
exec > >(tee -a "experiment.log") 2>&1
```

At the end of the analysis, the log file and generated sample list are moved to:

```text
logs/
├── experiment.log
└── samples
```

The log file can be used to inspect individual processing steps and troubleshoot failed or incomplete analyses.

---

# Important assumptions

The current implementation assumes that:

1. sequencing data consist of paired-end Illumina reads
2. FASTQ files follow the `<sample>_R1.fastq.gz` / `<sample>_R2.fastq.gz` naming convention
3. sample identifiers do not contain underscores
4. an appropriate reference set included in the repository is selected when running `hicov.sh`
5. corresponding FASTA, GFF3, and primer BED files use the same reference basename
6. the script is executed from the working/output directory
7. the `reference/` directory is located directly inside the working directory
8. the required bioinformatics software is available through the HiCov Singularity container
9. sufficient CPU, memory, and disk space are available for BAM sorting and intermediate files
10. the `reference/` directory is writable because reference index files are generated during the analysis

When repeating an analysis, using a new working directory is recommended to avoid mixing files generated by different runs or parameter configurations.

---

# Citation

Ošep, A., Goldswain, H., Suljič, A. et al. Development of type-specific amplicon schemes for whole-genome sequencing of seasonal human coronaviruses from clinical samples. Sci Rep (2026). https://doi.org/10.1038/s41598-026-63549-1

---
