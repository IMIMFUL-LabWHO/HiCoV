## Prerequisites

Before running HiCov, install:

* **Singularity**
  Installation instructions: https://singularity-tutorial.github.io/01-installation/

The following files are also required:

* `hicov.sh` — main HiCov pipeline script
* `20240710_hicov.def` — Singularity definition file used to build the container
* `reference/` directory containing:

  * `<reference>.fasta`
  * `<reference>.gff3`
  * `<reference>.primer.bed`

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

| Script                     | Purpose                                                                                                                     |
| -------------------------- | --------------------------------------------------------------------------------------------------------------------------- |
| `hicov.sh`                 | Main sequencing pipeline: trimming, mapping, primer trimming, mapping statistics, consensus generation, and variant calling |
| Coverage analysis R script | Calculates coverage metrics and generates genome-wide coverage visualizations for hCoV-229E, NL63, OC43, and HKU1           |
| `variants_data_prep.R`     | Cleans, classifies, and enriches the consolidated iVar variant table for downstream analysis                                |


## Usage

### 1. Build the Singularity container

Build the `hicov.sif` container from the supplied Singularity definition file:

```bash
sudo singularity build hicov.sif 20240710_hicov.def
```

If you do not have `sudo` privileges, build the container using `--fakeroot`:

```bash
singularity build --fakeroot hicov.sif 20240710_hicov.def
```

After a successful build, the working directory should contain:

```text
working_directory/
├── hicov.sh
├── hicov.sif
├── 20240710_hicov.def
└── reference/
    ├── <reference>.fasta
    ├── <reference>.gff3
    └── <reference>.primer.bed
```

### 2. Prepare the reference files

For a reference named, for example:

```text
hcov_229e
```

the `reference/` directory must contain:

```text
reference/hcov_229e.fasta
reference/hcov_229e.gff3
reference/hcov_229e.primer.bed
```

All reference files must use the same basename.

### 3. Prepare the input FASTQ files

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

### 4. Run HiCov

Run the pipeline inside the Singularity container:

```bash
singularity exec hicov.sif bash hicov.sh <reference_name> <fastq_directory>
```

For example:

```bash
singularity exec hicov.sif bash hicov.sh hcov_229e /home/imi.mf/asuljic/experiments/seasonal_coronas/run_id201384/data/e229
```

The first positional argument is the reference basename:

```text
hcov_229e
```

The second positional argument is the directory containing the paired-end FASTQ files.

## Main script overview

The main `hicov.sh` script automates the complete processing of paired-end Illumina sequencing data from raw FASTQ files to consensus sequences, variant calls, and sequencing statistics.

For each sample, the pipeline performs the following main steps:

1. **Sample detection**
   Identifies samples from paired-end FASTQ filenames.

2. **Read trimming and QC**
   Uses `fastp` to remove low-quality bases, adapters, poly-G/poly-X sequences, and short reads. HTML and JSON QC reports are generated for each sample.

3. **Reference preparation**
   Indexes the selected reference genome using `bwa` and `samtools`.

4. **Read mapping**
   Aligns trimmed reads to the reference genome using `bwa mem`.

5. **Primer trimming and BAM processing**
   Removes primer-derived sequence using `ivar trim`, fixes mate information, sorts alignments, removes duplicates, and produces the final indexed BAM files.

6. **Mapping and coverage statistics**
   Calculates alignment statistics, genome coverage, and per-position sequencing depth using `samtools`.

7. **Consensus generation**
   Generates a consensus sequence for each sample using `samtools mpileup` and `ivar consensus`. Low-confidence and ambiguous positions are represented as `N`.

8. **Variant calling and annotation**
   Calls variants using `ivar variants` and annotates them using the supplied reference FASTA and GFF3 files.

9. **Result consolidation**
   Combines individual sample outputs into summary files containing mapping statistics, consensus sequences, variants, and coverage data.

The main output files are written to the `results/` directory:

```text
results/
├── mapstats.tsv
├── consensus_sequences.fasta
├── sleek_variants.tsv
└── coverage.csv
```

In short, the pipeline converts:

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

