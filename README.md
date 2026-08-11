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

