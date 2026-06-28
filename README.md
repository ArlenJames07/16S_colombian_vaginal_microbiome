# 16S colombian vaginal microbiome

**QIIME 2 workflows for 16S rRNA amplicon analysis of the Colombian vaginal microbiome, with emphasis on Gardnerella-dominant community states, Fannyhessea-associated transitions, and bacterial vaginosis-related microbial profiles.**

---

## Repository information

| Field               | Description                                                         |
| ------------------- | ------------------------------------------------------------------- |
| **Pipeline author** | Arlen James Mosquera-Ruiz                                           |
| **Institution**     | Pontificia Universidad Javeriana Cali                               |
| **ORCID**           | [0009-0008-0796-9099](https://orcid.org/0009-0008-0796-9099)        |
| **Contact**         | [arlen22@javerianacali.edu.co](mailto:arlen22@javerianacali.edu.co) |

---

## Overview

This repository contains the analysis scripts supporting the manuscript:

**“A clinically heterogeneous Gardnerella-dominant vaginal microbiome state in Colombian women reveals a Fannyhessea-associated transition to bacterial vaginosis.”**

The workflow was developed for **16S rRNA amplicon sequencing data** from Colombian vaginal microbiome samples, integrating:

* paired-end FASTQ organization
* metadata validation
* QIIME 2 data import
* demultiplexed read quality inspection
* DADA2 denoising
* feature-table generation
* representative-sequence extraction
* denoising-statistics visualization
* taxonomic profiling
* vaginal community-state analysis
* Gardnerella-dominant microbiome-state characterization
* Fannyhessea-associated bacterial vaginosis transition analysis
* downstream diversity and clinical-metadata integration

The repository is intended to support transparent reuse of the custom scripts used to analyse vaginal microbiome structure in Colombian women. Raw sequencing files and generated QIIME 2 artifacts are not stored in this GitHub repository.

---

## Scientific scope

This pipeline investigates clinically heterogeneous vaginal microbiome profiles in Colombian women, with particular focus on Gardnerella-dominant states and their relationship with bacterial vaginosis-associated taxa.

The main analyses include:

1. import and quality control of paired-end 16S rRNA amplicon data
2. denoising and amplicon sequence variant reconstruction using QIIME 2/DADA2
3. taxonomic assignment of vaginal bacterial communities
4. identification of Gardnerella-dominant microbiome states
5. evaluation of Fannyhessea-associated microbial transitions toward bacterial vaginosis
6. integration of microbiome profiles with clinical and metadata variables
7. generation of reproducible outputs for downstream statistical analysis and manuscript figures

---

## Repository structure

Current repository layout:

```text
16S_colombian_vaginal_microbiome/
├── README.md
├── .gitignore
├── data/
│   └── input_files/
│       ├── metadata.tsv
│       └── manifest.tsv
├── envs/
│   └── setup_qiime2_rachis_2026.4.sh
├── results/
│   └── .gitkeep
└── scripts/
    └── qiime_pipeline.py
```

Raw FASTQ files are expected locally under:

```text
data/input_files/
```

but are intentionally excluded from GitHub tracking.

---

## Pipeline modules

| Module             | Script / file                   | Purpose                                                                                |
| ------------------ | ------------------------------- | -------------------------------------------------------------------------------------- |
| `data/input_files` | `metadata.tsv`                  | Sample metadata used for QIIME 2 summaries and downstream clinical/microbiome analysis |
| `data/input_files` | `manifest.tsv`                  | QIIME 2 manifest file linking sample IDs with forward and reverse FASTQ paths          |
| `envs`             | `setup_qiime2_rachis_2026.4.sh` | Installation instructions for the QIIME 2 Conda environment                            |
| `scripts`          | `qiime_pipeline.py`             | Main Python wrapper for executing QIIME 2 commands                                     |
| `results`          | `.gitkeep`                      | Placeholder directory for local QIIME 2 outputs; generated files are not tracked       |

---

## Software requirements

This workflow requires:

* Linux or Windows WSL2
* Git
* Miniconda or Anaconda
* QIIME 2
* Python 3
* sufficient local storage for FASTQ, `.qza`, `.qzv`, and intermediate result files

The recommended QIIME 2 environment is:

```text
rachis-qiime2-2026.4
```

---

## Installing QIIME 2

Install Miniconda and initialize Conda:

```bash
conda init
```

Close and reopen the terminal.

Update Conda:

```bash
conda update conda
```

Create the QIIME 2 environment:

```bash
conda env create \
  --name rachis-qiime2-2026.4 \
  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.4/qiime2/released/rachis-qiime2-linux-64-conda.yml
```

Activate the environment:

```bash
conda activate rachis-qiime2-2026.4
```

Verify the installation:

```bash
qiime info
```

---

## Input data

The expected input data are demultiplexed paired-end FASTQ files.

Expected naming format:

```text
SampleID_R1.fastq.gz
SampleID_R2.fastq.gz
```

Example:

```text
CGV001_R1.fastq.gz
CGV001_R2.fastq.gz
```

The manifest file should be located at:

```text
data/input_files/manifest.tsv
```

Expected QIIME 2 paired-end manifest format:

```text
sample-id	forward-absolute-filepath	reverse-absolute-filepath
CGV001	/home/rare/arlen/16S_colombian_vaginal_microbiome/data/input_files/CGV001_R1.fastq.gz	/home/rare/arlen/16S_colombian_vaginal_microbiome/data/input_files/CGV001_R2.fastq.gz
CGV003	/home/rare/arlen/16S_colombian_vaginal_microbiome/data/input_files/CGV003_R1.fastq.gz	/home/rare/arlen/16S_colombian_vaginal_microbiome/data/input_files/CGV003_R2.fastq.gz
```

The metadata file should be located at:

```text
data/input_files/metadata.tsv
```

---

## Importing FASTQ files into QIIME 2

QIIME 2 requires sequencing data to be imported as `.qza` artifacts before downstream analysis.

Import paired-end reads:

```bash
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path data/input_files/manifest.tsv \
  --output-path results/demux-paired-end.qza \
  --input-format PairedEndFastqManifestPhred33V2
```

Summarize demultiplexed reads:

```bash
qiime demux summarize \
  --i-data results/demux-paired-end.qza \
  --o-visualization results/demux-paired-end.qzv
```

The `.qzv` file can be inspected using QIIME 2 View.

---

## Denoising with DADA2

After inspecting the demultiplexed quality plots, denoise the paired-end reads with DADA2.

Initial truncation values used as a starting point:

```text
--p-trunc-len-f 240
--p-trunc-len-r 200
```

Run DADA2:

```bash
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs results/demux-paired-end.qza \
  --p-trunc-len-f 240 \
  --p-trunc-len-r 200 \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-n-threads 0 \
  --o-table results/table.qza \
  --o-representative-sequences results/rep-seqs.qza \
  --o-denoising-stats results/denoising-stats.qza
```

Summarize denoising statistics:

```bash
qiime metadata tabulate \
  --m-input-file results/denoising-stats.qza \
  --o-visualization results/denoising-stats.qzv
```

Summarize the feature table:

```bash
qiime feature-table summarize \
  --i-table results/table.qza \
  --o-visualization results/table.qzv \
  --m-sample-metadata-file data/input_files/metadata.tsv
```

Summarize representative sequences:

```bash
qiime feature-table tabulate-seqs \
  --i-data results/rep-seqs.qza \
  --o-visualization results/rep-seqs.qzv
```

---

## Running the Python pipeline

The main pipeline script is located at:

```text
scripts/qiime_pipeline.py
```

Activate the QIIME 2 environment:

```bash
conda activate rachis-qiime2-2026.4
```

Run the pipeline from the repository root:

```bash
python3 scripts/qiime_pipeline.py
```

---

## Expected outputs

Main local QIIME 2 outputs:

```text
results/
├── demux-paired-end.qza
├── demux-paired-end.qzv
├── table.qza
├── table.qzv
├── rep-seqs.qza
├── rep-seqs.qzv
├── denoising-stats.qza
└── denoising-stats.qzv
```

These files are generated locally and should not be committed to GitHub.

---

## Data availability and privacy

Raw sequencing data are **not included** in this GitHub repository because of file size, data-governance considerations, and the need to preserve participant confidentiality.

This repository may include:

* analysis scripts
* environment setup files
* manifest templates
* metadata templates or non-identifiable metadata
* documentation
* non-identifiable processed summaries, when applicable

This repository does **not** include:

* raw FASTQ files
* QIIME 2 `.qza` artifacts
* QIIME 2 `.qzv` visualizations
* large intermediate results
* participant-identifiable metadata
* private clinical information

Raw and processed data should be stored on the local server, institutional storage system, or an approved public/private repository according to the project’s ethics and data-sharing approvals.

---


## Reproducibility

This repository is intended to store the code, environment instructions, metadata structure, and workflow documentation needed to reproduce the analysis.

Raw FASTQ files and generated QIIME 2 artifacts should remain outside GitHub and be managed through controlled local or institutional storage.

---

## Project status

Current workflow components:

* paired-end FASTQ organization
* QIIME 2 environment setup
* FASTQ import into QIIME 2
* demultiplexed read quality summary
* DADA2 denoising
* feature-table summary
* representative-sequence summary
* denoising-statistics visualization

Planned or downstream components:

* taxonomic assignment
* Gardnerella-dominant state characterization
* Fannyhessea-associated bacterial vaginosis transition analysis
* alpha diversity analysis
* beta diversity analysis
* community state type analysis
* association analysis with clinical and nutritional metadata
* figure generation and manuscript reporting

---

## Citation

If you use this pipeline, code, or processed summary outputs, please cite:

Mosquera-Ruiz A. (2026). *A clinically heterogeneous Gardnerella-dominant vaginal microbiome state in Colombian women reveals a Fannyhessea-associated transition to bacterial vaginosis*. Manuscript in preparation / under review.

---

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.

---

## Contact

**Arlen James Mosquera-Ruiz**
Doctorate in Engineering and Applied Sciences
Pontificia Universidad Javeriana Cali
Cali, Colombia

Email: [arlen22@javerianacali.edu.co](mailto:arlen22@javerianacali.edu.co)
ORCID: [0009-0008-0796-9099](https://orcid.org/0009-0008-0796-9099)
