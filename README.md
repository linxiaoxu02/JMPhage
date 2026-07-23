<p align="center">
  <img src="img/JMPhage_logo.png" alt="JMPhage logo" height="300">
</p>

<h1 align="center">JMPhage</h1>

<p align="center">
  <b>A modular one-stop pipeline for annotation, characterization, and comparative analysis of tailed phages (<i>Caudoviricetes</i>)</b>
</p>

<p align="center">
  <a href="https://doi.org/10.5281/zenodo.20838733">
    <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.20838733.svg" alt="Zenodo DOI">
  </a>
  <a href="https://www.gnu.org/licenses/gpl-3.0">
    <img src="https://img.shields.io/badge/License-GPL%20v3-blue.svg" alt="License: GPL v3">
  </a>
  <img src="https://img.shields.io/badge/version-v1.0-blue" alt="Version">
  <img src="https://img.shields.io/badge/platform-Linux-lightgrey" alt="Platform">
</p>

---

## Overview

JMPhage is a command-line workflow for one-stop analysis of tailed bacteriophage (Caudoviricetes) genomes. Starting from complete or near-complete assemblies derived from cultured isolates or metagenomic data, with optional paired-end sequencing reads, JMPhage connects three analytical modules in a single reproducible workflow: a mapping module for read alignment, sequencing-depth assessment, and genome-terminus inference; an annotation module for ORF prediction, functional annotation, and circular genome visualization; and a characterization module for protein-sharing network construction, automated reference selection, gene synteny comparison, intergenomic similarity analysis, and phylogenetic placement. The pipeline is designed for single-genome analyses and small multi-genome batches.

<p align="center">
  <img src="img/Figure 1.png" alt="JMPhage workflow" width="850">
</p>

## Main Features

- **One-stop tailed-phage workflow**: run mapping, annotation, and characterization with a single command.
- **Flexible input modes**: analyze one genome or a directory containing up to five phage genomes per run.
- **Optional sequencing-read support**: map paired-end reads, generate depth profiles, and infer packaging mechanisms with PhageTerm.
- **Functional annotation and genome maps**: predict ORFs with Pharokka and annotate protein functions using VOG-based DIAMOND/HMMER searches.
- **Comparative characterization**: build protein-sharing networks with vConTACT2 and visualize related genomes with Clinker.
- **Species-level and phylogenetic context**: compute intergenomic similarity using VIRIDIC and construct GLUVAB-based phylogenetic trees.
- **Joint multi-phage characterization**: for batch input, JMPhage builds a shared network once, then performs per-phage and joint downstream analyses.
- **Restart-aware execution**: completed steps are detected and skipped on re-run.

## Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Input Requirements](#input-requirements)
- [Usage](#usage)
- [Output Structure](#output-structure)
- [JMPhage Database](#jmphage-database)
- [Integrated Software and References](#integrated-software-and-references)
- [Re-running and Resuming](#re-running-and-resuming)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [Data and Code Availability](#data-and-code-availability)
- [License](#license)
- [Contact](#contact)

## Installation

JMPhage requires Linux and Conda. The repository provides six Conda environment files under `jm_envs/` and a built-in `install` subcommand for environment deployment and database installation.

### System Requirements

| Requirement | Recommendation |
|:--|:--|
| Operating system | Linux x86_64; tested on Ubuntu 18.04 LTS and CentOS 7 |
| Python for launcher | Python 3 |
| Conda | Miniconda3, Anaconda3, or compatible Conda distribution |
| Disk space | ~16 GB for Conda environments; ~16 GB for `jm_db` after extraction |
| Database download | ~4.9 GB compressed archive |
| CPU | Default `-t 4`; increase according to available cores |
| Memory | Depends on batch size and database search steps; vConTACT2 and HMMER are typically the most memory-intensive |

### Step 1: Clone and Install the Launcher

The `jmphage` launcher is a lightweight Python command and does not require a dedicated Conda environment. You can install and run it directly from your `base` environment, or from any Python 3 environment you prefer. To keep things isolated, you may optionally create a dedicated environment first:

```bash
conda create -n JMPhage python=3
conda activate JMPhage
```

Then clone the repository and install the launcher:

```bash
git clone https://github.com/linxiaoxu02/JMPhage.git
cd JMPhage
pip install .
```

After installation, the command-line entry point is:

```bash
jmphage -h
jmphage --help
jmphage --version
```

Expected version output:

```text
JMPhage v1.0
```

### Step 2: Install Conda Environments

JMPhage uses six isolated environments:

| Environment | Purpose | Main tools |
|:--|:--|:--|
| `JM_mapping` | Read mapping and depth profiling | Bowtie2, SAMtools |
| `JM_annotation` | ORF prediction and functional annotation | Pharokka, DIAMOND, HMMER |
| `JM_characterization` | Protein network and comparative genomics | vConTACT2, Clinker, Pharokka |
| `jmp_Phageterm` | Packaging mechanism and terminus detection | PhageTerm |
| `jmp_VIRIDIC` | Intergenomic similarity analysis | VIRIDIC, R |
| `jmp_GLUVAB` | Phylogenetic tree construction and visualization | GLUVAB, ggtree |

For users in mainland China, we strongly recommend configuring the Tsinghua Conda mirrors before installing the environments. In our tests across multiple servers, default Conda channels were frequently slow or unstable and often caused environment installation to hang or fail. To add the mirrors:

```bash
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/pytorch/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
conda config --add channels bioconda
```

You can verify the resulting channel order with `conda config --show channels`.

Install all environments:

```bash
jmphage install -e ./jm_envs
```

If an environment already exists, JMPhage skips it. To reinstall an environment, remove it first with `conda env remove -n <environment_name>`.

Depending on your network conditions and local Conda configuration, one or more environments may occasionally fail to install. In that case, we recommend using [mamba](https://github.com/mamba-org/mamba) to install the failed environment manually from its YAML file under `jm_envs/`. Mamba resolves dependencies much faster and more reliably than Conda.

If mamba is not yet available, install it into your `JMPhage` launcher environment (or any Conda environment you prefer):

```bash
conda install -n JMPhage -c conda-forge mamba
```

Then create the failed environment from its YAML file:

```bash
mamba env create -n <environment_name> -f jm_envs/<environment_name>.yml
```

The environment name must match exactly one of the six names listed in the table above (`JM_mapping`, `JM_annotation`, `JM_characterization`, `jmp_Phageterm`, `jmp_VIRIDIC`, `jmp_GLUVAB`). JMPhage locates each Conda environment by its hardcoded name at runtime, so renaming an environment will cause the corresponding pipeline step to fail.

For example, to reinstall `JM_characterization`:

```bash
mamba env create -n JM_characterization -f jm_envs/JM_characterization.yml
```

After the failed environments are installed manually, you can re-run `jmphage install -e ./jm_envs` to confirm that all six environments are ready. JMPhage will skip the ones that are already present.

### Step 3: Install the Reference Database

The JMPhage database, `jm_db`, is hosted on Zenodo.

```bash
jmphage install -d /path/to/jm_db
```

This command downloads `jm_db.tar.gz`, verifies its MD5 checksum, extracts the archive, and writes an installation marker file.

- Zenodo DOI: [10.5281/zenodo.20838733](https://doi.org/10.5281/zenodo.20838733)
- Compressed size: ~4.9 GB
- Extracted size: ~16 GB
- Expected MD5: `20f8db131a251b150c5598ac192f6c15`

For offline servers, download `jm_db.tar.gz` manually from Zenodo. You can grab it through your browser from the [Zenodo record page](https://doi.org/10.5281/zenodo.20838733), or use a multi-connection accelerator such as [axel](https://github.com/axel-download-accelerator/axel) to speed up the download on slow networks. If axel is not already available, install it into your `JMPhage` launcher environment (or any Conda environment you prefer):

```bash
conda install -n JMPhage -c conda-forge axel
axel -n 12 https://zenodo.org/records/20838733/files/jm_db.tar.gz
```

Then verify the MD5 checksum and extract the archive:

```bash
md5sum jm_db.tar.gz     # Expected: 20f8db131a251b150c5598ac192f6c15
tar -zxvf jm_db.tar.gz
```

### One-step Installation

```bash
jmphage install -e ./jm_envs -d /path/to/jm_db
```

## Quick Start

### Full Pipeline with Paired-end Reads

```bash
jmphage all \
  -i my_phage.fasta \
  -1 my_phage.R1.fastq \
  -2 my_phage.R2.fastq \
  -d /path/to/jm_db \
  -o results \
  -t 8
```

By default, PhageTerm reports are generated but the downstream annotation and characterization steps continue to use the original input FASTA. To use PhageTerm-reoriented sequences when they are detected:

```bash
jmphage all \
  -i my_phage.fasta \
  -1 my_phage.R1.fastq \
  -2 my_phage.R2.fastq \
  -d /path/to/jm_db \
  -o results \
  -t 8 \
  --use-phageterm
```

### Genome-only Analysis

For genomes without sequencing reads:

```bash
jmphage no_mapping \
  -i my_phage.fasta \
  -d /path/to/jm_db \
  -o results \
  -t 8
```

### Batch Analysis

Pass a directory containing multiple `.fasta` files:

```bash
jmphage all \
  -i input_dir/ \
  -d /path/to/jm_db \
  -o results \
  -t 8
```

If reads are available, JMPhage searches for paired files with matching names in the input directory or in the directory supplied by `-r/--reads-dir`.

```bash
jmphage all \
  -i assemblies/ \
  -r reads/ \
  -d /path/to/jm_db \
  -o results \
  -t 8
```

JMPhage currently limits each batch to a maximum of **5 phages per run**. Split larger collections into multiple runs.

## Input Requirements

JMPhage uses strict naming rules so that FASTA files, read files, and downstream tool outputs remain synchronized.

### File Names

| File type | Required suffix | Example |
|:--|:--|:--|
| Genome assembly | `.fasta` | `myphage01.fasta` |
| Forward reads | `.R1.fastq` | `myphage01.R1.fastq` |
| Reverse reads | `.R2.fastq` | `myphage01.R2.fastq` |

Suffixes are case-sensitive. Files ending in `.fa`, `.fna`, `.fq`, `.fastq.gz`, or `.FASTA` are not recognized by the current strict-mode parser.

### Name Consistency

For each phage, these identifiers must match:

1. FASTA file name without the `.fasta` suffix
2. Paired-read prefix before `.R1.fastq` and `.R2.fastq`
3. Primary FASTA header ID immediately after `>`

Correct example:

```text
File:          myphage01.fasta
FASTA header:  >myphage01
Reads:         myphage01.R1.fastq / myphage01.R2.fastq
```

Incorrect examples:

```text
File:          myphage_v1.fasta
FASTA header:  >Sequence_1
Reads:         myphage.R1.fastq / myphage.R2.fastq
```

```text
File:          Phage_A.fasta
FASTA header:  >Phage_A genome assembled by SPAdes
```

The FASTA header should contain only the phage name, with no spaces or trailing description.

### Mixed Batches

JMPhage supports mixed batches containing newly sequenced isolates and reference genomes. In `all` mode, genomes with paired reads run through mapping, while genomes without paired reads skip the mapping module and continue to annotation and characterization. In batch runs, a per-phage status report is written to `results/summary/pipeline_report.tsv`.

## Usage

JMPhage provides six subcommands.

| Subcommand | Purpose | Reads required |
|:--|:--|:--:|
| `install` | Install Conda environments and/or `jm_db` | No |
| `mapping` | Read mapping, depth profile, and PhageTerm analysis | Yes |
| `annotation` | ORF prediction, functional annotation, and genome map | No |
| `characterization` | Network, collinearity, ANI, and phylogenetic analyses | No |
| `no_mapping` | Run `annotation` then `characterization` | No |
| `all` | Run `mapping`, `annotation`, and `characterization` | Yes (optional in batch mode) |

Common arguments for analysis subcommands:

| Argument | Required | Description |
|:--|:--:|:--|
| `-i, --input` | Yes | Input `.fasta` file or directory containing `.fasta` files |
| `-d, --database` | Yes | Path to installed `jm_db` |
| `-o, --output` | Yes | Output directory |
| `-t, --threads` | No | Number of threads; default `4` |

### `mapping`

Runs Bowtie2 and SAMtools to map paired-end reads to the input genome, calculates per-base depth, plots a coverage profile, and runs PhageTerm.

```bash
jmphage mapping \
  -i my_phage.fasta \
  -1 my_phage.R1.fastq \
  -2 my_phage.R2.fastq \
  -d /path/to/jm_db \
  -o results \
  -t 8
```

For directory input, reads are auto-paired by name:

```bash
jmphage mapping -i input_dir/ -d /path/to/jm_db -o results -t 8
```

With reads in a separate directory:

```bash
jmphage mapping -i assemblies/ -r reads/ -d /path/to/jm_db -o results -t 8
```

### `annotation`

Runs Pharokka for ORF prediction, annotates proteins using VOG-based searches, and produces a circular genome map.

```bash
jmphage annotation \
  -i my_phage.fasta \
  -d /path/to/jm_db \
  -o results \
  -t 8
```

### `characterization`

Builds a protein-sharing network, selects related reference genomes, visualizes genome collinearity, computes VIRIDIC intergenomic similarity, and constructs GLUVAB phylogenetic trees.

```bash
jmphage characterization \
  -i my_phage.fasta \
  -d /path/to/jm_db \
  -o results \
  -t 8
```

For directory input, JMPhage uses joint characterization mode: it builds a shared vConTACT2 network for all input phages and then performs per-phage collinearity, ANI, and tree analyses, plus joint ANI and joint phylogenetic analyses in `results/summary/`.

### `no_mapping`

Runs annotation and characterization without read mapping:

```bash
jmphage no_mapping \
  -i my_phage.fasta \
  -d /path/to/jm_db \
  -o results \
  -t 8
```

### `all`

Runs the complete workflow:

```bash
jmphage all \
  -i my_phage.fasta \
  -1 my_phage.R1.fastq \
  -2 my_phage.R2.fastq \
  -d /path/to/jm_db \
  -o results \
  -t 8
```

Optional PhageTerm reorientation:

```bash
jmphage all \
  -i my_phage.fasta \
  -1 my_phage.R1.fastq \
  -2 my_phage.R2.fastq \
  -d /path/to/jm_db \
  -o results \
  -t 8 \
  --use-phageterm
```

## Output Structure

JMPhage organizes outputs differently for single-phage and multi-phage runs.

- **Single-phage mode** writes all module outputs directly under the output directory.
- **Multi-phage (joint) mode** creates a subdirectory for each input phage, plus a shared `summary/` directory containing the joint vConTACT2 network, joint ANI heatmap, and joint phylogenetic tree.

### Single-phage Output

```text
output/
├── 1.read_mapping/                ← only present if mapping was run
├── 2.ORF_prediction/
├── 3.function_annotation/
├── 4.genome_visualization/
├── 5.shared_network/
├── 6.collinearity_analysis/
├── 7.ANI_analysis/
├── 8.phylogenetic_tree/
└── logs/
```

### Multi-phage Output

```text
output/
├── _joint_manifest.json           ← internal manifest of the joint run
├── <phage_1>/
│   ├── 1.read_mapping/            ← only present if reads were provided
│   ├── 2.ORF_prediction/
│   ├── 3.function_annotation/
│   ├── 4.genome_visualization/
│   ├── 5.shared_network/          ← per-phage handle to the shared network
│   ├── 6.collinearity_analysis/
│   ├── 7.ANI_analysis/            ← per-phage ANI vs its top neighbors
│   ├── 8.phylogenetic_tree/       ← per-phage small tree only
│   └── logs/
├── <phage_2>/
│   └── ...
└── summary/
    ├── pipeline_report.tsv        ← per-phage status table
    ├── shared_network/            ← joint vConTACT2 network (built once)
    ├── ANI_analysis/              ← joint ANI heatmap across all inputs and their neighbors
    ├── phylogenetic_tree/         ← joint small tree, joint big tree, and joint ggtree visualization
    └── logs/
```

In multi-phage mode, the **big phylogenetic tree** is computed only once at the `summary/` level (using all input phages plus the genus-balanced backbone). Each per-phage `8.phylogenetic_tree/` directory contains only its own small tree to keep runtime manageable.

### Module Outputs

The contents of each numbered subdirectory below are the same in both single-phage and multi-phage modes, with the few exceptions noted.

#### `1.read_mapping/`

```text
1.read_mapping/
├── <phage>.sorted.bam
├── base_depth.tsv
├── depth_plot.pdf
├── depth_plot.png
├── <phage>_phageterm.fasta        ← present only when PhageTerm identifies a fixed terminus
└── phageterm/
    ├── <phage>_PhageTerm_report.pdf
    ├── <phage>_sequence.fasta
    └── <phage>_direct-term-repeats.fasta
```

#### `2.ORF_prediction/`

```text
2.ORF_prediction/
├── combined_output.faa            ← predicted proteins (cleaned headers)
└── pharokka_output/               ← full Pharokka native output (GFF, GBK, tsv, etc.)
```

#### `3.function_annotation/`

```text
3.function_annotation/
├── <phage>_vog.tsv                ← DIAMOND blastp full hits against VOG
└── <phage>_top_vog.tsv            ← best VOG hit per protein
```

#### `4.genome_visualization/`

```text
4.genome_visualization/
├── genome_features_classified.tsv ← per-ORF: ID, start, end, strand, functional category
├── protein_function_summary.tsv   ← per-protein: VOG ID + functional annotation
└── genome_plot.pdf                ← circular genome map colored by functional category
```

The genome map uses a fixed color scheme for ten functional categories: structure, DNA replication and repair, transcriptional regulation, lysogeny/lysis, host interaction, metabolism, packaging, hypothetical protein, other, and unclassified.

#### `5.shared_network/`

In single-phage mode this directory holds the complete vConTACT2 input and output:

```text
5.shared_network/
├── pharokka_output/               ← Pharokka output for the input phage
├── all.faa                        ← merged proteins (input + reference database)
├── out_map.csv                    ← merged gene-to-genome mapping
└── vcontact2_output/
    └── c1.ntw                     ← protein-sharing network (used by downstream steps)
```

In multi-phage mode, each phage's `5.shared_network/` contains only a copy of `c1.ntw` together with the per-phage Pharokka output. The full `all.faa` and `out_map.csv` are stored once under `summary/shared_network/`.

#### `6.collinearity_analysis/`

```text
6.collinearity_analysis/
├── <phage>.html                   ← interactive Clinker visualization (query vs top related genomes)
├── gene_functions.csv             ← gene functional annotations used by Clinker
├── gbk/                           ← GenBank files of the query and selected references
└── faa/                           ← protein FASTA files used during annotation
```

#### `7.ANI_analysis/`

```text
7.ANI_analysis/
├── ANI.fasta                      ← merged FASTA of the query and selected neighbors
└── ANI_result/
    └── 04_VIRIDIC_out/
        ├── Heatmap.PDF            ← VIRIDIC heatmap
        └── sim_MA_genCol.csv      ← pairwise intergenomic similarity matrix
```

#### `8.phylogenetic_tree/`

In single-phage mode, both small and big trees are computed:

```text
8.phylogenetic_tree/
├── <phage>_tree.pdf               ← ggtree visualization (based on the small tree)
├── small_tree_in.fasta            ← input FASTA used for the small tree
├── big_tree_in.fasta              ← input FASTA used for the big tree
├── small_tree/
│   └── <phage>_Tree_With_node_IDs.newick
└── big_tree/
    └── <phage>_Tree_With_node_IDs.newick
```

In multi-phage mode, each per-phage `8.phylogenetic_tree/` directory contains only `small_tree/` (no big tree). The joint big tree is generated once under `summary/phylogenetic_tree/big_tree/`.

#### `summary/` (multi-phage mode only)

```text
summary/
├── pipeline_report.tsv            ← per-phage status: name, FASTA, R1/R2, has_reads, mapping_status
├── shared_network/
│   ├── pharokka/<phage_name>/     ← Pharokka outputs for each input phage
│   ├── all.faa
│   ├── out_map.csv
│   └── vcontact2_output/c1.ntw    ← joint network across all input phages
├── ANI_analysis/
│   ├── ANI.fasta
│   └── ANI_result/04_VIRIDIC_out/
│       ├── Heatmap.PDF            ← joint VIRIDIC heatmap
│       └── sim_MA_genCol.csv
├── phylogenetic_tree/
│   ├── joint_tree.pdf             ← joint ggtree visualization
│   ├── small_tree_in.fasta
│   ├── big_tree_in.fasta
│   ├── small_tree/
│   │   └── joint_Tree_With_node_IDs.newick
│   └── big_tree/
│       └── joint_Tree_With_node_IDs.newick
└── logs/
```

## JMPhage Database

The installed `jm_db` directory contains reference data, scripts, and bundled third-party resources used by the pipeline:

| Database component | Used by | Purpose |
|:--|:--|:--|
| `mapping_profile/` | `mapping` | Depth-profile plotting and PhageTerm wrapper resources |
| `pharokka_db/` | `annotation`, `characterization` | Pharokka database for phage gene prediction and annotation |
| `VOG_HMM/` | `annotation`, `characterization` | VOG DIAMOND/HMMER annotation and functional categories |
| `coline_profile/` | `characterization` | Reference proteins/genomes and scripts for vConTACT2 and Clinker analyses |
| `VIRIDIC/` | `characterization` | Stand-alone VIRIDIC scripts |
| `GLUVAB_tree/` | `characterization` | GLUVAB scripts, backbone references, taxonomy table, and tree plotting resources |
| `cluster_one-1.0.jar` | `characterization` | ClusterONE dependency used by vConTACT2 |

The database is versioned separately from the Python package and should be cited via its Zenodo DOI when used in published analyses.

## Integrated Software and References

JMPhage is a wrapper pipeline. It does not replace the original software packages it integrates; instead, it coordinates them in a reproducible workflow for tailed-phage analysis. If you use JMPhage in a publication, please cite JMPhage and the original tools that contributed to your reported results.

### Core Tools Used by JMPhage

| Module | Tool | Role | Reference |
|:--|:--|:--|:--|
| `mapping` | Bowtie2 | Read alignment | Langmead & Salzberg, 2012, *Nat Methods*, [10.1038/nmeth.1923](https://doi.org/10.1038/nmeth.1923) |
| `mapping` | SAMtools | BAM processing & depth calculation | Danecek et al., 2021, *GigaScience*, [10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008) |
| `mapping` | PhageTerm | Terminus & packaging mechanism detection | Garneau et al., 2017, *Sci Rep*, [10.1038/s41598-017-07910-5](https://doi.org/10.1038/s41598-017-07910-5) |
| `annotation` / `characterization` | Pharokka | ORF prediction & phage annotation | Bouras et al., 2023, *Bioinformatics*, [10.1093/bioinformatics/btac776](https://doi.org/10.1093/bioinformatics/btac776) |
| `annotation` / `characterization` | DIAMOND | Protein similarity search | Buchfink et al., 2015, *Nat Methods*, [10.1038/nmeth.3176](https://doi.org/10.1038/nmeth.3176) |
| `annotation` / `characterization` | HMMER | HMM-based functional classification | http://hmmer.org |
| `annotation` / `characterization` | VOG | HMM profile database | Trgovec-Greif et al., 2024, *Viruses*, [10.3390/v16081191](https://doi.org/10.3390/v16081191) |
| `characterization` | vConTACT2 | Protein-sharing network construction | Bin Jang et al., 2019, *Nat Biotechnol*, [10.1038/s41587-019-0100-8](https://doi.org/10.1038/s41587-019-0100-8) |
| `characterization` | Clinker | Gene collinearity visualization | Gilchrist & Chooi, 2021, *Bioinformatics*, [10.1093/bioinformatics/btab007](https://doi.org/10.1093/bioinformatics/btab007) |
| `characterization` | VIRIDIC | Intergenomic similarity (ANI) heatmap | Moraru et al., 2020, *Viruses*, [10.3390/v12111268](https://doi.org/10.3390/v12111268) |
| `characterization` | GLUVAB | Whole-genome phylogenetic tree construction | Coutinho et al., 2019, *BMC Biol*, [10.1186/s12915-019-0723-8](https://doi.org/10.1186/s12915-019-0723-8) |
| `characterization` | ggtree | Tree visualization in R | Yu et al., 2017, *Methods Ecol Evol*, [10.1111/2041-210X.12628](https://doi.org/10.1111/2041-210X.12628) |


## Re-running and Resuming

JMPhage is restart-aware. If key output files are already present, the corresponding step is skipped on re-run.

- `mapping` skips the read-mapping pipeline when a depth-profile PDF exists in `1.read_mapping/`. PhageTerm will still run if it has not been executed before.
- `annotation` skips Pharokka if `2.ORF_prediction/combined_output.faa` already exists.
- `annotation` skips HMMscan if `protein_function_summary.tsv` and `genome_features_classified.tsv` already exist, and only regenerates the genome plot.
- `characterization` skips vConTACT2, VIRIDIC, tree, and Clinker steps when their expected output files are present.

To force a step to run again, remove the corresponding output folder before re-running JMPhage.

## Troubleshooting

- **`conda` command not found**: install Miniconda/Anaconda and make sure Conda is available in your shell.
- **Input FASTA not recognized**: rename files to end exactly with `.fasta`.
- **Reads not found**: check that read files are named `<phage>.R1.fastq` and `<phage>.R2.fastq`, and use `-r/--reads-dir` if reads are stored separately.
- **FASTA header mismatch**: ensure the first token after `>` matches the FASTA file name without `.fasta`.
- **Database file missing**: verify that `jm_db` was fully extracted and contains `pharokka_db/`, `VOG_HMM/`, `coline_profile/`, `VIRIDIC/`, `GLUVAB_tree/`, `mapping_profile/`, and `cluster_one-1.0.jar`.
- **Batch size exceeded**: JMPhage limits each run to five phages; split larger collections into multiple runs.

## Citation

If JMPhage is useful in your work, please cite the JMPhage software and database. A manuscript citation will be added here when available.

JMPhage integrates multiple third-party tools. Please also cite the original software papers listed in [Integrated Software and References](#integrated-software-and-references), according to the modules and results used in your analysis.

## Data and Code Availability

The JMPhage source code is available from this GitHub repository. The accompanying `jm_db` reference database is available from Zenodo at [10.5281/zenodo.20838733](https://doi.org/10.5281/zenodo.20838733). Example data, benchmarking datasets, and manuscript source data should be deposited with stable identifiers before publication.

## License

JMPhage is distributed under the GNU General Public License v3.0. See the [GPL-3.0 license](https://www.gnu.org/licenses/gpl-3.0.en.html) for details.

## Contact

Please use the GitHub issue tracker for bug reports, installation problems, feature requests, and questions about JMPhage output interpretation.
