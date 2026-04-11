# Creating and Using Module Overlays

📦 **CondaTainer** allows you to create [module overlays](./concepts.md#-overlay-types) for your project.

They are stackable, read-only, and highly compressed overlays that contain app or reference data, which is ideal for HPC environments where inode usage is a concern.

Please read [Concepts](concepts.md) before proceeding.

## 👀 Quick Look

```bash
condatainer avail # List available recipes
condatainer list  # List installed

# Install and Run
condatainer create samtools/1.16
condatainer exec -o samtools/1.16 samtools --version

# Create a read-only project overlay from a YAML file
condatainer create -f environment.yml -p my_analysis
condatainer exec -o my_analysis.sqf bash

# Automatically install the dependencies
condatainer check analysis.sh -a

# Print helpful info when running with overlays
condatainer exec -o grch38/cellranger/2024-A bash
# [CNT] Overlay envs:
#   CELLRANGER_REF_DIR: cellranger reference dir
#   GENOME_FASTA      : genome fasta
#   ANNOTATION_GTF_GZ : 10X modified gtf
#   STAR_INDEX_DIR    : STAR index dir
```

## 📥 Installing Module Overlays

Use `condatainer avail` to browse available build scripts, then `condatainer create` to install.

### Normal Build Script Module

A normal build script targets a fixed name and version. Install it directly:

```bash
condatainer avail cellranger           # search available scripts
condatainer create cellranger/9.0.1    # app module
condatainer create grch38/genome/gencode  # data module
```

### Template Script Module

A template script covers many versions through placeholders. It appears in `avail` output as a collapsed group:

```
grch38/salmon-gencode  [486 variants]
  → grch38/salmon/{salmon_version}/gencode{gencode_version}
  - salmon_version:   1.0.0-1.11.4  (18 values)
  - gencode_version:  23-49  (27 values)
```

There are two ways to install:

**1. Use the template name** — **CondaTainer** prompts for each placeholder interactively:

```bash
condatainer create grch38/salmon-gencode
# [CNT] Placeholder template: grch38/salmon-gencode
# [CNT] Salmon GRCh38 GENCODE{gencode_version} index for transcript quantification
# Target: grch38/salmon/{salmon_version}/gencode{gencode_version}
#   salmon_version [1.0.0-1.11.4] (default: 1.11.4): 1.10.2
#   gencode_version [23-49] (default: 49):
#   → Creating grch38/salmon/1.10.2/gencode49
```

**2. Specify the target directly** — skip the prompts by providing the resolved name:

```bash
condatainer create grch38/salmon/1.10.2/gencode47
```

```{tip}
When entering placeholder values, you can hit <Tab> to see the available options and autocomplete.
```

## 🚀 Dependencies Automation

### 🧬 Data Overlay Installation

To Install a Salmon Index Overlay. You don't need to:

- Load modules or install Salmon manually.
- Manually download genome FASTA and transcript FASTA files.
- Create decoy FASTA.
- Build the Salmon index and submit scheduler jobs.

**Condatainer** will handle all these steps for you automatically!

```bash
condatainer create grch38/salmon/1.10.2/gencode47
# This command will:
# - Create Salmon 1.10.2 module overlay
# - Download GRCh38 genome FASTA as overlay
# - Download Gencode 47 transcript FASTA as another overlay
# - Submit scheduler jobs to build the Salmon index using these overlays
```

### 🏷️ Declaring Dependencies in Scripts

Declare dependencies with `#DEP:` tags at the top of your script.

Example Script (`analysis.sh`):

```bash
#!/bin/bash
#DEP: salmon/1.10.2
#DEP: grcm39/salmon/1.10.2/gencodeM33

salmon quant -i $SALMON_INDEX_DIR ...
```

````{note}
**CondaTainer** will automatically set environment variables like `SALMON_INDEX_DIR`.

If you don't know the variable names, you can use `info` to check:

```bash
condatainer info grcm39/transcript-gencode/M9
# Environment
#  - TRANSCRIPT_FASTA=/cnt/grcm39/transcript-gencode/M9/gencode.vM9.transcripts.fa
#    # GRCm39 GENCODE vM9 transcript
```
````

Check and auto install dependencies:

```bash
# Print dependencies and install status
condatainer check analysis.sh

# Auto install missing overlays
condatainer check analysis.sh -a
```

Execute the script with CondaTainer:

```bash
# Run with CondaTainer
condatainer run analysis.sh
```

## 🤖 Scheduler Automation

When you request a reference or an environment that requires significant computation to prepare, **CondaTainer will automatically submit scheduler jobs** (SLURM, PBS, LSF, or HTCondor) to handle the heavy lifting for you.

Example Script (`analysis.sh`):

```bash
#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#DEP:samtools/1.22.1

samtools --version
```

Install missing dependencies first:

```bash
condatainer check -a salmon_quant.sh
```

After all dependencies are installed, submit as a job:

```bash
condatainer run salmon_quant.sh
```

If no scheduler directives are found or job submission is disabled, the script will run immediately in the current shell.

## 🧫 Case Study: Cellranger Count

### 📜 Count Script

The following is an example scheduler script (SLURM) for running `cellranger count` using the cellranger overlays.

```bash
#!/bin/bash
#SBATCH --job-name=cellranger-quant
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB
#DEP: cellranger/9.0.1
#DEP: grch38/cellranger/2024-A

cellranger count --id=sample1 \
  --transcriptome=$CELLRANGER_REF_DIR \
  --fastqs=/path/to/fastqs \
  --sample=sample1 \
  --localcores=$NCPUS \
  --localmem=$MEM_GB
```

### 📥 Install required overlays

You can check the dependencies and automatically install them using:

```bash
condatainer check cellranger_quant.sh -a
```

or you can explicitly create the module overlays using:

```bash
condatainer create cellranger/9.0.1 grch38/cellranger/2024-A
```

Since the download link for cellranger is only valid for one day, you will be prompted to provide the download link during the build process.

```
[CNT][NOTE] Build script requires input: ⚠️ 10X links only valid for one day. Please go to the link below and get tar.gz link.
https://www.10xgenomics.com/support/software/cell-ranger/downloads/previous-versions
Enter here:
```

You need to paste the valid download link and press Enter to continue the build.

Since cellranger references are prebuilt, **CondaTainer** will download and extract the reference files and create overlays on the login node.

### 📤 Load and use overlays

Then you can submit the script using **CondaTainer**.

```bash
condatainer run cellranger_quant.sh
```

## 🔗 Related Resources

- [CondaTainer Manual](../manuals/condatainer.md)
- [Workspace Overlays: Writable Project-Level](./workspace_overlays.md)
- [Bundle Overlays: Read-Only Project-Level](./bundle_overlays.md)
- [Scheduler Integration](../qa/scheduler.md)

Tutorials using CondaTainer:

- [RStudio Server on HPC](../tutorials/rstudio-server_on_HPC.md)
- [VS Code Server on HPC](../tutorials/vscode-server_on_HPC.md)
- [VS Code Tunnel on HPC](../tutorials/vscode-tunnel_on_HPC.md)
