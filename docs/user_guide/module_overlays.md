# Creating and Using Module Overlays

📦 **CondaTainer** allows you to create [module overlays](./concepts.md#-overlay-types) for your project.

They are stackable, read-only, and highly compressed overlays that contain pre-installed modules (apps and references), which is ideal for HPC environments where inode usage is a concern.

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
# [CondaTainer] Overlay envs:
#   CELLRANGER_REF_DIR: cellranger reference dir
#   GENOME_FASTA      : genome fasta
#   ANNOTATION_GTF_GZ : 10X modified gtf
#   STAR_INDEX_DIR    : STAR index dir
```

## 🚀 Dependencies Automation

### 🧬 Reference Overlay Installation

To Install a Salmon Index Overlay. You don't need to:

- Load modules or install Salmon manually.
- Manually download genome FASTA and transcript FASTA files.
- Create decoy FASTA.
- Build the Salmon index and submit SLURM jobs.

**Condatainer** will handle all these steps for you automatically!

```bash
condatainer install grch38/salmon/1.10.2/gencode47
# This command will:
# - Install Salmon 1.10.2 via bioconda
# - Download GRCh38 genome FASTA
# - Download Gencode 47 transcript FASTA
# - Build Salmon decoy index (SLURM job automatically submitted)
```

### 🏷️ Declaring Dependencies in Scripts

**CondaTainer** can recognize both `#DEP:` tags and `module load` commands to automatically install required modules.

Yes, it can read `module load` and `ml` commands inside bash scripts!

Example Script (`analysis.sh`):

```bash
#!/bin/bash
#DEP: salmon/1.10.2
#DEP: grcm39/salmon/1.10.2/gencodeM33

salmon quant -i $SALMON_INDEX_DIR ...
```

It is also compatible with `module load` and `ml` commands:

```bash
#!/bin/bash
module purge
module load salmon/1.10.2
ml grcm39/salmon/1.10.2/gencodeM33

salmon quant -i $SALMON_INDEX_DIR ...
```

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

## 🤖 SLURM Automation

When you request a reference or an environment that requires significant computation to prepare, **CondaTainer will automatically submit SLURM jobs** to handle the heavy lifting for you.

You can follow the example below to define SLURM parameters in your script and submit it to SLURM as usual.

```bash
#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#DEP:samtools/1.22.1

if [ -z "$IN_CONDATAINER" ] && command -v condatainer >/dev/null 2>&1; then
    condatainer run "$0" "$@"
    exit $?
fi

samtools --version
```

Then you can

```bash
sbatch analysis.sh
```

## 🧫 Case Study: Cellranger Count

### 📜 Count Script

The following is an example SLURM script for running `cellranger count` using the cellranger overlays.

```bash
#!/bin/bash
#SBATCH --job-name=cellranger-quant
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB

if [ -z "$IN_CONDATAINER" ] && command -v condatainer >/dev/null 2>&1; then
    condatainer run "$0" "$@"
    exit $?
fi

module purge
module load cellranger/9.0.1
module load grch38/cellranger/2024-A

cellranger count --id=sample1 \
  --transcriptome=$CELLRANGER_REF_DIR \
  --fastqs=/path/to/fastqs \
  --sample=sample1 \
  --localcores=$SLURM_CPUS_PER_TASK \
  --localmem=$((SLURM_MEM_PER_NODE/1024))
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

Then you can submit the script as usual.

```bash
sbatch cellranger_quant.sh
```

## 🔗 Related Resources

- [Installation](./installation.md)
- [CondaTainer Manual](../manuals/condatainer.md)
- [Workspace Overlays: Writable Project-Level](./workspace_overlays.md)
- [Bundle Overlays: Read-Only Project-Level](./bundle_overlays.md)

Tutorials using CondaTainer:

- [RStudio Server on HPC](../tutorials/rstudio-server_on_HPC.md)
- [code-server on HPC](../tutorials/code-server_on_HPC.md)
- [VS Code Tunnel on HPC](../tutorials/vscode-tunnel_on_HPC.md)
