# Use ModGen to Manage System Level Modules

🏭 **ModGen** is a tool designed to automate the creation and management of **Lmod/Environment Modules** on HPC systems. 

It is the ideal choice for users who want to take full advantage of tools available on the HPC system rather. After the installation, you can use `module load` from Lmod/Environment Modules.

Please read [Concepts](concepts.md) before proceeding.

## ⚡ Quick Look

```bash
# 1. Initialize shell hooks (only needed once)
modgen init

# 2. Install samtools and generate the modulefile
modgen create samtools/1.16

# 3. Load the module (standard HPC way)
module load samtools/1.16
samtools --version

# 4. Search for packages to install
modgen avail salmon grcm M36

# 5. Auto install dependencies from a script
modgen check analysis.sh -a

# 6. Helpful info example: Cellranger reference
module load grch38/cellranger/2024-A
# Available environment variables:
#   CELLRANGER_REF_DIR (cellranger reference dir)
#   GENOME_FASTA       (genome fasta)
#   ANNOTATION_GTF_GZ  (10X modified gtf)
#   STAR_INDEX_DIR     (STAR index dir)
```

## 🚀 Dependencies Automation

ModGen can recognize both `#DEP:` tags and `module load` commands to automatically install required modules.

Example Script (`analysis.sh`):

```bash
#!/bin/bash
module purge
module load salmon/1.10.2
ml grcm39/salmon/1.10.2/gencodeM33

salmon quant -i $SALMON_INDEX_DIR ...
```

Check and auto install dependencies:

```bash
# Print modules and install status
modgen check analysis.sh

# Auto install missing modules
modgen check analysis.sh -a
```

## 🤖 SLURM Automation

When you request a reference or an environment that requires significant computation to prepare, **ModGen will automatically submit SLURM jobs** to handle the heavy lifting for you.

## 🧫 Case Study: Cellranger Count

### 📜 Count Script

The following is an example SLURM script for running `cellranger count` using the cellranger modules.

```bash
#!/bin/bash
#SBATCH --job-name=cellranger-quant
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB

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

### 📥 Install required modules

You can check the dependencies and auto install them using:

```bash
modgen check cellranger_quant.sh -a
```

or you can explicitly create the modules using:

```bash
modgen create cellranger/9.0.1 grch38/cellranger/2024-A
```

Since the download link for cellranger is only valid for one day, you will be prompted to provide the download link during the build process.

```
[ModGen][NOTE] Build script requires input: ⚠️ 10X links only valid for one day. Please go to the link below and get tar.gz link.
https://www.10xgenomics.com/support/software/cell-ranger/downloads/previous-versions
Enter here: 
```

You need to paste the valid download link and press Enter to continue the build.

Since cellranger references are prebuilt, **ModGen** will download the extract the reference files and generate the modulefile on the login node.

### 📤 Load and use modules

Then you can submit the script as usual.

```bash
sbatch cellranger_quant.sh
```

## 🔗 Related Resources

- [ModGen Installation](../installation/modgen.md)
