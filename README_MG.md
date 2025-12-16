# 🏭 ModGen

**ModGen** streamlines HPC software management by automatically converting apps and genome references into modules. It allows users to access Conda built and non-Conda software (e.g., 10X Cell Ranger, Illumina ORAD) but also genome references via standard `module load` commands without manual configuration.

- **Auto-Generation**: Installs packages and builds Lua/Tcl modulefiles automatically.
- **Smart Dependencies**: Scans scripts to detect and install missing modules. (`#DEP:` and `module load`)
- **Native Integration**: seamless `module avail` and `module load` experience.
- **Context-Aware**: Manages references/indexes and displays helpful info when not in a SLURM job.
- **SLURM Integration**: Submits index generation jobs automatically when needed.

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
ml grch38/cellranger/2024-A
# Available environment variables:
#   CELLRANGER_REF_DIR (cellranger reference dir)
#   GENOME_FASTA       (genome fasta)
#   ANNOTATION_GTF_GZ  (10X modified gtf)
#   STAR_INDEX_DIR     (STAR index dir)
```

[Read the full ModGen Manual](assets/MANUAL_MG.md)

> [!NOTE]
> Make sure to initialize shell by running `modgen init` once. Then source your shell configuration file (e.g., `source ~/.bashrc`).

## 📂 Naming Convention

ModGen utilizes a standardized naming schema to organize software and reference data.

- Apps: `name/version` (e.g., `bcftools/1.16`)
- References: `assembly/datatype/version`
  - `grcm39/genome/gencode`: mouse genome GRCm39 with Gencode style naming
  - `grcm39/salmon/1.10.2/gencodeM33`: salmon index for mouse genome GRCm39 with Gencode M33 annotation

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

## Acknowledgements

This project used computational resources provided by

<div>
<a href="https://www.mcmaster.ca" target="_blank" rel="noopener noreferrer"><img src="assets/logo_mcmaster.svg" height="45" alt="McMaster University Logo"></a>
&nbsp;&nbsp;&nbsp;
<a href="https://alliancecan.ca" target="_blank" rel="noopener noreferrer"><img src="assets/logo_alliance.svg" height="40" alt="Digital Research Alliance of Canada Logo"></a>
</div>
