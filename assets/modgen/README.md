# 🏭 ModGen

<img src="./logo_mg.png" height="100" alt="ModGen Logo">

**ModGen** streamlines HPC software management by automatically converting apps and genome references into modules. It allows users to access Conda-built and non-Conda software (e.g., 10X Cell Ranger, Illumina ORAD), as well as genome references, via standard `module load` commands without manual configuration.

- **Auto-Generation**: Installs packages and builds Lua/Tcl modulefiles automatically.
- **Smart Dependencies**: Scans scripts to detect and install missing modules. (`#DEP:` and `module load`)
- **SLURM Integration**: Submit jobs automatically when needed.

## 🛠️ Installation

Run the following command to install **ModGen**:

```bash
curl -fsSL https://raw.githubusercontent.com/Justype/condatainer/main/assets/install.modgen.sh | bash
```

## 👀 Quick Look

```bash
modgen init # Initialize shell hooks (run once)

modgen create samtools/1.16 # create module from Conda package
module load samtools/1.16 # load the module
samtools --version

modgen avail salmon grcm M36 # list available modules
modgen check analysis.sh -a # Auto-install missing dependencies in script

# Helpful info example: Cellranger reference
ml grch38/cellranger/2024-A
# Available environment variables:
#   CELLRANGER_REF_DIR (cellranger reference dir)
#   GENOME_FASTA       (genome fasta)
#   ANNOTATION_GTF_GZ  (10X modified gtf)
#   STAR_INDEX_DIR     (STAR index dir)
```

- 📜 [Read the full ModGen Manual](./manual.md)
- 📁 [Naming Conventions](../../docs/user_guide/concepts.md#-naming-convention)

> [!NOTE]
> Make sure to initialize shell by running `modgen init` once. Then source your shell configuration file (e.g., `source ~/.bashrc`).

## 🚀 Automation

**ModGen** can recognize both `#DEP:` tags and `module load` commands to automatically install required modules. If `#SBATCH` are present, submit to Slurm automatically.

Example Script (`analysis.sh`):

```bash
#!/bin/bash
#SBATCH --job-name=analysis
#SBATCH --output=analysis.out
#SBATCH --cpus-per-task=16
#SBATCH --mem=30G
module purge
module load salmon/1.10.2
ml grcm39/salmon/1.10.2/gencodeM33

salmon quant -i $SALMON_INDEX_DIR ...
```

Check and automatically and Run with auto-install:

```bash
modgen run analysis.sh -a
```

## 🗑️ Uninstallation

If you just want to remove **ModGen**, open your `~/.bashrc` or `~/.zshrc` file and remove the lines related to **ModGen**. Look for and delete the following block:

```bash
# >>> MODGEN MODULES >>>
modgen configs
# <<< MODGEN MODULES <<<
```

And remove the **ModGen** executable and the modules directory.

```bash
INSTALL_DIR="$SCRATCH/condatainer"  # or your custom installation path
rm -rf "$INSTALL_DIR/bin/modgen"
rm -rf "$INSTALL_DIR/{apps,apps-modules,ref,ref-modules}"
```

## Acknowledgements

This project used computational resources provided by

<div>
<a href="https://www.mcmaster.ca" target="_blank" rel="noopener noreferrer"><img src="assets/logo_mcmaster.svg" height="45" alt="McMaster University Logo"></a>
&nbsp;&nbsp;&nbsp;
<a href="https://alliancecan.ca" target="_blank" rel="noopener noreferrer"><img src="assets/logo_alliance.svg" height="40" alt="Digital Research Alliance of Canada Logo"></a>
</div>
