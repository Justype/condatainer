# 🧩 Concepts: Modules

**Modular, Reproducible, and Self-Contained Analysis Environments.**

To make analyses reproducible, we use a **modular architecture** to manage software and reference data. This approach allows us to dynamically load tools and genomes without worrying about installation conflicts, file paths, or complex dependencies.

## 💡 The Main Idea

Modules are self-contained packages. When loaded, both **Lmod** and **CondaTainer** can modify environment variables (like `PATH`, `LD_LIBRARY_PATH`, etc.) to make software and data easily accessible.

## 📦 Module Types

We distinguish between logic (app) and data (references) to ensure flexibility and reproducibility.

By combining **App** modules with **Ref** modules, we create a complete, version-controlled environment. This guarantees that the software version matches the reference data structure, preventing common analysis errors.

### 🛠️ 1. App Modules (Logic)

App modules encapsulate the software application, pipeline, or binaries.

- **Role**: Provides the executable tools and manages their runtime dependencies.
- **Action**: Updates system PATH so we can run commands (e.g., cellranger) directly.

**Example**: `salmon/1.10.2`

### 🧬 2. Ref Modules (Data)

Ref (Reference) modules provide version-specific genome indices, annotations, or reference files.

- **Role**: Ensures the analysis uses the exact reference build required by the software version.
- **Action**: Sets specific environment variables pointing to the data or the directories, removing the need to hardcode paths.

**Example**: `grch38/salmon/1.10.2/gencode47`

### 📂 Naming Convention

Both tools utilize a standardized naming schema to organize software and reference data.

- System Apps: `name` (e.g., `rstudio-server`, created by `name` or `-n name` flag) (**CondaTainer** only)
- Project Env: `name` (e.g., `env`, `sci_rna`, created by `-p prefix` flag) (**CondaTainer** only)
- Apps: `name/version` (e.g., `bcftools/1.16`)
- References: `assembly/datatype/version`
  - `grcm39/genome/gencode`: GRCm39 genome with Gencode style naming
  - `grcm39/salmon/1.10.2/gencodeM33`: salmon index for GRCm39 Gencode M33 transcripts


## 🐍 Leveraging Conda Resources

Thanks to `conda-forge` and `bioconda`, most bioinformatics software can be easily installed with conda. Both **ModGen** and **CondaTainer** leverage conda to create isolated environments for each app module, ensuring that dependencies are managed without conflicts.

If a required package is not available in conda, custom build scripts can be created to compile and install the software from source.

e.g. 10X `cellranger` and Illumina `orad`

## 🔄 Workflow

### 📥 Install

Use **ModGen** or **CondaTainer** to install the desired app and reference modules.

```bash
condatainer install grch38/salmon/1.10.2/gencode47
# It will install
# - salmon 1.10.2
# - grch38 genome fasta
# - gencode47 transcript fasta
# - Use these to build salmon index (SLURM job auto submitted)
```

Same as the **ModGen** way:

```bash
modgen install grch38/salmon/1.10.2/gencode47
```

### 🚀 Load and Use

Load the required app and ref modules and run your analysis.

**CondaTainer** way:

```bash
condatainer exec \
  -o salmon/1.10.2 \
  -o grch38/salmon/1.10.2/gencode47 \
  salmon quant -i $SALMON_INDEX_DIR ...
```

Environment Modules (Lmod) way:

```bash
module load salmon/1.10.2
module load grch38/salmon/1.10.2/gencode47
salmon quant -i $SALMON_INDEX_DIR ...
```
