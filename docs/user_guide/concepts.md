# Concepts: 🧩 Overlays

**Modular, Reproducible, and Self-Contained Analysis Environments.**

Reproducible analyses require careful management of both software tools and reference data. **CondaTainer** addresses this challenge through its overlay system.

## 💡 The Main Idea

Overlays are stackable, self-contained files that encapsulate executables and data. When loaded, **CondaTainer** modifies environment variables (such as `PATH`) and establishes path bindings to make software and data accessible.

## 📦 Overlay Types

| Term | Ext | R/W | Content Path | Purpose |
|------|-----------|-----|--------------|-------------|
| OS | `.sqf` | R/O | `/bin`, `/lib` etc. | System Foundation - Run standalone. |
| Module | `.sqf` | R/O | `/cnt/name/version` | Individual Tool - Run on top of OS. |
| Bundle | `.sqf` | R/O | `/cnt/<env_name>` | Frozen Conda Env - Run on top of OS. |
| Workspace | `.img` | R/W | `/ext3/env` | Writable Conda Env - Run on top of OS. |

```{note}
- `sqf` cannot be renamed. The name is used for the mount path.
- `img` can be renamed freely.
```

````{tip}
OS overlays are stackable when based on the same distribution version.

```bash
condatainer e rstudio-server build-essential
```
````

You can use `condatainer info` to check overlay types and contents.

### 🖥️ OS Overlays

OS overlays are built from Apptainer definition files. They expose system-level paths (`/bin`, `/lib`, etc.) and can run standalone as a base image.

### 🧩 Module Overlays

Module overlays are categorized into two types: **Apps** and **Data**.

- **Apps**: Software packages, binaries, pipelines. Mount at `/cnt/<name>/<version>`.
- **Data**: Project data, genome indices, annotations. Mount at `/cnt/<assembly|project>/<datatype>/...`.

Examples:
- App: `cellranger/9.0.1`
- Data: `grch38/salmon/1.10.2/gencode47`

### 📂 Naming Convention

- OS: `<distro>/<name>` (e.g., `ubuntu24/igv`)
- Bundle / Env: `<name>` — no slash (e.g., `env`, `sci_rna`)
- App (Module): `<name>/<version>` (e.g., `cellranger/9.0.1`)
- Data (Module): `<assembly|project>/<datatype>/<version>`
  - `grcm39/genome/gencode`: GRCm39 genome with Gencode style naming
  - `grcm39/salmon/1.10.2/gencodeM33`: salmon index for GRCm39 Gencode M33 transcripts

#### Version Delimiters

The following delimiters are accepted for version specification: `/`, `--`, `=`, `@`

**Important**: Because `--` serves as a delimiter, it should not be used within overlay names themselves.

## 🧱 Stacking Overlays

Overlays can be stacked in a specific order to create a layered environment.

- **Base Image**: Apptainer image (e.g. `ubuntu24--base_image.sif`) or OS overlay
- **OS Overlay**: Provides system-level libraries and tools (should have the same distro version as the base image)
- **Module Overlays**: Individual software packages and data (e.g. `cellranger/9.0.1`, `grch38/cellranger/2024-A`)
- **Bundle Overlay**: Frozen conda environment (e.g. `env.sqf`)
- **Workspace Overlay**: Writable conda environment (e.g. `env.img`)

```{note}
The later overlays will overwrite the previous ones if there are conflicts. For example, if both the base image and an OS overlay contain `/bin/bash`, the one from the OS overlay will be used.

**CondaTiainer** will always put workspace overlays on the top.
```

```{note}
Module and bundle overlays mounted later will appear earlier in the `$PATH` (i.e., be prepended).

Data overlays will not be added to the `$PATH`, but you can access them through the mount path (e.g., `/cnt/grch38/salmon/1.10.2/gencode47`). They have custom environment variables as well, use `info <overlay>` to check.
```

### Stacking Example: Override cutadapt version in trim-galore

By default `trim_galore` will come with the latest version of `cutadapt`, (5.2 as of Feb 2026). But you want to use `cutadapt` 5.0, to replicate a previous analysis.

By stacking the `cutadapt/5.0` module overlay on top of the `trim_galore` module overlay, you can easily override the version of `cutadapt` without rebuilding the entire `trim_galore` conda environment.

```bash
condatainer e trim-galore/0.6.11 -- trim_galore
# Multicore support not enabled. Proceeding with single-core trimming.
# Path to Cutadapt set as: 'cutadapt' (default)
# Cutadapt seems to be working fine (tested command 'cutadapt --version')
# Cutadapt version: 5.2
condatainer e trim-galore/0.6.11 -- which cutadapt
# /cnt/trim-galore/0.6.11/bin/cutadapt
```

```bash
condatainer e trim-galore/0.6.11 cutadapt/5.0 -- trim_galore
# ...
# Cutadapt version: 5.0
condatainer e trim-galore/0.6.11 cutadapt/5.0 -- which cutadapt
# /cnt/cutadapt/5.0/bin/cutadapt
condatainer e trim-galore/0.6.11 cutadapt/5.0 -- bash -c 'echo $PATH'
# ...:/cnt/cutadapt/5.0/bin:/cnt/trim-galore/0.6.11/bin:/usr/local/...

# Change to a different order:
condatainer e cutadapt/5.0 trim-galore/0.6.11 -- bash -c 'echo $PATH'
# ...:/cnt/trim-galore/0.6.11/bin:/cnt/cutadapt/5.0/bin:/usr/local/...
```

## 🐍 Leveraging Conda Resources

**CondaTainer** leverages the extensive `conda-forge` and `bioconda` ecosystems, which provide most bioinformatics software as conda packages. **CondaTainer** will automatically create module overlays for these packages.

For software unavailable through conda, custom build scripts can be created to download and install the software.

**Examples**: 10X [cellranger/9.0.1](https://github.com/Justype/cnt-scripts/blob/main/build-scripts/cellranger/9.0.1) and Illumina [orad/2.7.0](https://github.com/Justype/cnt-scripts/blob/main/build-scripts/orad/2.7.0)

## 🔄 Module Workflow

- [Manage Module Overlays](./module_overlays.md)
- [Bundle Overlays: Read-only project environment](./bundle_overlays.md)
- [Workspace Overlays: Writable project environment](./workspace_overlays.md)
- [Fakeroot](../qa/overlayfs.md#fakeroot)
