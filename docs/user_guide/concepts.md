# Concepts: 🧩 Overlays

**Modular, Reproducible, and Self-Contained Analysis Environments.**

Reproducible bioinformatics analyses require careful management of both software tools and reference data. **CondaTainer** addresses this challenge through its overlay system.

## 💡 The Main Idea

Overlays are self-contained files that encapsulate executables and data. When loaded, **CondaTainer** modifies environment variables (such as `PATH`) and establishes path bindings to make software and data immediately accessible.

## 📦 Overlay Types

| Term | Ext | R/W | Content Path | Purpose |
|------|-----------|-----|--------------|-------------|
| OS | `.sqf` | R/O | `/bin`, `/lib` etc. | System Foundation - Run standalone. |
| Module | `.sqf` | R/O | `/cnt/name/version` | Individual Tool - Run on top of OS. |
| Bundle | `.sqf` | R/O | `/cnt/<env_name>` | Frozen Conda Env - Run on top of OS. |
| Workspace | `.img` | R/W | `/ext3/env` | Writable Conda Env - Run on top of OS. |

```{note}
- `sqf` files cannot be renamed. **CondaTainer** uses the name to identify the mount path.
- `img` files can be renamed freely.
```

```{tip}
OS overlays are stackable when based on the same distribution version.

Like `condatainer -o rstudio-server -o build-essential`
```

### 🧩 Module Overlays

Module overlays are categorized into two main types: **App** modules and **Ref** modules.

- **App Modules**: Contain the software applications, pipelines, or binaries needed for analysis.
- **Ref Modules**: Contain version-specific genome indices, annotations, or reference files.

Examples:
- App Module: `salmon/1.10.2`
- Ref Module: `grch38/salmon/1.10.2/gencode47`

### 📂 Naming Convention

- System Apps: `name` (e.g., `rstudio-server`)
- Workspace: `name` (e.g., `env`, `sci_rna`)
- Apps: `name/version` (e.g., `bcftools/1.16`)
- References: `assembly/datatype/version`
  - `grcm39/genome/gencode`: GRCm39 genome with Gencode style naming
  - `grcm39/salmon/1.10.2/gencodeM33`: salmon index for GRCm39 Gencode M33 transcripts

#### Version Delimiters

The following delimiters are accepted for version specification: `/`, `--`, `=`, `@`

**Important**: Because `--` serves as a delimiter, it should not be used within overlay names themselves.

## 🐍 Leveraging Conda Resources

**CondaTainer** leverages the extensive `conda-forge` and `bioconda` ecosystems, which provide most bioinformatics software as conda packages. **CondaTainer** will automatically create module overlays for these packages.

For software unavailable through conda, custom build scripts can be created to download and install the software.

**Examples**: 10X [cellranger/9.0.1](https://github.com/Justype/condatainer/blob/main/build-scripts/cellranger/9.0.1) and Illumina [orad/2.7.0](https://github.com/Justype/condatainer/blob/main/build-scripts/orad/2.7.0)

## 🔄 Module Workflow

- [Manage Module Overlays](./module_overlays.md)
- [Bundle Overlays: Read-only project](./bundle_overlays.md)
- [Workspace Overlays: Writable project](./workspace_overlays.md)
