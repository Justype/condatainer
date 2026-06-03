# RStudio Server

Two variants are available depending on how R is installed:

| Variant | R Installation | Best For |
|---------|----------------|----------|
| `rstudio-server` | Posit R image overlays | Normal `install.packages()`, GitHub, Bioconductor from source |
| `rstudio-server-conda` | Conda (`mm-install r-base`) | All packages via conda-forge/bioconda, reproducible env export |

```bash
condatainer helper -u
condatainer helper rstudio-server -w         # Posit R
condatainer helper rstudio-server-conda -w   # Conda R
```

See [Helper Scripts](./helpers.md) for SSH port forwarding setup, resource flags, configuration, and reuse mode.

## rstudio-server (Posit R)

Uses [Posit R image overlays](https://hub.docker.com/r/posit/r-base) — a separate `r<version>` overlay is auto-installed alongside `rstudio-server` and `build-essential`.

### Writable Overlay

```bash
# Create a 30G writable overlay in the current directory
condatainer o -s 30g

# With Python (needed for reticulate):
condatainer o -s 30g -- python=3.11 conda
```

### Flags

```bash
condatainer helper rstudio-server --help
```

| Flag | Description | Default |
|------|-------------|---------|
| `--rversion,-r` | R version (e.g. `4.5.2`, `4.4`) | latest available |
| `--password,-p` | RStudio password (empty = no auth) | (empty) |

Partial versions are accepted: `-r 4.4` uses the latest 4.4.x available. See [Helper Scripts](./helpers.md#flags) for resource and overlay flags.

### Installing R Packages

The `build-essential` overlay provides common system libraries. Install packages normally:

```R
install.packages("tidyverse")
BiocManager::install("DESeq2")
pak::pak("user/repo")
```

#### Rprofile Setup

When no `*.Rproj` file exists in the working directory, CondaTainer creates:

- A `*.Rproj` file with default settings
- A `.Rprofile` configured for Posit Public Package Manager (P3M) binary packages

If `*.Rproj` already exists, `.Rprofile` is not modified. To use only CRAN:

```R
.set_repository_options(repo = "cran", latest_cran = TRUE)
```

#### Missing System Libraries

Create a custom overlay for extra system libraries — see [Custom OS Overlays](../advanced_usage/custom_os.md):

```bash
condatainer helper rstudio-server -o additional-deps.sqf
```

### Run R Without RStudio

```bash
condatainer exec -o r4.4.3 -o build-essential -o env.img Rscript script.R
condatainer exec -o r4.4.3 -o build-essential -o r-deps.sqf -o env.img Rscript script.R
```

## rstudio-server-conda (Conda R)

R is installed in the writable overlay via `mm-install r-base`. No separate R overlay is needed.

```{note}
Always use `mm-install` (Conda) to install R packages in this variant. Do not use `install.packages()` — it may conflict with the Conda-managed environment.
```

### Writable Overlay with R

```bash
# Create a 30G overlay with R pre-installed
condatainer o -s 30g -- r-base=4.4 r-tidyverse

# Enter the overlay to install more packages
condatainer e
```

Inside the overlay:

```bash
mm-pin r-base   # pin R version to prevent accidental updates
mm-install r-seurat r-patchwork bioconductor-clusterprofiler
mm-install python=3.11 conda   # add Python for reticulate
```

### Package Management

Conda naming conventions:
- CRAN: `r-<name>` (e.g. `r-ggplot2`, `r-dplyr`)
- Bioconductor: `bioconductor-<name>` (e.g. `bioconductor-deseq2`)

```bash
mm-search r-presto
mm-install r-ggplot2 bioconductor-deseq2
mm-export --no-builds > conda-env.yml   # export full environment
```

### Run R Without RStudio

```bash
condatainer exec -o env.img Rscript script.R
```

Once running, open the project in R:

```R
rstudioapi::openProject("/scratch/user/my-project")
```

## Common Issues

### UID conflict (UID = 1000)

The Posit R container image reserves UID 1000 for its internal user. If your account UID is 1000, `rstudio-server` will not start — use `rstudio-server-conda` instead.
