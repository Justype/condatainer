# Custom OS Overlays

OS Overlays can be created using Apptainer definition files or pulled from Docker Hub, allowing you to load os-level tools and libraries.

```bash
# Directly pull from a remote Docker image
condatainer create -p <prefix> -s docker://<docker_image>

# Build from a custom definition file
condatainer create -p <prefix> -f <path_to_def_file>
```

Examples:

- [Pulling PyTorch Docker Image](#example-pulling-pytorch-docker-image)
- [R Package Dependencies](#example-r-package-dependencies) (use with `env.img` for development)
- [Read-only R Package Environment](#example-read-only-r-package-environment) (for production)

## Change the Base Image

If an OS overlay has a different distro version than the base image used by other overlays, you may run into compatibility issues when loading multiple overlays together.

To avoid this, you can:

```bash
# Use a new base image matching the OS overlay distro
condatainer exec -b <custom_base_image.sqf> -o <os_overlay.sqf> bash

# or directly use the OS overlay as the base
condatainer exec -b <os_overlay.sqf> bash
```

## What is Included in a Base Image?

A base image typically includes:

- `micromamba` and `mm-*` scripts for managing conda environments.
- `apptainer` for nested container support.

For example, you have already launched the `code-server`. The base image allows you to:

1. Directly manage conda env under the `/ext3/env` directory. (if in writable mode)
2. Launch nested containers (like read module sqf like `grch38/gtf-gencode/47`)

If you are not using one of the features above, you can directly use another os overlay as the base image.

## Example: Pulling PyTorch Docker Image

For example, you might want to directly pull the PyTorch Docker image (based on `ubuntu:22.04`):

When you use that overlay directly, you may encounter missing library errors in Python.

```bash
condatainer create -p pytorch -s docker://pytorch/pytorch:2.9.1-cuda13.0-cudnn9-devel

# This may produce an error when you run:
condatainer exec -o pytorch.sqf bash
```

In the inner shell, you can see the OS version with:

```bash
cat /etc/os-release
```

It is Ubuntu 22.04, which differs from the default base image (`ubuntu:24.04`).

You can either:
- Use the PyTorch image as the base image
- Build an Ubuntu 22.04 base image

Use the first approach:

```bash
condatainer exec -b pytorch.sqf bash
```

Or create an Ubuntu 22.04 base image:

```bash
# Use available base image definition file
condatainer create ubuntu22/base_image

# Then load both the base image and the PyTorch overlay:
condatainer exec -b ubuntu22/base_image -o pytorch.sqf bash
```

```{note}
Available base image definition files:

- [ubuntu20/base_image](https://github.com/Justype/condatainer/tree/main/build-scripts/ubuntu20/base_image.def)
- [ubuntu22/base_image](https://github.com/Justype/condatainer/tree/main/build-scripts/ubuntu22/base_image.def)
- [ubuntu24/base_image](https://github.com/Justype/condatainer/tree/main/build-scripts/ubuntu24/base_image.def)
```

## Example: R Package Dependencies

Let's say when building an R package from source, you find that some libraries are missing, e.g., `libxml2-dev`, `libcurl4-openssl-dev`, and `libssl-dev`.

```{tip}
You can use `pak::pkg_sysreqs()` in R to check which libraries are required for specific packages.
```

### 1. Create a Custom Definition File

You can create a custom Apptainer definition file to include these dependencies.

`r-deps.def`

```
Bootstrap: docker
From: ubuntu:24.04
# Make sure distro version matches the base image

%post
    # Avoid interactive prompts
    export DEBIAN_FRONTEND=noninteractive

    # Put your custom installation commands here
    apt -y update && apt -y install \
        libxml2-dev libcurl4-openssl-dev libssl-dev

    apt clean && rm -rf /var/lib/apt/lists/*
```

### 2. Build the Overlay

Use the `condatainer create` command to build the overlay from your custom definition file.

```bash
condatainer create -p r-deps -f r-deps.def
```

This will generate an overlay file named `r-deps.sqf` in the current directory.

### 3. Use the Custom Overlay

You can now load this custom overlay when running analyses that require these additional system libraries.

```bash
condatainer exec -o r-deps.sqf your_command_here
```

In `rstudio-server-helper`:

```bash
condatainer helper rstudio-server -o r-deps.sqf
```

### 4. Share the File

Share the `r-deps.def` file with your team or include it in your project repository. This allows others to build the same overlay and ensures consistency across different environments.

Or you can directly share the `r-deps.sqf` overlay file for immediate use.

## Example: Read-only R Package Environment

**CondaTainer** is good for managing writable environments during development. If you have already developed a pipeline and fixed package versions, you can directly build a read-only Apptainer (Singularity) image with all dependencies included.

```{tip}
- Use Posit Package Manager CRAN mirror for precompiled R packages to speed up installation.
- Use `pak::pak()` to directly run `apt` or other system package managers for you.
```

Example: `sc-run-standalone.def`

```
Bootstrap: docker
From: posit/r-base:4.4.3-noble-amd64

%post
    # Avoid interactive prompts
    export DEBIAN_FRONTEND=noninteractive

    R -e '
options(repos = c(CRAN = sprintf("https://packagemanager.posit.co/cran/latest/bin/linux/noble-%s/%s", R.version["arch"], substr(getRversion(), 1, 3))))
install.packages("pak")
pak::pak(c(
    "Seurat@5.3.0",
    "chris-mcginnis-ucsf/DoubletFinder@1b244d8",
    "navinlabcode/copykat@b7a4763"
))
'
```

```bash
apptainer build sc-run-standalone.sif sc-run-standalone.def
```

Then run:

```bash
apptainer exec sc-run-standalone.sif \
    Rscript --vanilla -e "library(DoubletFinder); library(copykat); sessionInfo()"
```

```bash
# Or you can let condatainer set --nv --bind --env for you
condatainer exec -b sc-run-standalone.sif \
    Rscript --vanilla -e "library(DoubletFinder); library(copykat); sessionInfo()"
```

```
R version 4.4.3 (2025-02-28)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 24.04.3 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0

locale:
 [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
 [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
 [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
[10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   

time zone: UTC
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] copykat_1.1.0       DoubletFinder_2.0.6
```
