# Custom OS Overlays

OS Overlays can be created using Apptainer definition files or pulled from Docker Hub, allowing you to load os-level tools and libraries.

```bash
# Directly pull from a remote Docker image
condatainer create -p <prefix> -s docker://<docker_image>

# Create system wide overlay
condatainer create -n <name> -s docker://<docker_image>

# Build from a custom definition file
condatainer create -p <prefix> -f <path_to_def_file>
# -f cannot be used with -n
```

After creating the OS overlay, you may need to clean the Apptainer cache to save disk space:

```bash
apptainer cache clean
```

Examples:

- [Pulling PyTorch Docker Image](#example-pulling-pytorch-docker-image)
- [Pulling Posit R Docker Image](#example-pulling-posit-r-docker-image)
- [R Package Dependencies](#example-r-package-dependencies) (use with `env.img` for development)
- [Read-only R Package Environment](#example-read-only-r-package-environment) (for production)

## OS Overlay Compatibility

When loading multiple OS overlays together, all should share the same OS distro version. Mixing overlays built on different versions (e.g., Ubuntu 22.04 vs 24.04) can cause missing library errors or unexpected behavior.

To check an overlay's distro version, run:

```bash
condatainer info <overlay_file.sqf>
```

Example output:

```
File
  Name:          ubuntu22--igv.sqf
  Path:          path/images/ubuntu22--igv.sqf
  Size:          271.91 MB
  Type:          OS Overlay (Read-Only)
  Distro:        Ubuntu 22.04 (Jammy)
...
```

To get the base image distro version:

```bash
condatainer config get default_distro   # e.g. ubuntu24  →  24.04
# or
condatainer info base_image
```

The output is `ubuntu24`, which is different from the overlay's distro version `ubuntu22`. This means you cannot load this overlay together with other overlays built on `ubuntu24` base image.

## Change the Base Image

### Change the config (permanent)

If you want to permanently switch the default distro (e.g. to Ubuntu 22), run:

```bash
condatainer config set default_distro ubuntu22
```

This often happens when:

- The base image is end of life and you want to switch to a newer one.
- Most of your overlays target a different distro, so you want to switch the default to avoid compatibility issues.

### Use ENV variable (temporary)

Use `CONDATAINER_DEFAULT_DISTRO` to temporarily override the default distro for a single command without changing your global config:

```bash
CONDATAINER_DEFAULT_DISTRO=ubuntu22 condatainer exec -o myoverlay.sqf bash
```

### Use a different base image (temporary)

You can also directly specify a different base image that matches the overlay's distro:

```bash
condatainer exec -b ubuntu22/base_image.sqf -o myoverlay.sqf bash
```

You need to create it first. But with env, it will be auto-created if it does not exist.

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

It is Ubuntu 22.04, which differs from the default base image (Ubuntu 24.04 when `default_distro: ubuntu24`).

You can either:
- Use the PyTorch image as the base image
- Use the `ubuntu22--base_image.sqf` as the base image alongside the PyTorch overlay

Use the first approach:

```bash
condatainer exec -b pytorch.sqf bash
```

Second approach use Ubuntu 22 base image:

```bash
CONDATAINER_DEFAULT_DISTRO=ubuntu22 condatainer exec -o pytorch.sqf bash

# or explicitly specify the base image
condatainer exec -b ubuntu22/base_image.sqf -o pytorch.sqf bash
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

```{note}
For the `From` line, make sure to specify the same distro version as your base image (e.g., `ubuntu:24.04` if your `default_distro` is `ubuntu24`). Otherwise, you may encounter compatibility issues when loading the overlay.
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
