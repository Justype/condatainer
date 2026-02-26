# Run RStudio Server (Conda) on HPC

Simply run the following commands to set up and run RStudio Server with Conda-managed R on your HPC system:

```bash
# download/update helper scripts
condatainer helper -u
# create a 30G overlay to store conda env
condatainer o -s 30g
# install R in the conda env
condatainer e -- mm-install r-base=4.4 r-tidyverse -y
# pin R version to avoid accidental updates
condatainer e -- mm-pin r-base
# start RStudio Server
condatainer helper rstudio-server-conda
```

This tutorial uses **Conda-managed R** (via `mm-install r-base`). If you you need to build R libraries from source (e.g., GitHub packages), see [RStudio Server](./rstudio-server_on_HPC.md) instead.

## When to Use This Variant

Use `rstudio-server-conda` when:

- You want all R packages managed through Conda (CRAN packages as `r-*`, Bioconductor as `bioconductor-*`)
- You need tight integration between R and Python packages in the same environment
- You want to easily export and reproduce your entire environment with `mm-export`

Use the regular `rstudio-server` when:

- You prefer normal R package management via `install.packages()`, `BiocManager::install()`, etc.
- You need to build R packages from source (e.g., GitHub packages)

## Checklist

- [A supported scheduler is available on your HPC system](#scheduler)
- [Have SSH port forwarding set up to the login node](#ssh-port-forwarding)
- [Have CondaTainer installed](#install-condatainer)
- [Have R installed in a writable overlay image](#create-r-writable-overlay)
- [Check the Script Parameters](#rstudio-server-conda-helper-script)

Then you can run:

```bash
condatainer helper rstudio-server-conda
# Default
#   -p port_number=auto-selected if omitted
#   -c num_cpus=4
#   -m memory=32G
#   -t time_limit=12:00:00
```

You can set [Config](./helpers_on_HPC.md#configuration-and-defaults) and helpers also support [Reusing previous settings](./helpers_on_HPC.md#reuse-mode).

## Scheduler

Run the following command to check if a supported scheduler is available:

```bash
condatainer scheduler
```

```{note}
Slurm is fully tested. Other schedulers are experimentally supported. Please report any issues you encounter.
```

## SSH Port Forwarding

```
Your machine -----> HPC Login Node -----> HPC Compute Node
        (port forwarding)     (port forwarding)
```

`rstudio-server` is a web-based application, so you need to access it through a port on the HPC system.

But you cannot directly access compute nodes from your local machine. To accomplish this, you need to set up SSH port forwarding from your local machine to the HPC login node.

```bash
ssh -L <local_port>:localhost:<remote_port> your_username@hpc_address
```

After running this command, when your local machine accesses `localhost:<local_port>`, it will be forwarded to `localhost:<remote_port>` on the HPC login node.

```{warning}
HPC is a shared system! Please do not use a common port like `8787`.

Instead, choose an uncommon and unique port number like `13182`.
```

You can also edit the SSH config file (`~/.ssh/config`) to add a section for your HPC system:

```
Host hpc
    HostName hpc_address
    User your_username
    LocalForward <local_port> localhost:<remote_port>
    # Please change <local_port> to a unique port number
```

Then you can simply run:

```bash
ssh hpc
```

## Install CondaTainer

Run the following command to install CondaTainer if it is not installed:

```bash
curl -fsSL https://get-condatainer.justype.net | bash
```

## Create R Writable Overlay

Creating an ext3 overlay image with Conda-managed R:

```bash
# create a 30G ext3 image named `env.img` in the current working directory
condatainer o -s 30g

# Launch a shell within the overlay
condatainer e
```

Inside the overlay, install R using Conda:

```bash
# install R and commonly used packages
mm-install r-base=4.4 r-tidyverse

# pin the R version to avoid accidental updates
mm-pin r-base
```

If you need `reticulate`, you can also install conda:

```bash
# You can change python version and add packages as needed
mm-install python=3.11 conda
```

See [Launch a Shell within the Workspace Overlay](../user_guide/workspace_overlays.md#launch-a-shell-within-the-workspace-overlay) for more details.

```{note}
Always use Conda to install R packages. If you want to use `install.packages()` in R, use [rstudio-server](./rstudio-server_on_HPC.md) variant instead.
```

## RStudio Server Conda Helper Script

It will do the following steps for you:

1. Check if `rstudio-server-conda` is running on any compute node.
2. If yes, establish SSH port forwarding to that node.
3. If not,
   1. Check port and overlay integrity.
   2. Verify R is installed in the overlay.
   3. Submit the scheduler job to start `rstudio-server-conda`.
   4. When the job starts, record and set up SSH port forwarding.

```
Usage: rstudio-server-conda [options]

Options:
  -c <number>     Number of CPUs to allocate (default: 4)
  -m <memory>     Amount of memory to allocate (default: 32G)
  -t <time>       Time limit for the job (default: 12:00:00)
  -g <gpu>        GPU resources (e.g., 1, a100:2). Empty means no GPU.
  -v              View Mode NCPUS:1 MEM:8G TIME:02:00:00

  -p <port>       Port for rstudio-server (default: randomly picked). Valid range: 1024-65535.
  -e <overlay>    Environment overlay image file (default: env.img)
  -o <overlay>    Additional overlay files (can have multiple -o options)

  config          Show config file path and contents
```

## Configuration

The script saves its defaults to `$XDG_CONFIG_HOME/condatainer/helper/rstudio-server-conda` (defaults to `~/.config/condatainer/helper/rstudio-server-conda`) on first run. Subsequent runs load from this file.

```bash
# View current config
condatainer helper rstudio-server-conda config

# Reset to defaults (delete config file, next run recreates it)
rm ~/.config/condatainer/helper/rstudio-server-conda
```

## Running RStudio Server

Let's set up and run `rstudio-server-conda` on HPC:

```bash
# Download the helper scripts
condatainer helper -u
```

Then you can run the script:

```bash
condatainer helper rstudio-server-conda
```

After running the script, you will see output like this:

```
Job 1299123 is now running on node cm23.
rstudio-server-conda at http://localhost:<port>
You can run the following command in R to open the project directly:
  rstudioapi::openProject("/scratch/your_username/current_working_directory")
If you want to stop it, run: scancel 1299123
```

Now you can go to your local web browser and access the link. Then run the provided R command to open the project.

Don't forget to set up SSH port forwarding from your local machine to the HPC login node before accessing the link.

## Install Packages with Conda

Since R is managed by Conda, you can install R packages through Conda channels:

- CRAN packages: `r-<lowercase_name>` (e.g., `r-ggplot2`, `r-dplyr`)
- Bioconductor packages: `bioconductor-<lowercase_name>` (e.g., `bioconductor-deseq2`)

```bash
# Inside the overlay shell (condatainer e)

# Search for available R packages
mm-search r-presto

# Install packages
mm-install r-ggplot2 bioconductor-deseq2

# Pin packages to prevent accidental updates
mm-pin r-base
```

### Advantages of Conda-managed R Packages

1. **No compilation required** - Binary packages from conda-forge
2. **Automatic dependency management** - System libraries handled automatically
3. **Easy environment export** - Reproduce your environment exactly:

```bash
mm-export --no-builds > conda-env.yml
```

4. **Python integration** - Install Python packages alongside R:

```bash
mm-install python=3.11 numpy pandas
```

```{note}
Always use conda to install R packages within the overlay. Do not use `install.packages()` in R, as it may lead to inconsistencies.

If required packages are not available in conda-forge or bioconda, consider using the [Rstudio Server on HPC](./rstudio-server_on_HPC.md) variant instead.
```

## Run without RStudio Server

You can run R script directly without RStudio Server:

```bash
condatainer exec -o env.img Rscript your_script.R
```

## Using GPU

For GPU-accelerated R packages, specify GPU resources:

```bash
condatainer helper rstudio-server-conda -g 1        # Request 1 GPU
condatainer helper rstudio-server-conda -g a100:2   # Request 2 A100 GPUs
```

You can use the following command to check available GPU types on your HPC system:

```
condatainer scheduler
```
