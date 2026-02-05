# Run RStudio Server on HPC

Simply run the following commands to set up and run RStudio Server on your HPC system:

```bash
condatainer helper -u
condatainer helper rstudio-server
```

This tutorial uses **Posit R image overlays** for R installation. If you prefer using Conda-managed R (via `mm-install r-base`), see [RStudio Server (Conda)](./rstudio-server-conda_on_HPC.md) instead.

## Checklist

- [A supported scheduler is available on your HPC system](#scheduler--workload-manager)
- [Have SSH port forwarding set up to the login node](#ssh-port-forwarding)
- [Have CondaTainer installed](#install-condatainer)
- [Have a writable overlay image](#create-writable-overlay)
- [Check the Script Parameters](#rstudio-server-helper-script)

Then you can run:

```bash
condatainer helper rstudio-server
# Default
#   -r r_version=latest available (e.g., 4.5.2)
#   -p port_number=auto-selected if omitted
#   -c num_cpus=4
#   -m memory=32G
#   -t time_limit=12:00:00
```

You can create an alias in your shell config file (`~/.bashrc` or `~/.zshrc`):

```bash
# Change 13182 to your preferred port number
alias rstudio-start='condatainer helper rstudio-server -p 13182 -r 4.4'
```

## Scheduler / Workload Manager

Run the following command to check if a supported scheduler is available:

```bash
condatainer scheduler
```

Currently only **SLURM** is supported. PBS support is planned for future releases.

## SSH Port Forwarding

```
Your machine -----> HPC Login Node -----> HPC Compute Node
        (port forwarding)     (port forwarding)
```

RStudio Server is a web-based application, so you need to access it through a port on the HPC system.

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
curl -fsSL https://get-condatainer.justype.net/ | bash
```

## Create Writable Overlay

Creating an ext3 overlay image for your packages and project files:

```bash
# create a 30G ext3 image named `env.img` in the current working directory
condatainer o -s 30g
```

If you need `reticulate`, you need to set up a conda environment inside the overlay:

```bash
# Launch a shell within the overlay
condatainer e
```

```bash
# You can change python version and add packages as needed
mm-install python=3.11 conda
```

See [Launch a Shell within the Workspace Overlay](../user_guide/workspace_overlays.md#launch-a-shell-within-the-workspace-overlay) for more details.

## RStudio Server Helper Script

The script will do the following steps for you:

1. Check if `rstudio-server` is running on any compute node.
2. If yes, establish SSH port forwarding to that node.
3. If not,
   1. Check port and overlay integrity.
   2. Install R version overlay if needed.
   3. Submit the SLURM job to start `rstudio-server`.
   4. When the job starts, record and set up SSH port forwarding.

```
Usage: rstudio-server [options]

Options:
  -r <r_version>  R version to use (e.g., 4.4). If omitted, uses latest available.
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

### Available R Versions

The helper script uses [Posit R docker images](https://hub.docker.com/r/posit/r-base).

You can specify partial versions (e.g., `-r 4.4` will use the latest 4.4.x).

## Configuration

The script saves its defaults to `~/.config/condatainer/helper/defaults/rstudio-server` on first run. Subsequent runs load from this file.

```bash
# View current config
condatainer helper rstudio-server config

# Reset to defaults (delete config file, next run recreates it)
rm ~/.config/condatainer/helper/defaults/rstudio-server
```

## Running RStudio Server

Let's set up and run `rstudio-server` on HPC:

```bash
# Download the helper scripts
condatainer helper -u
```

Then you can run the script:

```bash
condatainer helper rstudio-server
```

After running the script, you will see output like this:

```
Waiting for job 1299123 to start running. Current status: PENDING.
Job 1299123 is now running on node cm23.
rstudio-server at http://localhost:<port>
You can run the following command in R to open the project directly:
  rstudioapi::openProject("/scratch/your_username/current_working_directory")
If you want to stop it, run: scancel 1299123
```

Now you can go to your local web browser and access the link. Then run the provided R command to open the project.

Don't forget to set up SSH port forwarding from your local machine to the HPC login node before accessing the link.

## Install Packages in CondaTainer RStudio Server

The `build-essential` overlay (auto-installed) contains common system libraries needed for building R packages. So you can directly install R packages from CRAN/Bioconductor/GitHub/source.

```R
install.packages("tidyverse")
BiocManager::install("DESeq2")
pak::pak("user/repo") # Or remotes::install_github("user/repo")
```

### Rprofile Setup

If `*.Rproj` file does not exist in the current working directory, **CondaTainer** will create:

- A `*.Rproj` file with default settings.
- A `.Rprofile` file to set up Posit Public Package Manager (P3M) CRAN and Bioconductor repositories.

See [Rprofile](https://github.com/Justype/condatainer/blob/main/helpers/slurm/.Rprofile) for more details. At the end, it will call `set_repository_options()` for you. So you can directly download binary packages from P3M.

If `*.Rproj` file already exists, `.Rprofile` will not be created or modified. You need to set up the repositories yourself if needed.

So if you do not want to use Posit repositories, you can remove or edit the `.Rprofile` file accordingly.

If you only use CRAN packages, you can change the function in the `.Rprofile` to:

```R
set_repository_options(repo = "cran", latest_cran = TRUE)
```

But if you plan to use Bioconductor packages, please keep the default settings to ensure compatibility.

### Install R Packages from Source

If a package is only available from source (e.g., GitHub), you need `remotes`, `pak` or other package managers. (Common system dependencies are included in `build-essential` overlay.)

I recommend using `pak`, it can tell you which system libraries are required.

```R
pak::pkg_sysreqs("user/repo@commit_hash") # or @tag
```

If system libraries are missing, you can create your `additional-deps` overlay with the required system libraries. See [Custom OS Overlays](../advanced_usage/custom_os.md) for more details.

After getting the `additional-deps.sqf` overlay, run `rstudio-server` with the `-o` option:

```bash
condatainer helper rstudio-server -o additional-deps.sqf
```

If you want to share your overlay with others, you should also provide the `def` file used to create it.

## Run without RStudio Server

You can run R scripts directly without RStudio Server:

```bash
# Basic execution
condatainer exec -o r4.4.3 -o build-essential -o env.img Rscript your_script.R
```

```bash
# If you created additional overlay `r-deps.sqf`:
condatainer exec -o r4.4.3 -o build-essential -o r-deps.sqf -o env.img Rscript your_script.R
```

## Using GPU

For GPU-accelerated R packages, specify GPU resources:

```bash
condatainer helper rstudio-server -g 1        # Request 1 GPU
condatainer helper rstudio-server -g a100:2   # Request 2 A100 GPUs
```

You can use the following command to check available GPU types on your HPC system:

```
condatainer scheduler
```

Example output:

```
Scheduler Information:
  Type:      SLURM
  Binary:    /opt/software/slurm/current/bin/sbatch
  Version:   24.11.7
  Status:    Available

The scheduler is available and ready for job submission.

Max Resource Limits:
  Max CPUs:   44
  Max Memory: 184000 MB (180 GB)
  Max Time:   7d

Available GPUs:
  h100: 2144 total, 416 available
  mi300a: 4 total, 0 available
  nvidia_h100_80gb_hbm3_1g.10gb: 2264 total, 728 available
  nvidia_h100_80gb_hbm3_2g.20gb: 768 total, 0 available
  nvidia_h100_80gb_hbm3_3g.40gb: 768 total, 0 available
  t4: 28 total, 0 available
```
