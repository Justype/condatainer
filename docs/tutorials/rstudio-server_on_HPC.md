# Run rstudio-server on HPC

Simply run the following commands to set up and run RStudio Server on your HPC system:

```bash
condatainer helper -u
condatainer helper rstudio-server
```

## Checklist

- [SLURM job scheduler is available on your HPC system](#slurm-job-scheduler)
- [Have SSH port forwarding set up to the login node](#ssh-port-forwarding)
- [Have CondaTainer installed](#install-condatainer)
- [Have required overlay images created](#install-required-overlays)
- [Have R installed in a writable overlay image](#create-r-writable-overlay)
- [Check the Script Parameters](#rstudio-server-helper-script)

Then you can run:

```bash
condatainer helper rstudio-server
# Default
# -p  port_number=auto-selected if omitted
# -c  num_cpus=4
# -m  memory=32G
# -t  time_limit=12:00:00
```

You can create alias in your shell config file (`~/.bashrc` or `~/.zshrc`):

```bash
# Change 13182 to your preferred port number
alias rstudio-server-start='condatainer helper rstudio-server -p 13182'
```

See [Install Packages in CondaTainer RStudio Server](#install-packages-in-condatainer-rstudio-server) for more details on installing R packages.

## SLURM Job Scheduler

Run the following command to check if SLURM is available on your HPC system:

```bash
sbatch --version
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

Run the following command to check if CondaTainer is installed:

```
condatainer --version
```

Run the following command to install CondaTainer if it is not installed:

```bash
curl -fsSL https://get-condatainer.justype.net/ | bash
```

Download the helper scripts:

```bash
condatainer helper --update
```

## Install Required Overlays

Creating required overlay images:

```bash
condatainer install rstudio-server build-essential
```

## Create R Writable Overlay

Creating an ext3 overlay image with R

```bash
# create a 30G ext3 image named `env.img` in the current working directory
condatainer overlay -s 30720

# go into the overlay
condatainer e
```

In side the overlay, you can use `mm-*` commands

```bash
# install R and required packages
mm-install r-base=4.4 r-tidyverse
```

See [Launch a Shell within the Workspace Overlay](../user_guide/workspace_overlays.md#launch-a-shell-within-the-workspace-overlay) for more details.

## RStudio Server Helper Script

It will do the following steps for you:

1. Check if `rstudio-server` is running on any compute node.
2. If yes, establish SSH port forwarding to that node.
3. If not,
   1. Check port and overlay integrity.
   2. Submit the SLURM job to start `rstudio-server`.
   3. When the job starts, record and set up SSH port forwarding.

```
Usage: rstudio-server [options]

Options:
  -c <number>     Number of CPUs to allocate (default: 4)
  -m <memory>     Amount of memory to allocate (default: 32G)
  -t <time>       Time limit for the job (default: 12:00:00)
  -p <port>       Port for rstudio-server (if omitted, an available port will be chosen)
  -o <overlay>    Additional overlay files (can have multiple -o options)
  -h              Show this help message
```

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
If you want to stop it, run: scancel 1299123
You can run the following command in R to open the project directly:
  rstudioapi::openProject("/scratch/your_username/current_working_directory")
```

Now you can go to your local web browser and access the link. Then run the provided R command to open the project.

Don't forget to set up SSH port forwarding from your local machine to the HPC login node before accessing the link.

## Install Packages in CondaTainer RStudio Server

There are two ways to install R packages in CondaTainer RStudio Server:

1. [Use Conda to manage R packages](#use-conda-to-manage-r-packages)
2. [The normal way](#install-r-packages-without-conda): `install.packages()`, `pak::pak()`, etc. (recommended)

```{note}
Both methods will install packages into the writable overlay image. `/ext3/env/lib/R/library/`
```

### Use Conda to Manage R Packages

If you only use R packages from CRAN/Bioconductor. Those packages are commonly available from conda-forge/bioconda channel.

- CRAN packages: `r-<lowercase_name>` (e.g., `r-ggplot2`)
- Bioconductor packages: `bioconductor-<lowercase_name>` (e.g., `bioconductor-deseq2`)

e.g. (In the overlay shell)

```bash
# find available R packages in conda-forge/bioconda
# Or go to https://anaconda.org
mm-search r-presto

# Install packages
mm-install r-base=4.4 r-ggplot2 bioconductor-deseq2

# Pin the R version to avoid accidental updates
mm-pin r-base
```

If you only install R packages from conda-forge/bioconda, you don't need to worry about Rprofile or system dependencies. You can easily export the conda environment later for reproducibility:

```bash
mm-export --no-builds > conda-env.yml
```

And you don't need to worry about system dependencies, as conda will handle them for you.

```bash
condatainer exec -o env.img Rscript your_script.R
```

### Install R Packages without Conda

The `build-essential` overlay installed earlier contains common system libraries needed for building R packages. So you can directly install R packages from CRAN/Bioconductor/GitHub/source.

```R
install.packages("tidyverse")
BiocManager::install("DESeq2")
pak::pak("user/repo") # Or remotes::install_github("user/repo")
```

#### Rprofile Setup

If `*.Rproj` file does not exist in the current working directory, **CondaTainer** will create:

- A `*.Rproj` file with default settings.
- A `.Rprofile` file to set up Posit Public Package Manager (P3M) CRAN and Bioconductor repositories.

See [Rprofile](https://github.com/Justype/condatainer/blob/main/helpers/slurm/.Rprofile) for more details. At the end, it will call `set_repository_options()` for you. So you can directly download the binary packages from P3M.

If `*.Rproj` file already exists, `.Rprofile` will not be created or modified. You need to set up the repositories yourself if needed.

If you only use CRAN packages, you can change the function in the `.Rprofile` to:

```R
set_repository_options(repo = "cran", latest_cran = TRUE)
```

But if you plan to use Bioconductor packages, please keep the default settings to ensure compatibility.

#### Install R Packages from Source

If a package is only available from source (e.g., GitHub), you need `remotes`, `pak` or other package managers.

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

If you only use R packages from Conda, you can run R directly without RStudio Server:

```bash
condatainer exec -o env.img Rscript your_script.R
```

If you have R packages built from GitHub or source, you need to load additional overlay too:

```bash
# If no additional overlay is needed, just run:
condatainer exec -o build-essential -o env.img Rscript your_script.R
```

```bash
# If you created additional overlay `r-deps.sqf`, run:
condatainer exec -o build-essential -o r-deps.sqf -o env.img Rscript your_script.R
```
