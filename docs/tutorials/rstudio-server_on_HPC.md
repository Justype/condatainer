# Run rstudio-server on HPC

## Checklist

- [SLURM job scheduler is available on your HPC system](#slurm-job-scheduler)
- [Have SSH port forwarding set up to the login node](#ssh-port-forwarding)
- [Have CondaTainer installed](#install-condatainer)
- [Have required overlay images created](#install-required-overlays)
- [Have R installed in a writable overlay image](#create-r-writable-overlay)
- [Have the `rstudio-server-helper` in your PATH](#rstudio-server-helper-script)

Then you can run:

```bash
rstudio-server-helper \
  -p <port_number> \
  -e <overlay_image> \
  -c <num_cpus> \
  -m <memory> \
  -t <time_limit>

# Default 
#   port_number=8787
#   overlay_image=env.img
#   num_cpus=4
#   memory=32G
#   time_limit=12:00:00
```

Please change the port number to a unique one to avoid conflicts with other users.

Always [Use Conda to Manage R Packages](#use-conda-to-manage-r-packages). For packages only available from source, see [Install R Packages from Source](#install-r-packages-from-source) section.

If you have preferred settings, you can modify the script directly. See [Change the default setting](#change-the-default-setting) section below.

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
ssh -L 8787:localhost:8787 your_username@hpc_address
```

After running this command, when your local machine accesses `localhost:8787`, it will be forwarded to `localhost:8787` on the HPC login node.

```{warning}
HPC is a shared system! Please do not use a common port like `8787`.

Instead, choose an uncommon and unique port number like `13182`.
```

You can also edit the SSH config file (`~/.ssh/config`) to add a section for your HPC system:

```
Host hpc
    HostName hpc_address
    User your_username
    LocalForward 8787 localhost:8787
    # Please change 8787 to a unique port number
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
curl -sL https://raw.githubusercontent.com/Justype/condatainer/main/assets/install.sh | bash
```

## Install Required Overlays

Creating required overlay images:

```bash
condatainer install rstudio-server
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

See [Launch a Shell with the Project Environment](../user_guide/condatainer_project_level_writable.md#launch-a-shell-with-the-project-environment) for more details.

## RStudio Server Helper Script

You should always run `rstudio-server` on compute nodes rather than login nodes.

You can use this script: [rstudio-server-helper](https://github.com/Justype/condatainer/blob/main/helpers/rstudio-server-helper)

It will do the following steps for you:

1. Check if `rstudio-server` is running on any compute node.
2. If yes, establish SSH port forwarding to that node.
3. If not,
   1. Check if the port is available on the login node.
   2. If so, check the overlay integrity and make sure R is installed.
   3. If so, use `sbatch` to submit a job which starts `rstudio-server`.
   4. Wait for the job to start and get the compute node name.
   5. Record the job ID, node name, port number, and working directory.
   6. Set up SSH port forwarding from login node to compute node.

```
Usage: rstudio-server-helper [options]

Options:
  -c <number>     Number of CPUs to allocate (default: 4)
  -m <memory>     Amount of memory to allocate (default: 32G)
  -t <time>       Time limit for the job (default: 12:00:00)
  -p <port>       Port for rstudio-server (default: 8787)
  -e <overlay>    Environment overlay image file (default: env.img)
  -o <overlay>    Additional overlay files (can have multiple -o options)
  -y              Accept all warnings and proceed
  -h              Show this help message
```

Let's set up and run `rstudio-server` on HPC:

```bash
# Download the helper script
# Please make sure $HOME/bin is in your PATH
mkdir -p $HOME/bin
wget https://raw.githubusercontent.com/Justype/condatainer/main/helpers/rstudio-server-helper -O $HOME/bin/rstudio-server-helper
chmod +x $HOME/bin/rstudio-server-helper
```

Then you can run the script: 

```bash
rstudio-server-helper
```

After running the script, you will see output like this:

```
Waiting for job 1299123 to start running. Current status: PENDING.
Job 1299123 is now running on node cm23.
rstudio-server at http://localhost:8787
If you want to stop it, run: scancel 1299123
You can run the following command in R to open the project directly:
  rstudioapi::openProject("/scratch/your_username/current_working_directory")
```

Then you can go to your local web browser and access [http://localhost:8787](http://localhost:8787).

Run the provided R command in the R console to open the project directly.

## Use Conda to Manage R Packages

It is highly recommended to use Conda (via `mm-install`) to install R packages whenever possible. They are pre-compiled.

- CRAN packages: `mm-install r-<package_name>`
- Bioconductor packages: `mm-install bioconductor-<package_name>`

e.g.

```bash
# in the overlay using condatainer e

# find available R packages in conda-forge/bioconda
# Or go to https://anaconda.org
mm-search r-presto

# Install packages
mm-install r-base=4.4 r-ggplot2 bioconductor-deseq2

# Pin the R version to avoid accidental updates
mm-pin r-base
```

## Install R Packages from Source

```{note}
Before installing R packages from source, please make sure the package is not available from conda-forge or bioconda.
```

If a package is only available from source (e.g., GitHub), you need `remotes`, `pak` or other package managers.

I recommend using `pak`, it can tell which system libraries are required.

```R
mm-install r-pak
```

Then in R:

```R
pak::pkg_sysreqs("user/repo@commit_hash") # or @tag
```

If system libraries are missing, you can create your `additional-deps` overlay with the required system libraries. (ignore pandoc missing warning)

see [Custom Apptainer Definition Files](../advanced_usage/condatainer_custom_def.md) for more details.

After getting the `additional-deps.sqf` overlay, run `rstudio-server-helper` with the `-o` option:

```bash
rstudio-server-helper -o additional-deps.sqf
```

If you want to share your overlay with others, you should also provide the `def` file used to create it.

## Run without RStudio Server

If you only use R packages from Conda, you can run R directly without RStudio Server:

```bash
condatainer exec -o env.img Rscript your_script.R
```

If you have R packages built from GitHub or source, you need to load your additional overlay too:

```bash
condatainer exec -o additional-deps.sqf -o env.img Rscript your_script.R
```

## Change the default setting

You can modify the default settings in the script, such as CPU, memory, time limit, overlay image path, port number.

They are at the beginning of the script:

```bash
#!/bin/bash

NCPUS=4
MEM=32G
TIME=12:00:00
PORT=8787
```

You can use editors like `vim` or `nano` to change these values.

or use `sed`

```bash
sed -i 's/NCPUS=4/NCPUS=8/' $HOME/bin/rstudio-server-helper
sed -i 's/MEM=32G/MEM=64G/' $HOME/bin/rstudio-server-helper
sed -i 's/TIME=12:00:00/TIME=24:00:00/' $HOME/bin/rstudio-server-helper
sed -i 's/PORT=8787/PORT=13182/' $HOME/bin/rstudio-server-helper
```
