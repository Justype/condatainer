# How to run rstudio-server on HPC with CondaTainer

Before getting started, ensure you have [CondaTainer installed](../../README.md).

## Checklist

- [ ] [SLURM job scheduler is available on your HPC system](#slurm-job-scheduler)
- [ ] [Have SSH port forwarding set up to the login node](#ssh-port-forwarding)
- [ ] [Have CondaTainer installed](#install-condatainer)
- [ ] [Have required overlay images created](#install-required-overlays)
- [ ] [Have R installed in a writable overlay image](#create-r-writable-overlay)
- [ ] [Have the `rstudio-server-helper` in your PATH](#rstudio-server-helper-script)

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

> [!WARNING]
> HPC is a shared system! 
> 
> Please do not use `8787` port. This is a common port and may be used by other users.
> 
> Instead, choose a strange and unique port number like `13182`.

You can also edit the SSH config file (`~/.ssh/config`) to add a section for your HPC system:

```
Host hpc
    HostName hpc_address
    User your_username
    LocalForward 8787 localhost:8787
    # change 8787 to a unique port number
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
condatainer install rstudio-server texlive build-essential
```

## Create R Writable Overlay

Creating an ext3 overlay image with R

```bash
# create 30G ext3 image named as env.img under current directory
condatainer overlay -s 30720

# go into the overlay
condatainer e

# INSIDE THE OVERLAY
# install R and required packages
mm-install r-base=4.4 r-tidyverse
# pin the R version to avoid accidental updates
mm-pin r-base
```

## RStudio Server Helper Script

You should always run `rstudio-server` on compute nodes rather than login nodes.

You can use this script: [rstudio-server-helper](./rstudio-server-helper)

It will do the following steps for you:

1. Check if `rstudio-server` is running on any compute node.
2. If yes, establish SSH port forwarding to that node.
3. If not,
   1. Use `sbatch` to submit a job that starts `rstudio-server` on a compute node.
   2. Wait for the job to start and get the node name if job status is `RUNNING`.
   3. Write down the node name and port number to a file.
   4. Set up SSH port forwarding from login node to compute node.

```
Usage: rstudio-server-helper [options]

Options:
  -c <number>     Number of CPUs to allocate (default: 4)
  -m <memory>     Amount of memory to allocate (default: 32G)
  -t <time>       Time limit for the job (default: 12:00:00)
  -p <port>       Port for rstudio-server (default: 8787)
  -e <overlay>    Environment overlay image file (default: env.img)
  -o <overlay>    Additional overlay files (can have -o options)
  -y              Accept all warnings and proceed
  -h              Show this help message
```

Let's set up and run `rstudio-server` on HPC:

```bash
# Download the helper script
# Please make sure $HOME/bin is in your PATH
mkdir -p $HOME/bin
wget https://raw.githubusercontent.com/Justype/condatainer/main/docs/tutorials/rstudio-server-helper -O $HOME/bin/rstudio-server-helper
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
You can run the following command in R to open the project directly
  rstudioapi::openProject("/scratch/your_username/current_working_directory")
```

Then you can go to your local web browser and access `http://localhost:8787`.

And run the provided R command in R console to open project directly.

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
sed -i 's/NCPUS=4/NCPUS=8/' rstudio-server-helper
sed -i 's/MEM=32G/MEM=64G/' rstudio-server-helper
sed -i 's/TIME=12:00:00/TIME=24:00:00/' rstudio-server-helper
sed -i 's/PORT=8787/PORT=13182/' rstudio-server-helper
```
