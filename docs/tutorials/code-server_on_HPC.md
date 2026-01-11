# Run code-server on HPC

## Checklist

- [SLURM job scheduler is available on your HPC system](#slurm-job-scheduler)
- [Have SSH port forwarding set up to the login node](#ssh-port-forwarding)
- [Have CondaTainer installed](#install-condatainer)
- [Have required overlay images created](#install-required-overlays)
- [Have a writable overlay image (optional)](#create-writable-overlay)
- [Have the `code-server-helper` in your PATH](#code-server-helper-script)

Then you can run:

```bash
code-server-helper \
  -p <port_number> \
  -a <auth_password> \
  -e <overlay_image> \
  -c <num_cpus> \
  -m <memory> \
  -t <time_limit>

# Default 
#   port_number=8080
#   auth_password="none" no authentication
#   overlay_image=env.img
#   num_cpus=4
#   memory=32G
#   time_limit=12:00:00
```

Please change the port number to a unique one to avoid conflicts with other users.

If you have preferred settings, you can modify the script directly. See [Change the default setting](#change-the-default-setting) section below.

If you have any issues, see [Common Issues](#common-issues) section below.

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

`code-server` is a web-based application, so you need to access it through a port on the HPC system. 

But you cannot directly access compute nodes from your local machine. To accomplish this, you need to set up SSH port forwarding from your local machine to the HPC login node.

```bash
ssh -L 8080:localhost:8080 your_username@hpc_address
```

After running this command, when your local machine accesses `localhost:8080`, it will be forwarded to `localhost:8080` on the HPC login node.

```{warning}
HPC is a shared system! Please do not use commont port like `8080`.

Instead, choose a strange and unique port number like `13182`.
```

You can also edit the SSH config file (`~/.ssh/config`) to add a section for your HPC system:

```
Host hpc
    HostName hpc_address
    User your_username
    LocalForward 8080 localhost:8080
    # Please change 8080 to a unique port number
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
condatainer install code-server
```

## Create Writable Overlay

Creating an ext3 overlay image

```bash
# create 20G ext3 image named as env.img under current directory
condatainer overlay
# You can -s 10240 to specify the size in MB

# go into the overlay
condatainer e

# INSIDE THE OVERLAY
# install R and required packages
mm-install python=3.11
# pin the Python version to avoid accidental updates
mm-pin python
```

## Code Server Helper Script

I always recommend running `code-server` on compute nodes rather than login nodes.

You can use this scripts: [code-server-helper](https://github.com/Justype/condatainer/blob/main/helpers/code-server-helper)

`code-server-helper` will do the following steps for you:

1. Check if `code-server` is running on any compute node.
2. If yes, establish SSH port forwarding to that node.
3. If not,
   1. Check if the port is available on the login node.
   2. If pass, check the overlay integrity.
   3. If pass, use `sbatch` to submit a job which starts `code-server`.
   4. Wait for the job to start and get the compute node name.
   5. Record the job ID, node name, port numer, and working directory.
   6. Set up SSH port forwarding from login node to compute node.

```
Usage: code-server-helper [options]

Options:
  -c <number>     Number of CPUs to allocate (default: 4)
  -m <memory>     Amount of memory to allocate (default: 16G)
  -t <time>       Time limit for the job (default: 12:00:00)
  -p <port>       Port for code-server (default: 8080)
  -a <auth>       Password for code-server authentication (default: none)
  -e <overlay>    Environment overlay image file (default: env.img)
  -o <overlay>    Additional overlay files (can have multiple -o options)
  -v              View Mode NCPUS:1 MEM:4G TIME:02:00:00
  -y              Accept all warnings and proceed
  -h              Show this help message
```

Let's set up and run `code-server` on HPC:

```bash
# Download the helper script
# Please make sure $HOME/bin is in your PATH
mkdir -p $HOME/bin
wget https://raw.githubusercontent.com/Justype/condatainer/main/helpers/code-server-helper -O $HOME/bin/code-server-helper
chmod +x $HOME/bin/code-server-helper
```

Then you can run the script: 

```bash
code-server-helper
```

After running the script, you will see output like this:

```
Waiting for job 1299123 to start running. Current status: PENDING.
Job 1299123 is now running on node cm23.
code-server at http://localhost:$PORT?folder=/scratch/your_username/current_working_directory
If you want to stop it, run: scancel 1299123
```

Then you can go to your local web browser and access `http://localhost:8080`. or click the link provided to open the project directly.

## Change the default setting

You can modify the default settings in the script, such as CPU, memory, time limit, overlay image path, port number.

They are at the beginning of the script:

```bash
#!/bin/bash

NCPUS=4
MEM=32G
TIME=12:00:00
PORT=8080
AUTH="none"
```

You can use editors like `vim` or `nano` to change these values.

or use `sed`

```bash
sed -i 's/NCPUS=4/NCPUS=8/' code-server-helper
sed -i 's/MEM=32G/MEM=64G/' code-server-helper
sed -i 's/TIME=12:00:00/TIME=24:00:00/' code-server-helper
sed -i 's/PORT=8080/PORT=13182/' code-server-helper
sed -i 's/AUTH="none"/AUTH="your_password"/' code-server-helper
```

## Common Issues

### Too many files under `.local/share/code-server/extensions`

If you have installed many extensions, the number of files under `.local/share/code-server/extensions` may exceed the quota.

To fix this, you can create a symbolic link to the writable overlay image.

```bash
# Inside the overlay
mkdir -p /ext3/home/.local/share/code-server
mv $HOME/.local/share/code-server/extensions /ext3/home/.local/share/code-server/extensions
ln -s /ext3/home/.local/share/code-server/extensions $HOME/.local/share/code-server/extensions
```

A new issue is that if you change the overlay image, you will lose all installed extensions. But at least you won't hit the quota limit.
