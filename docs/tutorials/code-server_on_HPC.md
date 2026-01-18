# Run code-server on HPC

Run following commands to start `code-server` on your HPC system:

```bash
condatainer helper -u
condatainer helper code-server
```

## Checklist

- [SLURM job scheduler is available on your HPC system](#slurm-job-scheduler)
- [Have SSH port forwarding set up to the login node](#ssh-port-forwarding)
- [Have CondaTainer installed](#install-condatainer)
- [Have required overlay images created](#install-required-overlays)
- [Have a writable overlay image (optional)](#create-writable-overlay)
- [Check the Script Parameters](#code-server-helper-script)

Then you can run:

```bash
condatainer helper code-server \
# Default
#   -p port_number=auto-selected if omitted
#   -a auth_password="none" or a password
#   -c num_cpus=4
#   -m memory=16G
#   -t time_limit=12:00:00
```

You can create alias in your shell config file (`~/.bashrc` or `~/.zshrc`):

```bash
# Change 13182 to your preferred port number
alias code-server-start='condatainer helper code-server -p 13182'
```

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
ssh -L <local_port>:localhost:<remote_port> your_username@hpc_address
```

After running this command, when your local machine accesses `localhost:<local_port>`, it will be forwarded to `localhost:<remote_port>` on the HPC login node.

```{warning}
HPC is a shared system! Please do not use a common port like `8080`.

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
condatainer install code-server
```

## Create Writable Overlay

Creating an ext3 overlay image

```bash
# create a 20G ext3 image named `env.img` in the current working directory
condatainer overlay
# You can use `-s 10240` to specify the size in MiB

# go into the overlay
condatainer e

# INSIDE THE OVERLAY
# install R and required packages
mm-install python=3.11
# pin the Python version to avoid accidental updates
mm-pin python
```

## Code Server Helper Script

`code-server` will do the following steps for you:

1. Check if `code-server` is running on any compute node.
2. If yes, establish SSH port forwarding to that node.
3. If not,
   1. Check port and overlay integrity.
   2. Submit the SLURM job to start `code-server`.
   3. When the job starts, record and set up SSH port forwarding.

```
Usage: code-server [options]

Options:
  -c <number>     Number of CPUs to allocate (default: 4)
  -m <memory>     Amount of memory to allocate (default: 16G)
  -t <time>       Time limit for the job (default: 12:00:00)
  -p <port>       Port for code-server (if omitted, an available port will be chosen)
  -a <auth>       Password for code-server authentication (default: none)
  -o <overlay>    Additional overlay files (can have multiple -o options)
  -v              View Mode NCPUS:1 MEM:4G TIME:02:00:00
  -h              Show this help message
```

Let's set up and run `code-server` on HPC:

```bash
# Download the helper scripts
condatainer helper -u
```

Then you can run the script: 

```bash
condatainer helper code-server
```

After running the script, you will see output like this:

```
Waiting for job 1299123 to start running. Current status: PENDING.
Job 1299123 is now running on node cm23.
code-server at http://localhost:<port>?folder=/scratch/your_username/current_working_directory
If you want to stop it, run: scancel 1299123
```

Now you can click the link provided to open the project.

Don't forget to set up SSH port forwarding from your local machine to the HPC login node before accessing the link.

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

### Pylance not available

Pylance is owned by Microsoft and is not open source. So it is not included in the [open-vsx registry](https://open-vsx.org/), which is used by `code-server`.

Solution: Use Pyright instead. see https://open-vsx.org/extension/ms-pyright/pyright

### Set up GitHub Copilot

Sadly, GitHub Copilot is not working with `code-server`.

Available workarounds:

- [sunpix/howto-install-copilot-in-code-server](https://github.com/sunpix/howto-install-copilot-in-code-server)
- [Code Server Discussion 5063](https://github.com/coder/code-server/discussions/5063)

I tried both methods, but none worked for me.

You can use [VS Code Tunnel](./vscode-tunnel_on_HPC.md) instead.
