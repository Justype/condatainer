# Run Code Server on HPC

Run following commands to start `code-server` on your HPC system:

```bash
# download/update helper scripts
condatainer helper -u
# create a 20G writable overlay (optional but recommended)
condatainer o -s 20g
# start code-server
condatainer helper code-server
```

```{note}
You cannot use Microsoft extensions like Pylance and Copilot with `code-server`. Use [VS Code Server](./vscode-server_on_HPC.md) instead for full VS Code experience.
```

## Checklist

- [A supported scheduler is available on your HPC system](#scheduler--workload-manager)
- [Have SSH port forwarding set up to the login node](#ssh-port-forwarding)
- [Have CondaTainer installed](#install-condatainer)
- [Have a writable overlay image (optional)](#create-writable-overlay)
- [Check the Script Parameters](#code-server-helper-script)

Then you can run:

```bash
condatainer helper code-server
# Default
#   -p port_number=auto-selected if omitted
#   -a auth_password="none" or a password
#   -c num_cpus=4
#   -m memory=16G
#   -t time_limit=12:00:00
```

You can set [Config](./helpers_on_HPC.md#configuration-and-defaults) and helpers also support [Reusing previous settings](./helpers_on_HPC.md#reuse-mode).

If you have any issues, see [Common Issues](#common-issues) section below.

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

Run the following command to install CondaTainer if it is not installed:

```bash
curl -fsSL https://get-condatainer.justype.net | bash
```

## Create Writable Overlay

Creating an ext3 overlay image (optional but recommended for persistent settings):

```bash
# create a 20G ext3 image named `env.img` in the current working directory
condatainer o -s 20g

# go into the overlay
condatainer e
```

```bash
# INSIDE THE OVERLAY
# install Python and required packages
mm-install python=3.11
```

See [Launch a Shell within the Workspace Overlay](../user_guide/workspace_overlays.md#launch-a-shell-within-the-workspace-overlay) for more details.

## VS Code Server Helper Script

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
  -g <gpu>        GPU resources (e.g., 1, a100:2). Empty means no GPU.
  -v              View Mode NCPUS:1 MEM:4G TIME:02:00:00

  -p <port>       Port for vscode-server (default: randomly picked). Valid range: 1024-65535.
  -a <auth>       Password for code-server authentication (default: none)
  -b <image>      Base image file
  -e <overlay>    Environment overlay image file (default: env.img)
  -o <overlay>    Additional overlay files (can have multiple -o options)

  config          Show config file path and contents
```

## Configuration

The script saves its defaults to `~/.config/condatainer/helper/defaults/code-server` on first run. Subsequent runs load from this file.

```bash
# View current config
condatainer helper code-server config

# Reset to defaults (delete config file, next run recreates it)
rm ~/.config/condatainer/helper/defaults/code-server
```

## Running `code-server`

Let's set up and run VS Code Server on HPC:

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

Now you can click the link provided to open the project in your browser.

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

You can use [VS Code Tunnel](./vscode-tunnel_on_HPC.md) or [VS Code Server](./vscode-server_on_HPC.md) instead.

### Port already in use

If your previously selected port is already in use, the script will notify you. You can either:
- Wait for the other process to finish
- Choose a different port with `-p <new_port>`

## Comparison with Other Options

| Feature | VS Code Server | VS Code Tunnel | code-server |
|---------|---------------|----------------|-------------|
| Port forwarding | Required | Not required | Required |
| Extensions | All (Pylance, Copilot) | All | Open VSX only |
| Authentication | Token-based | Microsoft/GitHub account | Password/none |
| Performance | Direct connection | Relay through Microsoft | Direct connection |

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
