# Run VS Code Server on HPC

Run following commands to start VS Code Server on your HPC system:

```bash
# download/update helper scripts
condatainer helper -u
# create a 20G writable overlay (optional but recommended)
condatainer o -s 20g
# start vscode-server
condatainer helper vscode-server
```

VS Code Server provides the **full VS Code experience** including all extensions (Pylance, Copilot, etc.) in your browser.

## Checklist

- [A supported scheduler is available on your HPC system](#scheduler--workload-manager)
- [Have SSH port forwarding set up to the login node](#ssh-port-forwarding)
- [Have CondaTainer installed](#install-condatainer)
- [Have a writable overlay image (optional)](#create-writable-overlay)
- [Check the Script Parameters](#vscode-server-helper-script)

Then you can run:

```bash
condatainer helper vscode-server
# Default
#   -p port_number=auto-selected if omitted
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

```{note}
Slurm is fully tested. Other schedulers are experimentally supported. Please report any issues you encounter.
```

## SSH Port Forwarding

```
Your machine -----> HPC Login Node -----> HPC Compute Node
        (port forwarding)     (port forwarding)
```

VS Code Server is a web-based application, so you need to access it through a port on the HPC system.

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

The script will do the following steps for you:

1. Download the VS Code CLI if not available.
2. Check if VS Code Server is running on any compute node.
3. If yes, establish SSH port forwarding to that node.
4. If not,
   1. Check port and overlay integrity.
   2. Submit the SLURM job to start VS Code Server.
   3. When the job starts, record and set up SSH port forwarding.

```
Usage: vscode-server [options]

Options:
  -c <number>     Number of CPUs to allocate (default: 4)
  -m <memory>     Amount of memory to allocate (default: 16G)
  -t <time>       Time limit for the job (default: 12:00:00)
  -g <gpu>        GPU resources (e.g., 1, a100:2). Empty means no GPU.
  -v              View Mode NCPUS:1 MEM:4G TIME:02:00:00

  -p <port>       Port for vscode-server (default: randomly picked). Valid range: 1024-65535.
  -a <token>      Connection token for the web UI (if empty, one is generated)
  -b <image>      Base image file
  -e <overlay>    Environment overlay image file (default: env.img)
  -o <overlay>    Additional overlay files (can have multiple -o options)

  config          Show config file path and contents
```

## Configuration

The script saves its defaults to `$XDG_CONFIG_HOME/condatainer/helper/vscode-server` (defaults to `~/.config/condatainer/helper/vscode-server`) on first run. Subsequent runs load from this file.

```bash
# View current config
condatainer helper vscode-server config

# Reset to defaults (delete config file, next run recreates it)
rm ~/.config/condatainer/helper/vscode-server
```

## Running VS Code Server

Let's set up and run VS Code Server on HPC:

```bash
# Download the helper scripts
condatainer helper -u
```

Then you can run the script:

```bash
condatainer helper vscode-server
```

After running the script, you will see output like this:

```
Waiting for job 1299123 to start running. Current status: PENDING.
Job 1299123 is now running on node cm23.
vscode-server at http://localhost:<port>?tkn=<token>&folder=/scratch/your_username/current_working_directory
If you want to stop it, run: scancel 1299123
```

Now you can click the link provided to open the project in your browser.

Don't forget to set up SSH port forwarding from your local machine to the HPC login node before accessing the link.

## Common Issues

### Too many files under `.vscode-server/extensions`

If you have installed many extensions, the number of files under `.vscode-server/extensions` may exceed the quota.

To fix this, you can create a symbolic link to the writable overlay image.

```bash
# Inside the overlay
mkdir -p /ext3/home/.vscode-server
mv $HOME/.vscode-server/extensions /ext3/home/.vscode-server/extensions
ln -s /ext3/home/.vscode-server/extensions $HOME/.vscode-server/extensions
```

A new issue is that if you change the overlay image, you will lose all installed extensions. But at least you won't hit the quota limit.

### Connection token issues

If you're having trouble connecting, ensure you're using the correct token from the output. The token is included in the URL query parameter `tkn=<token>`.

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
