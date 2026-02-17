# Run VS Code Tunnel on HPC

Run following commands to start VS Code Tunnel on your HPC system:

```bash
# download/update helper scripts
condatainer helper -u
# create a 20G writable overlay (optional but recommended)
condatainer o -s 20g
# start vscode-tunnel
condatainer helper vscode-tunnel
```

VS Code Tunnel connects your local VS Code to a remote compute node **without requiring SSH port forwarding**. It uses Microsoft's relay service to establish the connection.

## Checklist

- [A supported scheduler is available on your HPC system](#scheduler--workload-manager)
- [Have CondaTainer installed](#install-condatainer)
- [Have a writable overlay image (optional)](#create-writable-overlay)
- [Check the Script Parameters](#vs-code-tunnel-helper-script)

Then you can run:

```bash
condatainer helper vscode-tunnel
# Default
#   -a auth_option=github (or microsoft)
#   -n machine_name=<username>-<hostname> (truncated to 20 chars)
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

## SSH Port Forwarding is Not Required

VS Code Tunnel directly connects to the VS Code server instead of using SSH port forwarding.

```
Your machine <-----> VS Code server <-----> HPC Compute Node
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

## VS Code Tunnel Helper Script

The script will do the following steps for you:

1. Download the VS Code CLI if not available.
2. Check if a VS Code tunnel is running on any compute node.
3. If yes, print helpful information to connect to that tunnel.
4. If not,
   1. Check the overlay integrity.
   2. Submit the SLURM job to start `vscode-tunnel`.
   3. When the job starts, display authentication information.

Then you can either:
- Use [VS Code Remote - Tunnels extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode.remote-server) to connect
- Use https://vscode.dev/

```
Usage: vscode-tunnel [options]

Options:
  -c <number>     Number of CPUs to allocate (default: 4)
  -m <memory>     Amount of memory to allocate (default: 16G)
  -t <time>       Time limit for the job (default: 12:00:00)
  -g <gpu>        GPU resources (e.g., 1, a100:2). Empty means no GPU.
  -v              View Mode NCPUS:1 MEM:4G TIME:02:00:00

  -a <provider>   github or microsoft for authentication (default: github)
  -n <name>       Machine name for the tunnel (default: <username>-<hostname>)
  -b <image>      Base image file
  -e <overlay>    Environment overlay image file (default: env.img)
  -o <overlay>    Additional overlay files (can have multiple -o options)

  config          Show config file path and contents
```

## Configuration

The script saves its defaults to `$XDG_CONFIG_HOME/condatainer/helper/vscode-tunnel` (defaults to `~/.config/condatainer/helper/vscode-tunnel`) on first run. Subsequent runs load from this file.

```bash
# View current config
condatainer helper vscode-tunnel config

# Reset to defaults (delete config file, next run recreates it)
rm ~/.config/condatainer/helper/vscode-tunnel
```

Config files use simple `KEY="VALUE"` format and can be edited directly:

```bash
NCPUS="4"
MEM="16G"
TIME="12:00:00"
AUTH="github"
MACHINE_NAME="username-hpc"
OVERLAY="env.img"
```

## Running VS Code Tunnel

Let's set up and run VS Code Tunnel on HPC:

```bash
# Download the helper scripts
condatainer helper -u
```

Then you can run the script:

```bash
condatainer helper vscode-tunnel
```

After running the script, you will see output like this:

```
Waiting for job 6628611 to start running. Current status: PENDING.
[INFO] Waiting for VS Code Tunnel to initialize...

Please use the following message for authentication:
------------------------------------------------
To sign in, use a web browser to open the page https://microsoft.com/devicelogin and enter the code FE2G6WJQK to authenticate.
------------------------------------------------
Waiting for authentication...
```

Then you need to open the URL in your web browser and enter the code to authenticate. (Same as the GitHub authentication flow.)

After authentication, you will see output like this:

```
[PASS] Tunnel is ready!
To connect, you can either:
  - Use VS Code Remote - Tunnels extension
  - Use this link: https://vscode.dev/tunnel/username-hpc_name/scratch/username/playground
And use your Microsoft account to connect to username-hpc_name
CWD: /scratch/username/playground
```

Now you can either use the [VS Code Remote - Tunnels extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode.remote-server) or open the link in your web browser to connect to the VS Code server running on the compute node.

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

### Machine name already taken

If you see an error about the machine name being taken, edit your config file to use a different name:

```bash
condatainer helper vscode-tunnel config
# Note the config file path, then edit it
# Change MACHINE_NAME to something unique
```

### Authentication timeout

If authentication times out (after 10 minutes), the job will be cancelled. Run the script again to retry.

## Comparison with Other Options

| Feature | VS Code Tunnel | VS Code Server |
|---------|----------------|----------------|
| Port forwarding | Not required | Required |
| Extensions | All (including Copilot) | All (including Copilot) |
| Authentication | GitHub/Microsoft account | Token-based |
| Performance | Relay through Microsoft | Direct connection |
| Best for | When port forwarding is difficult | Full control, no external relay |

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
