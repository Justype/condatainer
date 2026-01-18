# Run VS Code Tunnel on HPC
Run following commands to start `vscode-tunnel` on your HPC system:

```bash
condatainer helper -u
condatainer helper vscode-tunnel
```

## Checklist

- [SLURM job scheduler is available on your HPC system](#slurm-job-scheduler)s
- [Have CondaTainer installed](#install-condatainer)
- [Have a writable overlay image (optional)](#create-writable-overlay)
- [Check the Script Parameters](#vs-code-tunnel-helper-script)

Then you can run:

```bash
condatainer helper vscode-tunnel \
    -a <auth_option> \
    -c <num_cpus> \
    -m <memory> \
    -t <time_limit>

# Default 
#   -a auth_option="microsoft" or "github"
#   -c num_cpus=4
#   -m memory=16G
#   -t time_limit=12:00:00
```

If you have any issues, see [Common Issues](#common-issues) section below.

## SLURM Job Scheduler

Run the following command to check if SLURM is available on your HPC system:

```bash
sbatch --version
```

## SSH port forwarding is not required

VS Code CLI will directly connect to the VS Code server instead of using SSH port forwarding.

```
Your machine <-----> VS Code server <-----> HPC Compute Node
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

## No required overlay images

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

## VS Code Tunnel Helper Script

`vscode-tunnel` will do the following steps for you:

1. Download the `code` CLI if not available.
2. Check if a VS Code tunnel is running on any compute node.
3. If yes, print helpful information to connect to that tunnel.
4. If not,
   2. Check the overlay integrity.
   2. Submit the SLURM job to start `vscode-tunnel`.
   3. When the job starts, record and print the authentication information.

Then you can either:
- Use [VS Code Remote - Tunnels extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode.remote-server) to connect
- Use https://vscode.dev/

```
Usage: vscode-tunnel [options]

Options:
  -c <number>     Number of CPUs to allocate (default: 4)
  -m <memory>     Amount of memory to allocate (default: 16G)
  -t <time>       Time limit for the job (default: 12:00:00)
  -a <provider>   microsoft or github for authentication (default: microsoft)
  -v              View Mode NCPUS:1 MEM:4G TIME:02:00:00
  -h              Show this help message
```

Let's set up and run `vscode-tunnel` CLI on HPC:

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
[SUCCESS] Tunnel is ready!
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
