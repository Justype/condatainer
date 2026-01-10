# How to use vs code on HPC with CondaTainer 

## Checklist

- [ ] [SLURM job scheduler is available on your HPC system](#slurm-job-scheduler)
- [ ] [Have CondaTainer installed](#install-condatainer)
- [ ] [Have a writable overlay image (optional)](#create-writable-overlay)
- [ ] [Have the `vscode-tunnel-helper` in your PATH](#vs-code-tunnel-helper-script)

Then you can run:

```bash
vscode-tunnel-helper \
  -a <auth_option> \
  -c <num_cpus> \
  -m <memory> \
  -t <time_limit>

# Default 
#   auth_option="microsoft" or "github"
#   overlay_image=env.img
#   num_cpus=4
#   memory=32G
#   time_limit=12:00:00
```

If you have preferred settings, you can modify the script directly. See [Change the default setting](#change-the-default-setting) section below.

If you have any issues, see [Common Issues](#common-issues) section below.

## SLURM Job Scheduler

Run the following command to check if SLURM is available on your HPC system:

```bash
sbatch --version
```

## SSH Port Forwarding Is not Required

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
curl -sL https://raw.githubusercontent.com/Justype/condatainer/main/assets/install.sh | bash
```

## No required overlay images

Creating required overlay images is not necessary.

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

## VS Code Tunnel Helper Script

Always run VS Code CLI on compute nodes rather than login nodes.

You can use this scripts: [vscode-tunnel-helper](./vscode-tunnel-helper)

`vscode-tunnel-helper` will do the following steps for you:

1. Download the `code` CLI if not available.
2. Check if a VS Code tunnel is running on any compute node.
3. If yes, print helpful information to connect to that tunnel.
4. If not,
   1. Check if the port is available on the login node.
   2. If pass, check the overlay integrity.
   3. If pass, use `sbatch` to submit a job which starts `code` CLI.
   4. Wait for the job to start.
   5. Print the authentication URL to connect to the VS Code server.

You can either:
- Have [VS Code Remote - Tunnels extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode.remote-server) installed
- Use https://vscode.dev/ and the tunnel

```
Usage: vscode-tunnel-helper [options]

Options:
  -c <number>     Number of CPUs to allocate (default: 4)
  -m <memory>     Amount of memory to allocate (default: 16G)
  -t <time>       Time limit for the job (default: 12:00:00)
  -a <provider>   microsoft or github for authentication (default: microsoft)
  -e <overlay>    Environment overlay image file (default: env.img)
  -o <overlay>    Additional overlay files (can have -o options)
  -y              Accept all warnings and proceed
  -h              Show this help message
```

Let's set up and run `code` CLI on HPC:

```bash
# Download the helper script
# Please make sure $HOME/bin is in your PATH
mkdir -p $HOME/bin
wget https://raw.githubusercontent.com/Justype/condatainer/main/docs/tutorials/vscode-tunnel-helper -O $HOME/bin/vscode-tunnel-helper
chmod +x $HOME/bin/vscode-tunnel-helper
```

Then you can run the script: 

```bash
vscode-tunnel-helper
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

Then you need to open the URL in your web browser and enter the code to authenticate. (Same as the github authentication flow.)

After authentication, you will see output like this:

```
[SUCCESS] Tunnel is ready!
To connect, you can either:
  - Use VSCode Remote - Tunnels extension
  - Use this link: https://vscode.dev/tunnel/username-hpc_name/scratch/username/playground
And use your microsoft account to connect to username-hpc_name
CWD: /scratch/username/playground
```

Then you can either use the [VSCode Remote - Tunnels extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode.remote-server) or open the link in your web browser to connect to the VS Code server running on the compute node.

Personally, I prefer use local extension.

## Change the default setting

You can modify the default settings in the script, such as CPU, memory, time limit, overlay image path, port number.

They are at the beginning of the script:

```bash
#!/bin/bash

NCPUS=4
MEM=16G
TIME=12:00:00
AUTH="microsoft"
```

You can use editors like `vim` or `nano` to change these values.

or use `sed`

```bash
sed -i 's/NCPUS=4/NCPUS=8/' vscode-tunnel-helper
sed -i 's/MEM=16G/MEM=64G/' vscode-tunnel-helper
sed -i 's/TIME=12:00:00/TIME=24:00:00/' vscode-tunnel-helper
sed -i 's/AUTH="microsoft"/AUTH="github"/' vscode-tunnel-helper
```

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
