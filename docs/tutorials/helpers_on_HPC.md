# Helper Scripts on HPC

Helper scripts automate the workflow of running web-based app (VS Code, RStudio, Jupyter, etc.) on HPC compute nodes.

This tutorial explains the core concepts behind helper scripts and how they simplify running interactive applications on HPC.

[Available Applications](#available-applications)

## Main Concepts

1. [Project-Based Workspace Overlays](#project-based-workspace-overlays) - Each project has its own isolated environment
2. [Helper Scripts](#helper-scripts) - Automate job submission and port forwarding
3. [SSH Port Forwarding](#ssh-port-forwarding) - Securely access the clusters

### Project-Based Workspace Overlays

**Best Practice:** Create an `env.img` overlay in each project directory.

```bash
cd $SCRATCH/projects/project-a
condatainer o -s 20g        # Creates env.img here

cd $SCRATCH/projects/project-b
condatainer o -s 30g        # Creates env.img here
```

This approach provides:
- **Isolated environments** for each project
- **Portable workspaces** - overlay stays with project files
- **Avoids inode limits** - one file instead of thousands

See [Workspace Overlays](../user_guide/workspace_overlays.md) for more details.

### Helper Scripts

Helper scripts automate the multi-step process of running web applications on HPC:

1. Check dependencies and overlay integrity
2. Auto-install required module overlays if missing
3. Submit jobs to the scheduler (SLURM, etc.) and set up port forwarding
4. Reconnect to existing sessions automatically

```bash
# Download/update helper scripts
condatainer helper -u

# Run a helper script (uses env.img in current directory by default)
condatainer helper vscode-server -p 13182
```

### SSH Port Forwarding

Helper scripts will establish SSH port forwarding from the login node to the compute node. To access the web app, you need to set up port forwarding from your local machine to the HPC login node.

**Setup:**

1. On your local machine, configure SSH port forwarding (`~/.ssh/config`):

```
Host hpc
    HostName hpc.university.edu
    User your_username
    LocalForward 13182 localhost:13182
```

2. Connect with: `ssh hpc`

3. On HPC, run the helper:

```bash
cd ~/projects/my-project
condatainer helper vscode-server -p 13182
```

4. Open in browser: `http://localhost:13182`

**How it works:**

```
Local Browser → localhost:13182 (local)
                     ↓ SSH tunnel
                localhost:13182 (HPC login node)
                     ↓ automatic SSH (by helper)
                localhost:13182 (compute node)
```

```{warn}
Make sure to choose a unique port number for each session to avoid conflicts. Avoid common ports like 8080, 8787, etc.
```

## Common Workflow

### Starting a New Project

```bash
# 1. Create project directory
mkdir $SCRATCH/my-analysis
cd $SCRATCH/my-analysis

# 2. Create overlay for this project
condatainer o -s 20g

# 3. Install dependencies
condatainer e -- mm-install python=3.11 numpy pandas scikit-learn

# 4. Start VS Code Server
# -w sets working directory to pwd
condatainer helper vscode-server -p 13182 -w
```

### Working on Existing Project

```bash
# Navigate to project (helper automatically finds env.img here)
cd $SCRATCH/my-analysis

# Start the web app
# -w sets working directory to pwd
condatainer helper vscode-server -p 13182 -w
```

```{note}
If a previous session is still active, the helper will reconnect you to it instead of starting a new one.

If the session has ended, it can reuse your previous settings based on your [REUSE_MODE](#reuse-mode).
```

## Resource Allocation

### Common Options

All helper scripts support:

```bash
-c <number>     # CPUs (default: 4)
-m <memory>     # Memory (default: 16G)
-t <time>       # Time limit (default: 12:00:00)
-g <gpu>        # GPU (e.g., 1, a100:2)
-v              # View mode: 1 CPU, 4G, 2h

-p <port>       # Port number (for port-forwarding apps)
-w              # Use current directory as working directory (clears additional overlays)

-e <.img>       # Writable overlay (default: env.img)
-o <overlay>    # Additional overlay (can use multiple -o)
```

```{note}
Using `-w` will set the working directory to your current directory and clear any additional overlays (`-o`) from previous runs. This is useful when you've `cd` to a new project directory and want to start fresh.
```

### Examples

```bash
# Light workload (view mode)
condatainer helper vscode-server -v

# Standard workload
condatainer helper vscode-server -c 4 -m 16G -p 13182 -t 12:00:00

# Heavy computation with GPU
condatainer helper vscode-server -c 14 -m 250G -g h100 -t 48:00:00
```

### Check Available Resources

```bash
condatainer scheduler
```

Shows scheduler type, max resources, and available GPUs.

Example output:

```
Scheduler Information:
  Type:      SLURM
  Binary:    /opt/software/slurm/current/bin/sbatch
  Version:   24.11.7
  Status:    Available

The scheduler is available and ready for job submission.

Max Resource Limits:
  Max CPUs:   44
  Max Memory: 184000 MB (180 GB)
  Max Time:   7d

Available GPUs:
  h100: 2144 total, 416 available
  mi300a: 4 total, 0 available
  nvidia_h100_80gb_hbm3_1g.10gb: 2264 total, 728 available
  nvidia_h100_80gb_hbm3_2g.20gb: 768 total, 0 available
  nvidia_h100_80gb_hbm3_3g.40gb: 768 total, 0 available
  t4: 28 total, 0 available
```

## Configuration and Defaults

Helper scripts save your preferences to `~/.config/condatainer/helper/defaults/<app-name>`:

```bash
# View current config
condatainer helper vscode-server config

# Edit config file directly
nano ~/.config/condatainer/helper/defaults/vscode-server

# Reset to defaults
rm ~/.config/condatainer/helper/defaults/vscode-server
```

Example config:
```bash
NCPUS="8"
MEM="32G"
TIME="24:00:00"
PORT="13182"
OVERLAY="env.img"
REUSE_MODE="ask"    # Control reuse behavior: 'always', 'ask', or 'never'
```

### Reuse Mode

When a previous job has ended, the helper can reuse its settings (CPUs, memory, time, GPU, port, working directory, overlays) so you don't need to re-type everything.

The `REUSE_MODE` setting controls this behavior:

- **`ask`** (default): Shows previous settings and prompts to reuse
- **`always`**: Automatically reuses without prompting
- **`never`**: Always uses defaults from config

```bash
# Edit config file
nano ~/.config/condatainer/helper/defaults/vscode-server

# Set to auto-reuse previous settings:
REUSE_MODE="never"
```

You can override specific reused settings with flags:

```bash
# Reuse previous settings but change GPU
condatainer helper vscode-server -g h100 -c 14 -m 200g

# Reuse but clear GPU
condatainer helper vscode-server -g ''

# Reuse but switch to current directory and clear previous overlays
condatainer helper vscode-server -w
```

```{note}
You will be still prompted to reuse even you enter `-w` or `-g ''`, since these flags modify the reused settings.

To skip the prompt entirely, set `REUSE_MODE="always"` in your config.
```

## Common Issues

### Port Already in Use

```bash
# Let script auto-select a port
condatainer helper vscode-server

# Or choose a different port
condatainer helper vscode-server -p 13183
```

### Cannot Connect

**Checklist:**
1. Job is running: `squeue -u $USER`
2. SSH port forwarding is active: `ssh -L 13182:localhost:13182 user@hpc`
3. Correct URL in browser: `http://localhost:13182`
4. Check job logs: `cat ~/logs/<helper-name>-<job_id>.log`
   - Example: `cat ~/logs/vscode-server-1234567.log`

### Overlay Integrity Issues

Helper scripts automatically check and repair overlays. If manual repair needed:

```bash
e2fsck -f env.img
```

## Available Applications

### VS Code

- [VS Code Server](./vscode-server_on_HPC.md) - Full VS Code experience (recommended)
- [VS Code Tunnel](./vscode-tunnel_on_HPC.md) - No port forwarding needed
- [code-server](./code-server_on_HPC.md) - Open-source alternative

### RStudio

- [RStudio Server](./rstudio-server_on_HPC.md) - Using Posit R docker images
- [RStudio Server (Conda)](./rstudio-server-conda_on_HPC.md) - Using Conda R
