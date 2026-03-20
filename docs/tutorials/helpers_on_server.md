# Helper Scripts on Headless Servers

Helper scripts automate the workflow of running web-based applications (VS Code, RStudio, noVNC, etc.) directly on a headless server without HPC schedulers.

This tutorial explains how to use helper scripts on servers where you have direct access (no job scheduler required).

[Available Applications](#available-applications)

## Main Concepts

1. [Project-Based Workspace Overlays](#project-based-workspace-overlays) - Each project has its own isolated environment
2. [Helper Scripts](#helper-scripts) - Launch applications directly on the server
3. [SSH Port Forwarding](#ssh-port-forwarding) - Securely access applications from your local machine

### Project-Based Workspace Overlays

**Best Practice:** Create an `env.img` overlay in each project directory.

```bash
cd ~/projects/project-a
condatainer o -s 20g        # Creates env.img here

cd ~/projects/project-b
condatainer o -s 30g        # Creates env.img here
```

This approach provides:
- **Isolated environments** for each project
- **Portable workspaces** - overlay stays with project files
- **Avoids inode limits** - one file instead of thousands

See [Workspace Overlays](../user_guide/workspace_overlays.md) for more details.

### Helper Scripts

Helper scripts automate the process of running web applications on your server:

1. Check dependencies and overlay integrity
2. Auto-install required module overlays if missing
3. Launch applications directly and handle port conflicts

```bash
# Download/update helper scripts
condatainer helper -u

# Run a helper script (uses -w to set working directory)
condatainer helper vscode-server -w -p 13182
```

### SSH Port Forwarding

To access web applications running on your server from your local browser, you need to set up SSH port forwarding from your local machine to the server.

**Setup:**

1. On your local machine, configure SSH port forwarding (`~/.ssh/config`):

```
Host myserver
    HostName server.example.com
    User your_username
    LocalForward 13182 localhost:13182
```

2. Connect with: `ssh myserver`

3. On the server, run the helper:

```bash
cd ~/projects/my-project
condatainer helper vscode-server -w -p 13182
```

4. Open in browser: `http://localhost:13182`

**How it works:**

```
Local Browser → localhost:13182 (local)
                     ↓ SSH tunnel
                localhost:13182 (server)
```

```{warning}
Make sure to choose a unique port number for each session to avoid conflicts. Avoid common ports like 8080, 8787, etc.
```

## Common Workflow

### Starting a New Project

```bash
# 1. Create project directory
mkdir ~/projects/my-analysis
cd ~/projects/my-analysis

# 2. Create overlay for this project
condatainer o -s 20g -- python=3.11 numpy pandas scikit-learn

# 3. Install dependencies (optional)
condatainer e -- mm-install seaborn scipy

# 4. Start VS Code Server
condatainer helper vscode-server -w -p 13182
```

```{note}
`o -- <packages>` will create the overlay on the node local SSD, which is much faster.

The latter `e -- mm-install <packages>` will be slower since it modifies the overlay on the network filesystem.

Therefore, it's best to install as many dependencies as possible in the initial `condatainer o` step. You can always add more later with `condatainer e` if needed.
```

### Working on Existing Project

```bash
# Navigate to project (helper automatically finds env.img here)
cd ~/projects/my-analysis

# Start the web app
condatainer helper vscode-server -w -p 13182
```

```{note}
If the port is already in use, the script will warn you and exit. Either kill the existing process or choose a different port with `-p`.
```

## Common Options

Helper scripts support:

```bash
-p <port>       # Port number (default: auto-selected)
-w              # Use current directory as working directory (clears previous overlays)

-e <.img>       # Writable overlay (default: env.img)
-o <overlay>    # Additional overlay (can use multiple -o)
-b <image>      # Base image file
```

When `-e` is `env.img` (the default), the helper searches for it in order:

1. `env.img` (current directory)
2. `overlay/env.img`
3. `src/overlay/env.img`

```{note}
Using `-w` will set the working directory to your current directory and clear any additional overlays (`-o`) from previous runs. This is useful when you've `cd` to a new project directory and want to start fresh.
```

### Examples

```bash
# Auto-select an available port
condatainer helper vscode-server -w

# Use specific port
condatainer helper vscode-server -w -p 13182

# Use custom overlay
condatainer helper vscode-server -w -p 13182 -e my-env.img

# Add additional overlays
condatainer helper jupyter -w -p 13182 -o pytorch -o scikit-learn
```

## Configuration and Defaults

Helper scripts save your preferences to `$XDG_CONFIG_HOME/condatainer/helper/<app-name>` (defaults to `~/.config/condatainer/helper/<app-name>`):

```bash
# View current config
condatainer helper vscode-server config

# Edit config file directly
nano ~/.config/condatainer/helper/vscode-server

# Reset to defaults
rm ~/.config/condatainer/helper/vscode-server
```

Example config:
```bash
PORT="13182"
OVERLAY="env.img"
BASE_IMAGE=""
TOKEN=""
REUSE_MODE="ask"   # 'always', 'ask' (default), or 'never'
```

### Reuse Mode

When a previous job has ended, the helper can reuse its settings (CPUs, memory, time, GPU, port, working directory, overlays) so you don't need to re-type everything.

The `REUSE_MODE` setting controls this behavior:

- **`ask`** (default): Shows previous settings and prompts:
  - `Y` — reuse all previous settings (CPUs, memory, time, GPU, port, working directory, overlays)
  - `n` — discard previous settings, use config defaults and **current directory** as working directory
  - `Ctrl+C` — cancel
- **`always`**: Automatically reuses without prompting
- **`never`**: Always uses config defaults and current directory as working directory

```bash
# Edit config file
nano ~/.config/condatainer/helper/vscode-server

# Set to never-reuse previous settings:
REUSE_MODE="never"
```

You can override specific reused settings with flags:

```bash
condatainer helper vscode-server -w   # reuse but switch to current directory
condatainer helper vscode-server -w -p 13182   # reuse but change port and directory
```

## Common Issues

### Port Already in Use

If you try to use a port that's already occupied:

```bash
# Let script auto-select a port
condatainer helper vscode-server

# Or specify a port explicitly
condatainer helper vscode-server -w -p 13182

# Or kill the process using the port (be careful with this!)
lsof -ti :13182 # list processes using the port
lsof -ti :13182 | xargs kill -9
```

### Cannot Connect

**Checklist:**
1. Application is running on server
2. SSH port forwarding is active: `ssh -L 13182:localhost:13182 user@server`
3. Correct URL in browser: `http://localhost:13182`
4. Check for firewall rules blocking the port

### Overlay In Use

If an overlay file is already in use by another process, the helper will:
1. Show you which processes are using the file
2. Ask if you want to kill them
3. Continue if you confirm, or abort otherwise

Manual check:
```bash
# Check what's using the overlay
lsof env.img

# Kill all processes using it
lsof -ti env.img | xargs kill -9
```

### Overlay Integrity Issues

Helper scripts automatically check and repair overlays. If manual repair needed:

```bash
e2fsck -f env.img
```

## Available Applications

### VS Code

**VS Code Server** - Full VS Code experience in browser

```bash
condatainer helper vscode-server -w -p 13182
```

Access at: `http://localhost:13182?tkn=<auto-generated-token>&folder=<current-dir>`

Options:
- `-p <port>` - Port number (default: auto-selected)
- `-a <token>` - Custom connection token (default: auto-generated)

**VS Code Tunnel** - No port forwarding needed (uses Microsoft's relay service)

```bash
condatainer helper vscode-tunnel -w
```

No `-p` flag needed - connects via Microsoft account authentication.

Options:
- `-n <name>` - Machine name for the tunnel (default: auto-generated)

**code-server** - Open-source VS Code alternative

```bash
condatainer helper code-server -w -p 13182
```

Options:
- `-p <port>` - Port number (default: auto-selected)
- `-a <auth>` - Password for authentication (default: none)

### RStudio

**RStudio Server (Docker-based)** - Using Posit R docker images

```bash
condatainer helper rstudio-server -w -r 4.4.3 -p 13182
```

Options:
- `-r <version>` - R version (e.g., 4.4.3, 4.4, or omit for latest)
- `-p <port>` - Port number (default: auto-selected)

**Note:** rstudio-server requires UID ≠ 1000 (reserved for Ubuntu default user in containers).

**RStudio Server (Conda-based)** - Using Conda R

```bash
condatainer helper rstudio-server-conda -w -p 13182
```

Must have R installed in your overlay first:
```bash
condatainer e -- mm-install r-base r-essentials rstudio
```

### Jupyter

**Jupyter Lab or Notebook**

```bash
# Jupyter Lab (default)
condatainer helper jupyter -w -p 13182

# Jupyter Notebook
condatainer helper jupyter -w -j notebook -p 13182
```

Options:
- `-j <lab|notebook>` - Run Jupyter Lab or Notebook (default: lab)
- `-p <port>` - Port number (default: auto-selected)
- `-a <token>` - Custom authentication token (default: auto-generated)

Access at: `http://localhost:13182?token=<auto-generated-token>`

### Desktop / GUI Applications

**XFCE4** - Full desktop environment via VNC/noVNC

```bash
condatainer helper xfce4 -w -p 13182
```

Options:
- `-p <port>` - Port for noVNC (default: auto-selected)
- `-a <passwd>` - VNC password (default: auto-generated)

**IGV** - XFCE desktop with IGV (Integrative Genomics Viewer) pre-installed

```bash
condatainer helper igv -w -p 13182
```

Options:
- `-p <port>` - Port for noVNC (default: auto-selected)
- `-a <passwd>` - VNC password (default: auto-generated)

## Differences from HPC Helpers

The headless server helpers are simpler than their HPC counterparts:

**What's the same:**
- Project-based overlay workflow
- Auto-install missing dependencies
- Configuration file management
- Overlay integrity checking

**What's different:**
- **No resource allocation flags** - no `-c`, `-m`, `-t`, `-g`, `-v` flags
- **Direct process management** - can kill processes using overlays on user confirmation

**When to use headless vs HPC helpers:**
- Use **headless helpers** when running on your own server or workstation
- Use **HPC helpers** when working on a cluster with SLURM/PBS/LSF/HTCondor
