# Run XFCE Desktop (VNC) on HPC

Run following commands to start an XFCE desktop session on your HPC system:

```bash
# download/update helper scripts
condatainer helper -u
# start XFCE desktop
condatainer helper xfce4
```

The desktop is accessible through a browser via noVNC. You can also use the `igv` helper to launch IGV (Integrative Genomics Viewer) inside the desktop.

## Checklist

- [A supported scheduler is available on your HPC system](#scheduler--workload-manager)
- [Have SSH port forwarding set up to the login node](#ssh-port-forwarding)
- [Have CondaTainer installed](#install-condatainer)
- [Check the Script Parameters](#helper-script-options)

Then you can run:

```bash
condatainer helper xfce4
# Default
#   -p port_number=auto-selected if omitted
#   -c num_cpus=3
#   -m memory=12G
#   -t time_limit=12:00:00
```

You can set [Config](./helpers_on_HPC.md#configuration-and-defaults) and helpers also support [Reusing previous settings](./helpers_on_HPC.md#reuse-mode).

If you have any issues, see [Common Issues](#common-issues) section below.

## xfce4 vs igv

| Helper | Description |
|--------|-------------|
| `xfce4` | General-purpose XFCE desktop session |
| `igv` | XFCE desktop with IGV auto-launched |

Both use the same VNC/noVNC setup. The `igv` helper additionally installs and launches IGV on startup.

```bash
# General desktop
condatainer helper xfce4

# Desktop with IGV
condatainer helper igv
```

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

noVNC is a web-based application, so you need to access it through a port on the HPC system.

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

## Helper Script Options

The script will do the following steps for you:

1. Check if a VNC desktop session is running on any compute node.
2. If yes, establish SSH port forwarding to that node.
3. If not,
   1. Check port and overlay integrity.
   2. Install `xfce4` overlay (and `igv` overlay for the igv helper) if missing.
   3. Set up VNC password and xstartup.
   4. Submit the SLURM job to start the desktop.
   5. When the job starts, set up SSH port forwarding.

```
Usage: xfce4 [options]

Options:
  -p <port>       Port for noVNC (default: randomly picked). Valid range: 1024-65535
  -a <passwd>     VNC password (if empty, one is generated)

  -c <number>     Number of CPUs to allocate (default: 3)
  -m <memory>     Amount of memory to allocate (default: 12G)
  -t <time>       Time limit for the job (default: 12:00:00)
  -g <gpu>        GPU resources (e.g., 1, a100:2). Use -g '' to clear
  -v              View Mode NCPUS:2 MEM:8G TIME:02:00:00
  -w              Use current directory as working directory

  -b <image>      Base image file
  -e <.img>       Writable overlay (default: env.img)
  -o <overlay>    Additional overlay files (can have multiple -o options)

  config          Show config file path and contents
```

## Configuration

The script saves its defaults to `$XDG_CONFIG_HOME/condatainer/helper/xfce4` (defaults to `~/.config/condatainer/helper/xfce4`) on first run. Subsequent runs load from this file.

```bash
# View current config
condatainer helper xfce4 config

# Reset to defaults (delete config file, next run recreates it)
rm ~/.config/condatainer/helper/xfce4
```

### Desktop Home Directory

By default, the desktop XDG directories (Desktop, Documents, Downloads, etc.) are redirected to `$SCRATCH/desktop-home/` to avoid cluttering your `$HOME` directory. You can change this by editing the `DESKTOP_HOME` setting in the config file.

## Running the Desktop

```bash
# Download the helper scripts
condatainer helper -u
```

Then you can run the script:

```bash
condatainer helper xfce4
```

After running the script, you will see output like this:

```
Waiting for job 1299123 to start running. Current status: PENDING.
Job 1299123 is now running on node cm23.
noVNC at http://localhost:<port>/vnc.html?autoconnect=1&resize=remote#password=<password>
If you want to stop it, run: scancel 1299123
```

Now you can click the link provided to open the desktop in your browser.

Don't forget to set up SSH port forwarding from your local machine to the HPC login node before accessing the link.

## Common Issues

### Port already in use

If your previously selected port is already in use, the script will notify you. You can either:
- Wait for the other process to finish
- Choose a different port with `-p <new_port>`

### Screen resolution

The default resolution is 1600x900. The noVNC URL includes `resize=remote` so the desktop will scale to your browser window.

## Using GPU

For GPU-accelerated applications, specify GPU resources:

```bash
condatainer helper xfce4 -g 1        # Request 1 GPU
condatainer helper xfce4 -g a100:2   # Request 2 A100 GPUs
```

You can use the following command to check available GPU types on your HPC system:

```
condatainer scheduler
```
