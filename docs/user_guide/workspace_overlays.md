# Creating and Using Workspace Overlays

📦 **CondaTainer** allows you to create [workspace overlays](./concepts.md#-overlay-types) for your project.

The main advantage of **CondaTainer** is that it packs a Conda environment inside a writable ext3 overlay for development, which significantly reduces inode usage.

## Table of Contents

- [Create a Workspace Overlay](#create-a-workspace-overlay)
- [Launch a Shell within the Workspace Overlay](#launch-a-shell-within-the-workspace-overlay)
- [Writable or Read-Only](#writable-or-read-only)
- [Use Workspace Overlay in a Script](#use-workspace-overlay-in-a-script)
- [Set Environment Variables for the Overlay](#set-environment-variables-for-the-overlay)
- [Export the Conda Environment](#export-the-conda-environment)
- [Share the Overlay with Others](#share-the-overlay-with-others)
- [Common Issues](#common-issues)

## Create a Workspace Overlay

By default, CondaTainer will create a 10GiB ext3 overlay named `env.img` in the current working directory.

```bash
condatainer o
```

Initialize with specific packages directly (no YAML file needed):

```bash
condatainer o -- python=3.11 numpy
condatainer o myenv.img -- python=3.11 r-base

# Add extra channels with -c; conda-forge is always appended if not listed
condatainer o -- pytorch torchvision torchaudio pytorch-cuda=12.4 -c pytorch -c nvidia
```

Full command with options:

```
Usage:
  condatainer o [flags] [image_path] [-- packages...]

Examples:
  condatainer o # 10G with default inode ratio
  condatainer o my_data.img -s 50G -t data
  condatainer o --fakeroot --sparse
  condatainer o -f environment.yml
  condatainer o myenv.img -- python=3.11

Flags:
      --fakeroot      Create a fakeroot-compatible overlay (owned by root)
  -f, --file string   Initialize with Conda environment file (.yml or .yaml)
      --no-tmp        Create directly at target path instead of local tmp (slower on network filesystems)
  -s, --size string   Set overlay size (e.g., 500M, 10G) (default "10G")
  -S, --sparse        Create a sparse overlay image (no pre-allocation)
  -t, --type string   Overlay profile: small/balanced/large files (default "balanced")
```

Default channels when no `-c` is given: `conda-forge` + `bioconda`. When `-c` is given, `conda-forge` is appended automatically if not already listed. The resolved channels are saved in `.condarc` inside the overlay and reused by `mm-install`, `mm-update`, and `mm-search`.

- For Python projects, 10GiB is usually sufficient.
- For R projects, you may want to increase the size to 20GiB or more, especially if you are working with bioconductor packages.

```{note}
By default, CondaTainer stages the overlay at a fast local tmp directory (e.g. `/tmp`) and moves it to the target path once ready.
This avoids slow random I/O on HPC network filesystems (LustreFS) during image creation and conda initialization.

Use `--no-tmp` to skip this and write directly to the target path, if the tmp size is limited.
```

```{note}
Only one workspace overlay can be mounted at a time, regardless of writable or read-only mode.

The mounting path is `/ext3/env`. So the img name does not matter.
```

### Change the Overlay Size

Don't worry if you run out of space in the overlay later. You can easily resize it with the following command:

```bash
condatainer overlay resize env.img -s 30G
```

Or you can shrink it as well:

```bash
condatainer overlay resize env.img -s 5g
```

## Launch a Shell within the Workspace Overlay

To activate the workspace overlay, simply run the following command under the directory where the overlay image is created:

```bash
condatainer e
```

```
Usage:
  condatainer e [flags] [overlays...] [-- command...]

Flags:
  -b, --base-image string   Base image to use instead of default
      --bind strings        Bind path 'HOST:CONTAINER' (can be used multiple times)
      --env strings         Set environment variable 'KEY=VALUE' (can be used multiple times)
  -f, --fakeroot            Run container with fakeroot privileges
  -h, --help                help for e
  -n, --no-autoload         Disable auto-loading env.img from current directory
  -r, --read-only           Mount .img overlays as read-only (default: writable)
```

This command will mount the overlay and set the `PATH` and `CONDA_PREFIX` variables accordingly.

Then you can use `mm-*` commands to manage your project environment.

```bash
mm-install r-base=4.4 r-tidyverse  # Install packages
mm-install pytorch-cuda=12.4 -c pytorch -c nvidia  # Install with extra channels
mm-pin r-base           # Pin a package version
mm-pin -r r-base        # Unpin a package
mm-list                 # List installed packages
mm-search r-ggplot2     # Search for a package (uses saved channels)
mm-remove r-tidyverse   # Remove a package
mm-update               # Update packages
mm-clean -a             # Clean the cache and unused
mm-export               # Export environment

mm-channels get          # Show configured channels
mm-channels prepend pytorch  # Move/add channel to highest priority
mm-channels append bioconda  # Move/add channel to lowest priority
mm-channels remove nvidia    # Remove a channel
```

You don't need to activate the environment because the **CondaTainer** sets the `CONDA_PREFIX` and `PATH` for you.

## Writable or Read-Only

Two commands can enter a container with a workspace overlay, and they have different defaults:

| Command | Default mode | Auto-loads `env.img` | Best for |
|---------|-------------|----------------------|----------|
| `condatainer e` | **Writable** | Yes | Interactive development |
| `condatainer exec -o env.img` | **Read-only** | No (must specify `-o`) | Scripts and parallel jobs |

**Writable mode** (`condatainer e`) lets you install Conda packages and modify files inside the overlay, but it exclusively locks the image — only one process can mount it writable at a time.

**Read-only mode** allows multiple processes to mount the same overlay simultaneously, which is required for parallel job execution. Use `condatainer exec` with an explicit overlay path:

```bash
condatainer exec -o env.img <command>
```

To make `.img` overlays writable with `exec`, add `-w`/`--writable`:

```bash
condatainer exec -w -o env.img bash
```

## Set Environment Variables for the Overlay

You may need to set some environment variables when using the overlay in scripts or jobs.

For example, you may need to set `GOROOT` and `GOPATH` for Go projects.

Create `env.img.env` and add the following lines:

```
GOROOT=/ext3/env/go
GOPATH=/ext3/home/go
```

Then when you run `e` or `exec`, these environment variables will be set automatically.

```bash
$ condatainer e
[CNT][NOTE] Autoload env.img at /path/condatainer/env.img
[CNT] Overlay envs:
  GOROOT: /ext3/env/go
  GOPATH: /ext3/home/go
  CNT_CONDA_PREFIX: /ext3/env
```

You can add note to the `env.img.env` file as well:

```
GOROOT=/ext3/env/go
#ENVNOTE:GOROOT=Go Installation Path
GOPATH=/ext3/home/go
#ENVNOTE:GOPATH=Go Workspace Path
```

Then when you run `e` or `exec`, the notes will be displayed:

```bash
$ condatainer e
[CNT][NOTE] Autoload env.img at /path/condatainer/env.img
[CNT] Overlay envs:
  GOROOT: Go Installation Path
  GOPATH: Go Workspace Path
  CNT_CONDA_PREFIX: /ext3/env
```

## Use Workspace Overlay in a Script

You can use the workspace overlay in your job scripts as follows:

For example, you have the following project structure:

```
project/
├── env.img
└── src/
    ├── run_job.sh
    └── test.py
```

`run_job.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=test_env
## Other SBATCH directives

condatainer exec \
    -o env.img \
    python src/test.py
```

Then run `sbatch src/run_job.sh` to submit the job.

But I recommend using `condatainer run` to let **CondaTainer** check the overlay and the file locking for you:

```bash
#!/bin/bash
#SBATCH --job-name=test_env
## Other SBATCH directives
#DEP:env.img

python src/test.py
```

You should run this from the project directory:

```bash
# under project/ directory
condatainer run src/run_job.sh
```

The `run` will check the `img` file lock, before submission:
- if read-only (default), it will make sure no exclusive lock exists.
- if writable (`-w`), it will make sure no other lock exists.

## Export the Conda Environment

To export a reproducible environment spec from a workspace overlay (without entering it interactively), use:

```bash
condatainer overlay export env.img > environment.yml
```

Or mount the overlay and run `mm-export` directly:

```bash
condatainer e -- mm-export --no-builds > environment.yml
```

## Share the Overlay with Others

You can share the `env.img` overlay image with your collaborators or other researchers who want to reproduce your analysis environment.

Please compress the `env.img` file to reduce the file size before sharing.

```bash
gzip -k env.img
```

Or use zstd for better compression and speed:

```bash
zstd env.img
```

Because `env.img` is a standard ext3 filesystem image and preserves ownership information, your collaborators need to take care of the permissions when they use it.

### Permissions inside the overlay

By default, the overlay image is created with the user ID (UID) and group ID (GID) of the user who created it.

When another user mounts the overlay, the files inside will retain their original UID and GID, which may cause permission issues.

To avoid this, when other users get the overlay, they need to:

```bash
condatainer overlay chown env.img
```

Or they can use `--fakeroot` when running `condatainer exec` to avoid permission issues:

```bash
# Due to clusters apptainer configuration, this may not work.
condatainer exec --fakeroot -o env.img <command>
```

See [Fakeroot](../qa/overlayfs.md#fakeroot) for details on when and how fakeroot works.

## Common Issues

### Read-only file system

If you see "Read-only file system" errors when trying to write to the overlay, it means the overlay is used in read-only mode.

You need to exit the current shell and run it in writable mode:

```bash
condatainer e <name>.img
# or condatainer exec -w -o <name>.img <command>
```

### No space left on device

It means the overlay image is full. You can try to:

- run `mm-clean -a` to clean up unused packages and caches.
- increase the size of the overlay image by running:

```bash
condatainer overlay resize env.img -s 30G
```

### Permission errors when executing commands

See the [Permissions inside the overlay](#permissions-inside-the-overlay) section above, or the [Exec Troubleshooting](../qa/exec.md) FAQ for detailed steps.

```bash
# change UID/GID inside to your own
condatainer overlay chown env.img
```

### The overlay is used by another process

Cancel the job holding the lock (check `squeue`, `qstat`, `bjobs`, or `condor_q` depending on your scheduler), then run a filesystem check before remounting:

```bash
e2fsck -p env.img
```

See [Exec Troubleshooting](../qa/exec.md) for full instructions including per-scheduler commands and local server steps.
