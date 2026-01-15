# Use CondaTainer to Manage Writable Env

📦 **CondaTainer** allows you to create isolated project environments using Apptainer and micromamba.

The main advantage of **CondaTainer** is that it packs the project environment inside a writable `ext3` overlay for development, which significantly reduces inode usage compared to traditional conda environments.

## Table of Contents

- [Create a New Project Environment](#create-a-new-project-environment)
- [Launch a Shell with the Project Environment](#launch-a-shell-with-the-project-environment)
- [Writable or Read-Only](#writable-or-read-only)
- [Use `env.img` in a script](#use-envimg-in-a-script)
- [Share the Overlay with Others](#share-the-overlay-with-others)
- [Common Issues](#common-issues)

## Create a New Project Environment

By default, CondaTainer will create a 20GiB ext3 overlay named `env.img` in the current working directory.

```bash
condatainer o
```

Full command with options:

```
usage: condatainer o [-h] [-s SIZE] [-f FILE] [--fakeroot] [--sparse] [image]

positional arguments:
  image                 Path to the output overlay image (default: env.img)

options:
  -h, --help            show this help message and exit
  -s SIZE, --size SIZE  Set overlay size (default: 10G). Accepts GB/MB suffixes; assumes MB if omitted
  -f FILE, --file FILE  Initialize with Conda environment file (.yaml/.yml)
  --fakeroot            Create a fakeroot-compatible overlay
  --sparse              Create a sparse overlay image (default: false)
```

- For Python projects, 10GiB is usually sufficient.
- For R projects, you may want to increase the size to 20GiB or more, especially if you are working with bioconductor packages.

```{note}
Only one overlay can be mounted at a time, regardless of writable or read-only mode.

The mounting path is `/ext3/env`. So the img name does not matter.
```

### Change the Overlay Size

Don't worry if you run out of space in the overlay later. You can easily resize it with the following command:

```bash
condatainer overlay resize env.img -s 30G
```

Or you can shrink it as well:

```bash
condatainer overlay resize env.img -s 5G
```

## Launch a Shell with the Project Environment

To activate the project environment, simply run the following command under the directory where the overlay image is created:

```bash
condatainer e
```

```
usage: condatainer e [-h] [-r] [-n] [--fakeroot] [overlays ...]

positional arguments:
  overlays           Overlay files to mount (env.img if not provided)

options:
  -h, --help         show this help message and exit
  -r, --read-only    Do not make .img overlays writable (default: writable)
  -n, --no-autoload  Do not autoload local env.img from current directory
  --fakeroot         Run command with fake root privileges inside the container
```

This command will mount the overlay and set the `PATH` and `CONDA_PREFIX` variables accordingly.

Then you can use `mm-*` commands to manage your project environment.

```bash
# Install numpy package
mm-install numpy

# Pin a package version
mm-pin numpy
# Unpin a package
mm-pin -d numpy

# List installed packages
mm-list

# Remove a package
mm-remove numpy

# Update
mm-update numpy

# Clean cache, tarballs, unused packages
mm-clean -a

# Export environment
mm-export > environment.yaml
```

## Writable or Read-Only

`condatainer e` mounts the overlay in writable mode by default.

The limitations of writable mode are:
- Only one writable overlay can be mounted at a time.
- If an overlay is mounted in writable mode, other processes cannot mount it in read-only mode.

So, if you need to use it to run parallel jobs, you can mount it in read-only mode:

```bash
condatainer exec -o env.img <command>
```

The `exec` command:

- Mounts the overlay in read-only mode.
- Does not automatically mount the overlay, so you need to specify the overlay image file with the `-o` option.

## Use `env.img` in a script

You can use the project environment overlay image in your job scripts as follows:

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

You can also write the `run_job.sh` in this way:

```bash
#!/bin/bash
#SBATCH --job-name=test_env
## Other SBATCH directives
#DEP:env.img

if [ -z "$IN_CONDATINER" ] && command -v condatainer >/dev/null 2>&1; then
    condatainer run "$0" "$@"
    exit $?
fi

python src/test.py
```

You should run this from the project directory:

```bash
sbatch src/run_job.sh
```

And make sure when you run the script, the `env.img` is not mounted in writable mode by another process. You can also make a copy of `env.img` and use that one instead.

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

Remember they need to run this command instead of you.

Or they can use `--fakeroot` when running `condatainer exec` to avoid permission issues:

```bash
# May not work on HPC systems
condatainer exec --fakeroot -o env.img <command>
```

## Common Issues

### No space left on device

It means the overlay image is full. You can try to:

- run `mm-clean -a` to clean up unused packages and caches. 
- increase the size of the overlay image by running:

```bash
condatainer overlay resize env.img -s 30G
```

### Permission errors when executing commands

See the [Permissions inside the overlay](#permissions-inside-the-overlay) section above.

```bash
# change UID/GID inside to your own
condatainer overlay chown env.img
```

### The overlay is used by another process

On the HPC system, you can use `squeue -u $USER` to check if there are other jobs running under your account that may be using the overlay in writable mode.

If you want to stop it immediately, you can use `scancel <job_id>` to cancel the job.

On the local machine, you can use `lsof env.img` to check if there are other processes using the overlay. Then you can use `kill <pid>` to stop the process. Or use `fuser -k env.img` to kill all processes using the overlay.

Make sure you run `e2fsck -p env.img` to check and fix the overlay image before mounting it again.

```bash
e2fsck -p env.img
```

If `-p` option does not work, you can run `e2fsck env.img` without `-p` to fix the image manually.

```bash
e2fsck env.img
```
