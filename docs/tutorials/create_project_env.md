# Create Project Environment

Like [pixi](https://pixi.prefix.dev/latest/) or other project environment managers, CondaTainer allows you to create isolated project environments using micromamba.

The main advantage of CondaTainer is it packs the project environment inside an overlay, which significantly reduces the inode usage on your HPC system. From 50k+ to a single ext3 file!

## Table of Contents

- [Prerequisites](#prerequisites)
- [Create a New Project Environment](#create-a-new-project-environment)
- [Launch a Shell with the Project Environment](#launch-a-shell-with-the-project-environment)
- [Writable or Read-Only](#writable-or-read-only)
- [Use `env.img` in A Script](#use-envimg-in-a-script)
- [Common Issues](#common-issues)

## Prerequisites

Have CondaTainer installed and set up on your HPC system.

## Create a New Project Environment

By default, CondaTainer will create a 10GiB ext3 overlay named as `env.img` under current working directory.

```bash
condatainer overlay
```

Full command with options:

```
usage: condatainer overlay [-h] [-s SIZE] [-f FILE] [image]

positional arguments:
  image                 Path to create the overlay image

options:
  -h, --help            show this help message and exit
  -s SIZE, --size SIZE  Size of the overlay image in MiB
  -f FILE, --file FILE  Conda env YAML file
```

- For python projects, 10GiB is usually sufficient.
- For R projects, you may want to increase the size to 20GiB or more especially if you are working with bioconductor packages.

## Launch a Shell with the Project Environment

To activate the project environment, simply run the following command under the directory where the overlay image is created:

```bash
condatainer e
```

This command will mount the overlay and set the `PATH` and `CONDA_PREFIX` variables accordingly.

Then you can use `mm-*` commands to manage your project environment.

```bash
# Install numpy package
mm-install numpy

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

`exec`

- Mounts the overlay in read-only mode.
- Does not auto mount the overlay, so you need to specify the overlay image file with `-o` option.

## Use `env.img` in A Script

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

SCRIPT_PATH=$(realpath "$0")
SCRIPT_DIR=$(dirname "$SCRIPT_PATH")
PROJECT_DIR=$(dirname "$SCRIPT_DIR")

condatainer exec \
    -o "$PROJECT_DIR/env.img" \
    python "$SCRIPT_DIR/test.py"
```

You can use absolute path if you want.

And make sure when you run the script, the `env.img` is not mounted in writable mode by another process. You can also make a copy of `env.img` and use that one instead.

## Common Issues

### No space left on device

It means the overlay image is full. You can try to:

- run `mm-clean -a` to clean up unused packages and caches. 
- increase the size of the overlay image by creating a new one with larger size and reinstall the packages.

### There is no other process using the overlay, but still cannot mount in read-only mode

You can:

1. Try to fix the overlay image `e2fsck -p env.img`
2. If that does not work, copy the overlay image to a new file and use that one.

```bash
e2fsck -p env.img
```

```bash
cp env.img env_copy.img

# And use env_copy.img instead
condatainer exec -o env_copy.img <command>
```
