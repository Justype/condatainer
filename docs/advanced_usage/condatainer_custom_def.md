# Custom System Overlays

You can extend the functionality of **Condatainer** by creating custom Apptainer (Singularity) definition files or by directly pulling from remote sources like Docker Hub. This allows you to build specific system-level overlays that can be loaded alongside your standard environment.

```bash
# Directly pull from a remote Docker image
condatainer create -p <prefix> -s docker://<docker_image>

# Build from a custom definition file
condatainer create -p <prefix> -f <path_to_def_file>
```

```{note}
Make sure you use the same distro version as the base image. Current default: `ubuntu:24.04`
```

## Change the Base Image

If a system overlay has a different distro version than the base image used by other overlays, you may run into compatibility issues when loading multiple overlays together.

To avoid this, you can create your own base image.

### Example: PyTorch image

For example, you might want to directly pull the PyTorch Docker image (it uses `ubuntu:20.04` as the base image):

When you use that overlay directly, you may encounter missing library errors in Python.

```bash
condatainer create -p pytorch -s docker://pytorch/pytorch:2.9.1-cuda13.0-cudnn9-devel

# This may produce an error when you run:
condatainer exec -o pytorch.sqf bash
```

In the inner shell, you can see the OS version with:

```bash
cat /etc/os-release
```

It is Ubuntu 22.04, which differs from the default base image (`ubuntu24`).

You can either:
- Build an Ubuntu 22 base image
- Use the PyTorch image as the base image

Use the second approach:

```bash
condatainer exec -b pytorch.sqf bash
```

Or create an Ubuntu 22 base image:

```bash
# Use available base image definition file
condatainer create ubuntu22/base_image

# Then load both the base image and the PyTorch overlay:
condatainer exec -b ubuntu22/base_image -o pytorch.sqf bash
```

```{note}
Available base image definition files:

- [ubuntu20/base_image](https://github.com/Justype/condatainer/tree/main/build-scripts/ubuntu20/base_image.def)
- [ubuntu22/base_image](https://github.com/Justype/condatainer/tree/main/build-scripts/ubuntu22/base_image.def)
- [ubuntu24/base_image](https://github.com/Justype/condatainer/tree/main/build-scripts/ubuntu24/base_image.def)
```

## Step-by-Step Guide: Additional Dependencies

Let's say when building an R package from source, you find that some system libraries are missing, e.g., `libxml2-dev`, `libcurl4-openssl-dev`, and `libssl-dev`.

```{tip}
You can use `pak::pkg_sysreqs()` in R to check which system libraries are required for specific packages.
```

### 1. Create a Custom Definition File

You can create a custom Apptainer definition file to include these dependencies.

`additional-deps.def`

```
Bootstrap: docker
From: ubuntu:24.04
# Make sure distro version matches the base image

%labels
    Version 1.0
    Description "Ubuntu24-based container with additional dependencies"

%environment
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8
    export PS1='CNT \w> '

%runscript
    exec /bin/bash "$@"	

%post
    # Avoid interactive prompts
    export DEBIAN_FRONTEND=noninteractive

    # Put your custom installation commands here
    apt -y update && apt -y install \
        libxml2-dev \
        libcurl4-openssl-dev \
        libssl-dev

    apt clean
```

### 2. Build the Overlay

Use the `condatainer create` command to build the overlay from your custom definition file.

```bash
condatainer create -p additional-deps -f additional-deps.def
```

This will generate an overlay file named `additional-deps.sqf` in the current directory.

### 3. Use the Custom Overlay

You can now load this custom overlay when running analyses that require these additional system libraries.

```bash
condatainer exec -o build-essential -o additional-deps.sqf your_command_here
```

### 4. Share the Definition File

If you use `additional-deps.sqf` in your projects, consider sharing the `additional-deps.def` file with your team or including it in your project repository.

This way, others can easily recreate the same environment by building the overlay from the provided definition file.
