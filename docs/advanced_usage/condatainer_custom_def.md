# Custom Apptainer Definition Files

You can extend the functionality of **Condatainer** by creating custom Apptainer (Singularity) definition files. This allows you to build specific system-level overlays that can be loaded alongside your standard environment.

## File Placement Requirement

Currently, **Condatainer** only supports building definition files that are located within the build-scripts directory of your **Condatainer** installation.

Location: `condatainer_path/build-scripts/`

You must place your `.def` file in this specific folder for the build command to recognize it.

## Step-by-Step Guide: Additional Dependencies

Let's say when you building a R packages from source, you found that some system libraries are missing, e.g., `libxml2-dev`, `libcurl4-openssl-dev`, and `libssl-dev`.

```{tip}
You can use `pak::pkg_sysreqs("package_name_or_repo")` in R to check which system libraries are required for a specific package.
```

### 1. Create a Custom Definition File

You can create a custom Apptainer definition file to include these dependencies.

`additional-deps.def`

```
Bootstrap: docker
From: ubuntu:24.04

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

### 2. Place the Definition File

Move the `additional-deps.def` file into the `build-scripts` directory of your **Condatainer** installation.

```bash
mv additional-deps.def $SCRATCH/condatainer/build-scripts/
```

### 3. Build the Overlay

Use the `condatainer create` command to build the overlay from your custom definition file.

```bash
condatainer create additional-deps
```

This will generate an overlay file named `additional-deps.sqf` in your overlays directory.

### 4. Use the Custom Overlay

You can now load this custom overlay when running your analyses that require these additional system libraries.

```bash
condatainer exec -o build-essential -o additional-deps your_command_here
```

### 5. Share the Definition File

If you use the `additional-deps` overlay in your projects, consider sharing the `additional-deps.def` file with your team or including it in your project repository. This way, others can easily recreate the same environment by building the overlay from the provided definition file.
