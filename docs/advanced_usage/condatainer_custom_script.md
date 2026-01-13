# Custom CondaTainer Build Scripts

If you have packages not available via Conda, you can create a custom build scripts to make the installation process reproducible.

See [Build Script Manual](../manuals/build_script.md) for detailed instructions on creating and using custom build scripts with CondaTainer.

## Step-by-Step Guide: Custom Build Script

Let's say you want to install `DoubletFinder`, which is not available via Conda, but can be installed from GitHub using R. And you also want the environment read-only.

So you need:

1. Use conda to install R and other packages.
2. Use `pak` to install `DoubletFinder` from GitHub.

```{note}
R will compile packages from source. If you need additional system libraries, check this: [Custom Definition File](condatainer_custom_def.md#step-by-step-guide-additional-dependencies).
```

### 1. Create a Custom Build Script

`doubletfinder.sh`

```bash
#!/bin/bash
#DEP:build-essential
# Set the dependencies may need build-essential for R package compilation

# Available variables:
# - target_dir:      /cnt/name (set by condatainer)
# - tmp_dir:     tmp directory (set by condatainer)

micromamba create -y \
    --root-prefix $tmp_dir \
    --prefix $target_dir \
    -c conda-forge -c bioconda \
    --quiet \
    r-base=4.4 r-pak r-seurat=5.3.0 r-tidyverse=2.0.0

export PATH="$target_dir/bin:$PATH"

# Then run R script to install DoubletFinder
Rscript --vanilla - <<'EOT'
pak::pak("chris-mcginnis-ucsf/DoubletFinder")
EOT
```

### 2. Create the Read-Only Environment

```bash
condatainer create -f doubletfinder.sh -p doubletfinder-env
```

### 3. Use the Environment

```bash
condatainer exec -o doubletfinder-env.sqf \
    Rscript --vanilla -e "library(DoubletFinder); sessionInfo()"
```

### 4. Share the Build Script

Consider sharing the `doubletfinder.sh` build script with your team or including it in your project repository. This way, others can easily recreate the same environment by building the overlay from the provided script.

Also in the build script, you may need to specify the GitHub commit hash for reproducibility.

```R
pak::pak("chris-mcginnis-ucsf/DoubletFinder@1b244d8")
```
