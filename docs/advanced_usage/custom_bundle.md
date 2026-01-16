# Custom Bundle Overlays using Build Scripts

This guide focuses on building [bundle overlays](../user_guide/concepts.md#-overlay-types) using custom build scripts. The build scripts aim to let you modify the conda packages.

If you don't need modifications, directly use Conda YAML file:

```bash
condatainer create -p <prefix> -f <conda_yaml>
```

```{warning}
If you have packages not available via Conda, please directly build a os overlay instead. See the [Read-only R Package Environment](./custom_os.md#example-read-only-r-package-environment).
```

See [Build Scripts Manual](../manuals/build_script.md) for more details about writing build scripts.

## Step-by-Step Guide: Custom Build Script

In the following steps, we will:

1. Use conda to install `pySCENIC`.
2. Fix the [issue 475](https://github.com/aertslab/pySCENIC/issues/475) by editing `pySCENIC` files.

### Pause a Moment: Why Not Just Build an OS Overlay?

**CondaTainer** is good for managing writable environments during development. If you have already developed a pipeline and fixed package versions, you can directly build a read-only Apptainer (Singularity) image with all dependencies included.

If target app already has docker or singularity images, consider using them as base instead. See [Custom System Overlays](./custom_os.md#example-pulling-pytorch-docker-image) for details.

### 1. Create a Custom Build Script

`pyscenic.sh`

```bash
#!/bin/bash
# Available variables:
# - target_dir:      /cnt/name (set by condatainer)
# - tmp_dir:     tmp directory (set by condatainer)

micromamba create -y \
    --root-prefix $tmp_dir \
    --prefix $target_dir \
    -c conda-forge -c bioconda \
    --quiet \
    python=3.10 pyscenic=0.12.1 setuptools=79

# Fix issue 475:
sed -i 's/auc_thresholds\.iteritems()/auc_thresholds.items()/g' \
    $target_dir/lib/python3.10/site-packages/pyscenic/cli/utils.py
```

### 2. Create the Read-only Environment

```bash
condatainer create -f pyscenic.sh -p pyscenic
```

### 3. Use the Environment

```bash
condatainer exec -o pyscenic.sqf \
    pyscenic -h
```

### 4. Share the File

Share the `pyscenic.sh` file with your team or include it in your project repository. This allows others to build the same overlay and ensures consistency across different environments.

Or you can directly share the `pyscenic.sqf` overlay file for immediate use.
