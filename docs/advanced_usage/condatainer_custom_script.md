# Custom CondaTainer Build Scripts

**CondaTainer** supports two types of overlays:

- **Default Overlays**: Contain only the `/ext3/env` directory to store Conda environments or `/cnt/name/version` to store apps or data.
- **System Overlays**: Standalone containers featuring standard system directories (`/bin`, `/lib`, `/root`).

This guide focuses on building **default overlays** using custom build scripts. The build scripts aim to let you modify the conda packages or download data and pre-compiled binaries.

See [Build Script Manual](../manuals/build_script.md) for detailed instructions on creating custom build scripts.

```{warning}
If you have packages not available via Conda and want the environment read-only, please directly build a system overlay instead. See the [Read Only R Package Environment](./condatainer_custom_def.md#example-read-only-r-package-environment).
```

## Step-by-Step Guide: Custom Build Script

In the following steps, we will:

1. Use conda to install `pySCENIC`.
2. Fix the [issue 475](https://github.com/aertslab/pySCENIC/issues/475) by editing `pySCENIC` files.

### Pasue a Moment: Why Not Just Build a System Overlay?

**CondaTainer** is good for managing writable environments during development. If you have already developed a pipeline and fixed package versions, you can directly build a read-only Apptainer (Singularity) image with all dependencies included.

If target app already has docker or singularity images, consider using them as system overlays instead. See [Custom System Overlays](./condatainer_custom_def.md) for details.

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

### 2. Create the Read-Only Environment

```bash
condatainer create -f pyscenic.sh -p pyscenic
```

### 3. Use the Environment

```bash
condatainer exec -o pyscenic.sqf \
    pyscenic -h
```

### 4. Share the Build Script or Overlay

Consider sharing the `pyscenic.sh` file with your team or including it in your project repository. This way, others can easily recreate the same environment by building the overlay from the provided script.
