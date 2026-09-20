# Custom Bundle Overlays

This guide focuses on building [bundle overlays](../user_guide/concepts.md#overlay-types) using custom build scripts. The scripts aim to let you modify the conda packages.

If you don't need modifications, directly use Conda YAML file:

```bash
condatainer create -p <prefix> -f <conda_yaml>
```

```{warning}
If you have packages not available via Conda, please directly build a os overlay instead. See the [Read-only R Package Environment](./custom_os.md#example-read-only-r-package-environment).

For a single app shipped as a vendor tarball or installer, write an [app recipe](./custom_app.md) instead.
```

See [Build Scripts Manual](../manuals/build_script.md) for more details about writing build scripts.

## Step-by-Step Guide: Custom Build Script

In the following steps, we will:

1. Use conda to install `pySCENIC`.
2. Fix the [issue 475](https://github.com/aertslab/pySCENIC/issues/475) by editing `pySCENIC` files.

### Pause a Moment: Why Not Build an OS Overlay or Apptainer Image?

**CondaTainer** is good for managing writable environments during development. But if you have already developed a pipeline and fixed package versions, you can also directly build a read-only Apptainer image or OS overlay with all dependencies included.

If target app already has docker or singularity images, consider using them as base instead. See [Custom System Overlays](./custom_os.md#example-pulling-pytorch-docker-image) for details.

### 1. Create a Custom Build Script

`pyscenic.sh`

```bash
#!/bin/bash
# Available variables:
# - CNT_PREFIX: /cnt/name, where the environment goes (set by condatainer)
# - CNT_TMP:    scratch directory (set by condatainer)

# $CNT_PREFIX already exists, so mark it as an environment before installing
mkdir -p "$CNT_PREFIX/conda-meta"
touch "$CNT_PREFIX/conda-meta/history"

micromamba install -y \
    --prefix "$CNT_PREFIX" \
    -c conda-forge -c bioconda \
    --quiet \
    python=3.10 pyscenic=0.12.1 setuptools=79

# Fix issue 475:
sed -i 's/auc_thresholds\.iteritems()/auc_thresholds.items()/g' \
    "$CNT_PREFIX/lib/python3.10/site-packages/pyscenic/cli/utils.py"
```

While a script builds, `micromamba` is on `$PATH` (CondaTainer installs it on first use), and `$MAMBA_ROOT_PREFIX` is the build's scratch and your `.condarc` is ignored (`MAMBA_NO_RC=true`), so package caches do not land in your home directory. Name the channels with `-c`, as above.

```{note}
The script is identified by its text, not by the packages micromamba resolves. Pin the versions you depend on, as above; the same script run later may resolve something newer.
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

Sharing the `pyscenic.sqf` overlay is self-contained: the build script is embedded inside it, so collaborators can use it immediately **and** recover the exact recipe with `condatainer export pyscenic.sqf`. There is no longer a need to hand over the `.sh` separately.

You can still share (or version-control) the `pyscenic.sh` script if you prefer others rebuild from source rather than copy the image.
