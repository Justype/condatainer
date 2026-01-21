# Creating and Using Bundle Overlays

📦 **CondaTainer** allows you to create read-only, highly-compressed [bundle overlays](./concepts.md#-overlay-types) for your project.

They are ideal for sharing pre-configured environments with collaborators.

## Table of Contents

- [Create a Bundle Overlay](#create-a-bundle-overlay)
- [Launch a Shell within the Bundle Overlay](#launch-a-shell-within-the-bundle-overlay)
- [Use Bundle Overlays in a Script](#use-bundle-overlays-in-a-script)

## Create a Bundle Overlay

You can create a read-only SquashFS bundle overlay using the following command:

```bash
condatainer create -p prefix_name -f environment.yml
```

Make sure the conda env YAML file is properly defined.

```{warning}
CondaTainer uses the `sqf` name to locate resources inside the overlay. The overlay should not be renamed after creation. See [Naming Convention](./concepts.md#-naming-convention) for more details.
```

```{note}
If your environment uses packages not available in conda, like using `pip` or `remotes::install_github()`, check [Custom OS Overlays](../advanced_usage/custom_os.md) for more details.

But you can also [share the writable workspace overlay with others](./workspace_overlays.md#share-the-overlay-with-others).
```

```{note}
If you only need a slight modification of a conda env, like editing Python package code, see [Custom Bundle Overlays using Build Scripts](../advanced_usage/custom_bundle.md).
```

## Launch a Shell within the Bundle Overlay

To activate the bundle overlay, run the following command:

```bash
condatainer exec -o prefix_name.sqf bash
```

Then you can use the applications installed in the overlay.

Bundle overlays are read-only and stackable. You can mount multiple overlays together:

```bash
condatainer exec -o \
  grch38/salmon/1.10.2/gencode47 \
  prefix_name.sqf \
  bash
```

```{note}
Overlays mounted later will appear earlier in the `$PATH` (i.e., be prepended).
```

## Use Bundle Overlays in a Script

You can use the bundle overlay in your job scripts as follows:

For example, you have the following project structure:

```
project/
├── env.sqf
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
    -o env.sqf \
    python src/test.py
```

You can also write the `run_job.sh` in this way:

```bash
#!/bin/bash
#SBATCH --job-name=test_env
## Other SBATCH directives
#DEP:env.sqf

if [ -z "$IN_CONDATAINER" ] && command -v condatainer >/dev/null 2>&1; then
    condatainer run "$0" "$@"
    exit $?
fi

python src/test.py
```

You should run this from the project directory:

```bash
sbatch src/run_job.sh
```
