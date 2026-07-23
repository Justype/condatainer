# Creating and Using Bundle Overlays

📦 **CondaTainer** allows you to create read-only, highly-compressed [bundle overlays](./concepts.md#-overlay-types) for your project.

They are ideal for sharing pre-configured conda environments with collaborators.

## Table of Contents

- [Create a Bundle Overlay](#create-a-bundle-overlay)
- [Launch a Shell within the Bundle Overlay](#launch-a-shell-within-the-bundle-overlay)
- [Use Bundle Overlays in a Script](#use-bundle-overlays-in-a-script)
- [Export the Environment Spec](#export-the-environment-spec)

## Create a Bundle Overlay

You can create a read-only SquashFS bundle overlay using the following command:

```bash
condatainer create -p prefix_name -f environment.yml
```

Make sure the conda env YAML file is properly defined.

```{warning}
CondaTainer uses the name to locate resources inside the overlay. The bundle overlay should not be renamed after creation. See [Naming Convention](./concepts.md#-naming-convention) for more details.
```

```{note}
If you only need a slight modification of a conda env, like editing Python package code, see [Custom Bundle Overlays using Build Scripts](../advanced_usage/custom_bundle.md).
```

## Launch a Shell within the Bundle Overlay

To use the bundle overlay, run the following command:

```bash
$ condatainer exec -o prefix_name.sqf bash
# or condatainer e prefix_name.sqf
```

Then you can use the applications installed in the overlay.

```
CNT> <command>
```

Bundle overlays are read-only and stackable. You can mount multiple overlays together:

```bash
condatainer exec \
  -o grch38/salmon/1.10.2/gencode47 \
  -o prefix_name.sqf
```

```{note}
Overlays mounted later will appear earlier in the `$PATH` (i.e., be prepended).
```

## Use Bundle Overlays in a Script

You can use the bundle overlay in your job scripts as follows:

For example, you have the following project structure:

```
project/
├── custom.sqf
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
    -o custom.sqf \
    python src/test.py
```

You can also write the `run_job.sh` in this way:

```bash
#!/bin/bash
#SBATCH --job-name=test_env
## Other SBATCH directives
#DEP:custom.sqf

python src/test.py
```

You should run this from the project directory:

```bash
condatainer run src/run_job.sh
```

## Export the Environment Spec

A bundle overlay is built from a conda environment file, so `condatainer export` recovers the `environment.yml` that describes it.

```bash
condatainer export prefix_name.sqf > environment.yml
```

For an exact, fully-pinned spec that rebuilds without re-solving, use `-e`/`--explicit`. It round-trips through `create -f`:

```bash
# will export exact 
condatainer export prefix_name.sqf -e > spec.txt
condatainer create -f spec.txt -p rebuilt
```

See [`condatainer export`](../manuals/condatainer.md#export) for all options (`--from-history`, `--no-build`, channel handling, and more).
