# Use CondaTainer to Manage Read-Only Env

📦 **CondaTainer** allows you to create isolated read-only project environments using Apptainer and micromamba. Once the development and testing phase is complete, you can pack the project environment into a highly-compressed SquashFS overlay for production use. 

**CondaTainer** supports two types of overlays:

- **Default Overlays**: Use `/cnt/name/version` directory to store apps or data.
- **System Overlays**: Standalone containers featuring standard system directories (`/bin`, `/lib`, `/root`).

This guide focuses on creating **Default Overlays** using Conda YAML files.

## Table of Contents

- [Create a New Read-Only Environment](#create-a-new-read-only-environment)
- [Launch a Shell with the Project Environment](#launch-a-shell-with-the-project-environment)
- [Use `sqf` in a script](#use-sqf-in-a-script)

## Create a New Read-Only Environment

You can create a read-only SquashFS overlay using the following command:

```bash
condatainer create -p prefix_name -f environment.yml
```

Make sure the conda env YAML file is properly defined.

```{warning}
CondaTainer uses the `sqf` name to locate resources inside the overlay. The overlay should not be renamed after creation. See [Naming Convention](./concepts.md#-naming-convention) for more details.
```

```{note}
If your environment uses packages not available in conda, like using `pip` or `remotes::install_github()`, check [Custom System Overlays](../advanced_usage/condatainer_custom_def.md) for more details.

But you can also [share the writable overlay with others](./condatainer_project_level_writable.md#share-the-overlay-with-others).
```

```{note}
If you only need a slight modification of a conda env, like editing Python package code, see [Custom Build Script](../advanced_usage/condatainer_custom_script.md).
```

## Launch a Shell with the Project Environment

To activate the project environment, run the following command:

```bash
condatainer exec -o prefix_name.sqf bash
```

Then you can use the applications installed in the overlay.

An `sqf` overlay can only be mounted in read-only mode and can be stacked with other overlays.

```bash
condatainer exec -o \
  grch38/salmon/1.10.2/gencode47 \
  prefix_name.sqf \
  bash
```

```{note}
Overlays mounted later will appear earlier in the `$PATH` (i.e., be prepended).
```

## Use sqf in a Script

You can use the project environment overlay image in your job scripts as follows:

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
