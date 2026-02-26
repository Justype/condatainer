# CondaTainer Installation

📦 **CondaTainer** is designed to manage tools/data/env on HPC systems by wrapping Conda environments into efficient SquashFS files or ext3 OverlayFS. It can:

- Pack tools into highly-compressed, read-only SquashFS overlays. e.g. `cellranger/9.0.1`
- Run rstudio-server, code-server, and other web tools on HPC.
- Generate project-wide writable images for development. e.g. `env.img`

## 📋 Prerequisites

- **Linux (x86_64 only)**: AArch64 is not supported yet.
- **Apptainer/Singularity**: Required for all core container operations.
- **squashfs-tools**, **e2fsprogs**: For overlay creation and management.

## 🛠️ Quick Installation

To install **CondaTainer**, run the following interactive script in your terminal:

```bash
curl -fsSL https://get-condatainer.justype.net | bash
```

The installation script is **interactive** and will guide you through the following:
- **Installation Path**: You will be asked to confirm or set the installation path (defaults to `$SCRATCH/condatainer/` or `$HOME/condatainer/`).
- **Shell Configuration**: The script will automatically **edit your shell configuration** (e.g., `.bashrc`) so that the `condatainer` command is ready to use.

## ✅ Verify Installation

Once the installation is complete, restart your terminal or run `source ~/.bashrc`.

```bash
source ~/.bashrc
```

Then, you can verify it is working by listing available recipes:

```bash
condatainer avail
```

## ⚙️ Initialize Configuration

Run the following command to create a default configuration file:

```bash
condatainer config init
```

If your system has Apptainer or Singularity modules instead of system-wide binaries, load the module first:

```bash
module load apptainer  # or: module load singularity
condatainer config init
```

This step ensures **CondaTainer** locates and saves the Apptainer path for future use, so future sessions work without reloading the module.

Compression is chosen automatically based on the runtime:
- **Singularity**: gzip (native default)
- **Apptainer >= 1.4**: zstd
- **Apptainer < 1.4**: lz4

## ⌨️ Shell Completion

CondaTainer supports shell completion for **Bash** and **Zsh**. Add the following lines to your shell configuration file:

- **Bash** (`~/.bashrc`): (must have `bash-completion` installed)
- **Zsh** (`~/.zshrc`):

```bash
source <(condatainer completion)
```

## 🗑️ Uninstallation

If you need to remove **CondaTainer** from your system, you can do so by reversing the steps taken during the interactive installation.

1.  **Remove the Installation Directory**: Delete the folder where **CondaTainer** was installed. eg. `rm -r $SCRATCH/condatainer/`
2.  **Clean Up Shell Configuration**: Open your `~/.bashrc` or `~/.zshrc` file and remove the lines related to **CondaTainer**. And find and remove the following lines:

```bash
# >>> CONDATAINER >>>
condatainer settings
# <<< CONDATAINER <<<
```

## Next Steps

- [Concepts: Overlays](./concepts.md) — Understand the overlay model before proceeding
- [Helpers on HPC](../tutorials/helpers_on_HPC.md) — Running Applications (e.g. VS Code, RStudio, IGV) with CondaTainer on HPC
- [Helpers on headless servers](../tutorials/helpers_on_server.md) — Running Applications without a scheduler
