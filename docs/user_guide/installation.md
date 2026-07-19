# CondaTainer Installation

📦 **CondaTainer** is a rootless CLI for managing *tools* / *data* / *project environments* and launching apps on HPC. It can:

- Pack tools into highly-compressed, read-only SquashFS overlays. e.g. `cellranger/9.0.1`
- Generate writable environment overlay for development. e.g. `env.img`
- Run rstudio-server, code-server, and other web tools on HPC.

## 📋 Prerequisites

- **Linux (x86_64 only)**: AArch64 is not supported yet.
- **Apptainer/Singularity**: Required for all core container operations.
- **squashfs-tools**, **e2fsprogs**: For overlay creation and management.

Most HPC systems have already met these requirements. If not, please contact your system administrator to install them.

## 🛠️ Quick Installation

To install **CondaTainer**, run the command in your terminal:

```bash
curl -fsSL https://get-condatainer.justype.net | bash
```

The installation script will guide you through the following:
- **Installation Path**: You will be asked to confirm or set the path (defaults to `$SCRATCH/condatainer/` or `$HOME/condatainer/`).
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

The installation script will also run `condatainer config init` for you. If you want to re-run it or customize the configuration, run:

```bash
condatainer config init
```

**CondaTainer** will search for Apptainer/Singularity in your `PATH` and via `module avail`. If you want to use a specific Apptainer/Singularity module, please run:

```bash
module load apptainer/1.4.5 # or the version you want
condatainer config init
```

Compression is chosen automatically based on the runtime:
- **Singularity**: gzip
- **Apptainer >= 1.4**: zstd level 8
- **Apptainer < 1.4**: lz4

## ⌨️ Shell Completion

CondaTainer supports shell completion for **Bash** and **Zsh**. The installation script will automatically add the necessary lines to your shell configuration file:

- **Bash** (`~/.bashrc`):
- **Zsh** (`~/.zshrc`):

```bash
source <(condatainer completion)
```

## 🗺️ Next Steps

- [Concepts: Overlays](./concepts.md) — Understand the overlay model before proceeding
- [Helper Scripts](../tutorials/helpers.md) — Running Applications (VS Code, RStudio, IGV, etc.) on HPC or headless servers

## 🗑️ Uninstallation

To uninstall **CondaTainer**, follow these steps:

1.  **Clean Up Shell Configuration**: Open your `~/.bashrc` or `~/.zshrc` file and remove the lines related to **CondaTainer**. And find and remove the following lines:

```bash
# >>> CONDATAINER >>>
condatainer settings
# <<< CONDATAINER <<<
```

2. **Remove Data Directories**: CondaTainer stores images and scripts in one or more of the following locations. Check which exist and remove them in the following order:

| Type | Path | Descripton |
|----|----|----|
| Standalone | `<install-dir>/` | Install directory containing `bin/condatainer`|
| User Scratch | `$SCRATCH/condatainer/` | If `$SCRATCH` is set (most HPC systems) |
| XDG Data | `~/.local/share/condatainer/` | or `$XDG_DATA_HOME/condatainer/`  |

3. **Remove Config File**: The config file is stored at one of the following locations:

| Type | Path | Descripton |
|----|----|----|
| Standalone | `<install-dir>/config.yaml` | Same directory as `bin/` |
| XDG Config | `~/.config/condatainer/config.yaml` | or `$XDG_CONFIG_HOME/condatainer/` |

4. **Remove State Files**: Server and helper state files are stored at:

```
~/.local/state/condatainer/  (or $XDG_STATE_HOME/condatainer/)
```

5. **Remove Cache Files**: Available build scripts and installed OS overlays.

```
~/.cache/condatainer/        (or $XDG_CACHE_HOME/condatainer/)
```
