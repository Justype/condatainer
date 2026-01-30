# CondaTainer Installation

📦 **CondaTainer** is designed to manage tools/data/env on HPC systems by wrapping Conda environments into efficient SquashFS files using Apptainer. It can:

- Pack tools into highly-compressed, read-only SquashFS overlays. e.g. `cellranger/9.0.1.sqf`
- Run rstudio-server, code-server, and other web tools on HPC.
- Generate project-wide writable images for development. e.g. `env.img`

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

If your system has Apptainer modules instead of system-wide binaries, load the module first:

```bash
ml apptainer
condatainer config init
```

Then the condatainer will save the detected apptainer path. So next time you start a new terminal, it will work without loading the module again.

If apptainer version >= 1.4, zstd compression will be used by default. Otherwise, lz4 compression will be used.

## ⌨️ Shell Completion

CondaTainer supports shell completion for **Bash** and **Zsh**.

Add one of the following lines to your shell configuration file:

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
