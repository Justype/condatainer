# CondaTainer Installation

📦 **CondaTainer** is designed to manage tools/data/env on HPC systems by wrapping Conda environments into efficient SquashFS files using Apptainer. It can:

- Pack tools into highly-compressed, read-only SquashFS overlays. e.g. `cellranger/9.0.1.sqf`
- Run rstudio-server, code-server, and other web tools on HPC.
- Generate project-wide writable images for development. e.g. `env.img`

## 🛠️ Quick Installation

To install **CondaTainer**, run the following interactive script in your terminal:

```bash
curl -fsSL https://raw.githubusercontent.com/Justype/condatainer/main/assets/install.py-condatainer.sh | bash
```

## 📝 What to Expect

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

## ⌨️ Shell Completion

To make the tool even more user-friendly, the **shell completion** feature is automatically set up during installation.

By default, the installation script adds the following line to your shell configuration file (e.g., `~/.bashrc`):

```bash
source <(condatainer completion)
```

## 🗑️ Uninstallation

If you need to remove **CondaTainer** from your system, you can do so by reversing the steps taken during the interactive installation.

1.  **Remove the Installation Directory**: Delete the folder where **CondaTainer** was installed. eg. `rm -r $SCRATCH/condatainer/`
2.  **Clean Up Shell Configuration**: Open your `~/.bashrc` or `~/.zshrc` file and remove the lines related to **CondaTainer**. And find and remove the following lines:

```bash
# >>> CONDATAINER >>>
condatainer configs
# <<< CONDATAINER <<<
```
