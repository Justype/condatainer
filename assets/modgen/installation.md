# ModGen Installation

🏭 **ModGen** is a tool designed to automate the creation of **Lmod/Environment Modules** on HPC systems. It is ideal for users who prefer working with existing system-available modules rather than isolated containers.

But you cannot manage project environment and run rstudio-server using **ModGen**. For those use cases, please use [CondaTainer](../../README.md).

## 🛠️ Quick Installation

**ModGen** is installed as part of the **CondaTainer** toolkit. To install it only, run the following interactive script in your terminal:

```bash
curl -fsSL https://raw.githubusercontent.com/Justype/condatainer/main/assets/install.modgen.sh | bash
```

## 📝 What to Expect During Installation

The installation script is **interactive** and will guide you through these steps:

1.  **Installation Path**: You will be prompted to confirm or set the installation path (defaults to `$SCRATCH/condatainer/` or `$HOME/condatainer/`).
2.  **Shell Configuration**: The script will automatically **edit your shell configuration** (such as `.bashrc`) so that the toolkit commands are accessible.

## ✅ Verify Installation

Once the installation is complete, restart your terminal or run `source ~/.bashrc`.

```bash
source ~/.bashrc
```

Then, you can verify it is working by listing available recipes:

```bash
modgen avail
```

## 🗑️ Uninstallation

If you want to remove both **CondaTainer** and **ModGen**, follow the [CondaTainer uninstallation instructions](../../docs/installation/condatainer.md#️-uninstallation) to remove **CondaTainer** and then proceed to remove **ModGen** as described below.

If you just want to remove **ModGen**, open your `~/.bashrc` or `~/.zshrc` file and remove the lines related to **ModGen**. Look for and delete the following block:

```bash
# >>> MODGEN MODULES >>>
modgen configs
# <<< MODGEN MODULES <<<
```

And remove the **ModGen** executable and the modules directory.

```bash
INSTALL_DIR="$SCRATCH/condatainer"  # or your custom installation path
rm -rf "$INSTALL_DIR/bin/modgen"
rm -rf "$INSTALL_DIR/apps"
rm -rf "$INSTALL_DIR/apps-modules"
rm -rf "$INSTALL_DIR/ref"
rm -rf "$INSTALL_DIR/ref-modules"
```
