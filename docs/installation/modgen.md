# ModGen Installation

🏭 **ModGen** is a tool designed to automate the creation of **Lmod/Environment Modules** on HPC systems. It is ideal for users who prefer working with existing system-available modules rather than isolated containers.

But you cannot run rstudio-server, code-server using **ModGen**. For those use cases, please use [CondaTainer](condatainer.md).

## 🛠️ Quick Installation

**ModGen** is installed as part of the **CondaTainer** toolkit. To install it, run the following interactive script in your terminal:

```bash
curl -sL https://raw.githubusercontent.com/Justype/condatainer/main/assets/install.sh | bash -s -- -m
```

If you want to install both **CondaTainer** and **ModGen**, simply run the installation script with `-a`:

```bash
curl -sL https://raw.githubusercontent.com/Justype/condatainer/main/assets/install.sh | bash -s -- -a
```

## 📝 What to Expect During Installation

The installation script is **interactive** and will guide you through these steps:

1.  **Installation Path**: You will be prompted to confirm where the toolkit should be installed. The default locations are typically `$SCRATCH/condatainer/` or `$HOME/condatainer/`.
2.  **Shell Configuration**: The script will automatically **edit your shell configuration** (such as `.bashrc`) so that the toolkit commands are immediately accessible in your environment.

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

To remove **ModGen** and the associated toolkit, follow these steps:

If you want to remove both **CondaTainer** and **ModGen**, follow the [CondaTainer uninstallation instructions](condatainer.md#-uninstallation).

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
