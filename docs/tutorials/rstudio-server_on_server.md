# Run rstudio-server on server

## Checklist

- [Have SSH port forwarding set up to the remote server](#ssh-port-forwarding)
- [Have CondaTainer installed](#install-condatainer)
- [Have required overlay images created](#install-required-overlays)
- [Have R installed in a writable overlay image](#create-r-writable-overlay)
- [Check the Script Parameters](#rstudio-server-helper-script)

Then you can run:

```bash
condatainer helper \
  rstudio-server \
  -e <overlay_image>

# Default
#   port_number=auto-selected if omitted
#   overlay_image=env.img
```

Please change the port number to a unique one to avoid conflicts with other users.

```bash
# You can create alias in your shell config file (~/.bashrc or ~/.zshrc):
# Change 13182 to your preferred port number
alias rstudio-server-start='condatainer helper rstudio-server -p 13182'
```

Always [Use Conda to Manage R Packages](#use-conda-to-manage-r-packages). For packages only available from source, see [Install R Packages from Source](#install-r-packages-from-source) section.

If you encounter `file name too long` error, see [File name too long ERROR](#file-name-too-long-system36-error) section below.

## SSH Port Forwarding

```
Your machine -----> Remote Server
        (port forwarding)
```

`rstudio-server` is a web-based application, so you need to access it through a port on the remote server. 

```bash
ssh -L <local_port>:localhost:<remote_port> your_username@remote_server_address
```

After running this command, when your local machine accesses `localhost:<local_port>`, it will be forwarded to `localhost:<remote_port>` on the remote server. If you omit `-p` when launching the helper, the helper will auto-select a port and print the port number and an `ssh -L` example you can copy.

```{warning}
Remote server is a shared system! Please do not use a common port like `8787`.

Instead, choose an uncommon and unique port number like `13182`.
```

You can also edit the SSH config file (`~/.ssh/config`) to add a section for your remote server:

```
Host server
    HostName remote_server_address
    User your_username
    LocalForward <local_port> localhost:<remote_port>
    # Please change <local_port> to a unique port number
```

Then you can simply run:

```bash
ssh server
```

## Install CondaTainer

Run the following command to check if CondaTainer is installed:

```
condatainer --version
```

Run the following command to install CondaTainer if it is not installed:

```bash
curl -fsSL https://get-condatainer.justype.net/ | bash
```

Download the helper scripts:

```bash
condatainer --local helper --update
```

## Install Required Overlays

Creating required overlay images:

```bash
condatainer install rstudio-server
```

## Create R Writable Overlay

Creating an ext3 overlay image with R

```bash
# create a 30G ext3 image named `env.img` in the current working directory
condatainer overlay -s 30720

# go into the overlay
condatainer e
```

In side the overlay, you can use `mm-*` commands

```bash
# install R and required packages
mm-install r-base=4.4 r-tidyverse
```

See [Launch a Shell within the Workspace Overlay](../user_guide/workspace_overlays.md#launch-a-shell-within-the-workspace-overlay) for more details.

## RStudio Server Helper Script

It will do the following steps for you:

- Check if the port is available
- Check the integrity of the writable overlay image
- Stop jobs using the overlay image
- Start `rstudio-server` using overlay images on that port

```
Usage: rstudio-server [options]

Options:
  -p <port>       Port for rstudio-server (if omitted, an available port will be chosen)
  -b <image>      Base image file
  -e <overlay>    Environment overlay image file (default: env.img)
  -o <overlay>    Additional overlay files (can have multiple -o options)
  -y              Accept all warnings and proceed
  -h              Show this help message
```

Let's set up and run `rstudio-server` on remote headless server:

```bash
# Download the helper scripts
condatainer --local helper -u
```

Then you can run the script: 

```bash
condatainer helper rstudio-server
```

After running the script, you will see output like this:

```
rstudio-server at http://localhost:<port>
You can run the following command in R to open the project directly:
  rstudioapi::openProject("/scratch/your_username/current_working_directory")
```

Then you can go to your local web browser and access `http://localhost:<port>`. The helper prints the actual port it chose when started.

Run the provided R command in the R console to open the project directly.

## Use Conda to Manage R Packages

It is highly recommended to use Conda (via `mm-install`) to install R packages whenever possible. They are pre-compiled.

- CRAN packages: `mm-install r-<package_name>`
- Bioconductor packages: `mm-install bioconductor-<package_name>`

e.g.

```bash
# in the overlay using condatainer e

# find available R packages in conda-forge/bioconda
# Or go to https://anaconda.org
mm-search r-presto

# Install packages
mm-install r-base=4.4 r-ggplot2 bioconductor-deseq2

# Pin the R version to avoid accidental updates
mm-pin r-base
```

## Install R Packages from Source

```{note}
Before installing R packages from source, please make sure the package is not available from conda-forge or bioconda.
```

If a package is only available from source (e.g., GitHub), you need `remotes`, `pak` or other package managers.

I recommend using `pak`, it can tell which system libraries are required.

```R
mm-install r-pak
```

Then in R:

```R
pak::pkg_sysreqs("user/repo@commit_hash") # or @tag
```

If system libraries are missing, you can create your `additional-deps` overlay with the required system libraries. (ignore pandoc missing warning)

see [Custom OS Overlays](../advanced_usage/custom_os.md) for more details.

After getting the `additional-deps.sqf` overlay, run `rstudio-server` with the `-o` option:

```bash
condatainer helper rstudio-server -o additional-deps.sqf
```

If you want to share your overlay with others, you should also provide the `def` file used to create it.

## Run without RStudio Server

If you only use R packages from Conda, you can run R directly without RStudio Server:

```bash
condatainer exec -o env.img Rscript your_script.R
```

If you have R packages built from GitHub or source, you need to load your additional overlay too:

```bash
condatainer exec -o r-deps.sqf -o env.img Rscript your_script.R
```

## File name too long system:36 ERROR

If you encounter the following error when starting `rserver`:

```
TTY detected. Printing informational message about logging configuration. Logging configuration loaded from '/etc/rstudio/logging.conf'. Logging to '/home/cic/devgab/.local/share/rstudio/log/rserver.log'.
2024-07-31T16:26:50.011828Z [rserver] ERROR Unexpected exception: File name too long [system:36]; LOGGED FROM: int main(int, char* const*) src/cpp/server/ServerMain.cpp:975
```

See [this issue](https://github.com/rstudio/rstudio/issues/15024) for more details.

This is caused by `rserver` trying to create runtime files with long UUID names.

For example, if your home/scratch directory is in a deep path like `/workspace/lab/long_last_name_lab/<your_username>/`, then RStudio will try to create socket files with very long paths like: 

```
/workspace/lab/long_last_name_lab/<your_username>/.local/share/rstudio/run/other_stuff/<very_long_uuid>/rstudio-server/session-server-rpc.socket
```

To fix this, you need to change this line:

```
        --server-data-dir=$SCRATCH/.local/run/rstudio-server \
```

to

```
        --server-data-dir=$SCRATCH/.r \
```

Hopefully, this can fix the issue. 

If that does not work, please contact your system administrator for assistance.
