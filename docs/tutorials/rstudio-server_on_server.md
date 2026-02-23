# Run RStudio Server on Headless Server

Simply run the following commands to set up and run RStudio Server on your remote server (without a job scheduler):

```bash
# download/update helper scripts
condatainer --local helper -u
# create a 30G overlay to store R packages and conda env
condatainer o -s 30g
# start RStudio Server
condatainer helper rstudio-server
```

This tutorial is for **headless servers** that do not have a job scheduler (SLURM, PBS, etc.). For HPC systems with a workload manager, see [RStudio Server on HPC](./rstudio-server_on_HPC.md) instead.

## Checklist

- [Have SSH port forwarding set up to the remote server](#ssh-port-forwarding)
- [Have CondaTainer installed](#install-condatainer)
- [Have a writable overlay image](#create-writable-overlay)
- [Check the Script Parameters](#rstudio-server-helper-script)

Then you can run:

```bash
condatainer helper rstudio-server
# Default
#   -r r_version=latest available (e.g., 4.5.2)
#   -p port_number=auto-selected if omitted
```

You can set [Config](./helpers_on_HPC.md#configuration-and-defaults) and helpers also support [Reusing previous settings](./helpers_on_HPC.md#reuse-mode).

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

After running this command, when your local machine accesses `localhost:<local_port>`, it will be forwarded to `localhost:<remote_port>` on the remote server.

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

Run the following command to install CondaTainer if it is not installed:

```bash
curl -fsSL https://get-condatainer.justype.net | bash
```

## Create Writable Overlay

Creating an ext3 overlay image for your packages and project files:

```bash
# create a 30G ext3 image named `env.img` in the current working directory
condatainer o -s 30g
```

If you need `reticulate`, you need to set up a conda environment inside the overlay:

```bash
# Launch a shell within the overlay
condatainer e
```

```bash
# You can change python version and add packages as needed
mm-install python=3.11 conda
```

See [Launch a Shell within the Workspace Overlay](../user_guide/workspace_overlays.md#launch-a-shell-within-the-workspace-overlay) for more details.

## RStudio Server Helper Script

The headless version works differently from the SLURM version:

- Runs the service directly on the current server (blocking)
- Checks if the port and overlay image are available
- Offers to kill processes using the overlay (with confirmation)
- Starts `rstudio-server` on the specified port

```
Usage: rstudio-server [options]

Options:
  -r <r_version>  R version to use (e.g., 4.4). If omitted, uses latest available.
  -p <port>       Port for rstudio-server (if omitted, an available port will be chosen)
  -e <overlay>    Environment overlay image file (default: env.img)
  -o <overlay>    Additional overlay files (can have multiple -o options)

  config          Show config file path and contents
```

### Available R Versions

The helper script uses [Posit R docker images](https://hub.docker.com/r/posit/r-base).

You can specify partial versions (e.g., `-r 4.4` will use the latest 4.4.x).

## Configuration

The script saves its defaults to `$XDG_CONFIG_HOME/condatainer/helper/rstudio-server` (defaults to `~/.config/condatainer/helper/rstudio-server`) on first run. Subsequent runs load from this file.

```bash
# View current config
condatainer helper rstudio-server config

# Reset to defaults (delete config file, next run recreates it)
rm ~/.config/condatainer/helper/rstudio-server
```

## Running RStudio Server

Let's set up and run `rstudio-server` on the headless server:

```bash
# Download the helper scripts (headless mode)
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
  rstudioapi::openProject("/home/your_username/current_working_directory")
```

Now you can go to your local web browser and access the link. Then run the provided R command to open the project.

Don't forget to set up SSH port forwarding from your local machine to the remote server before accessing the link.

```{note}
The script blocks while RStudio Server is running. Press `Ctrl+C` to stop the server.
```

```{tip}
You can run the script in the background using `nohup` or `tmux`/`screen`:
```

## Install Packages in CondaTainer RStudio Server

The `build-essential` overlay (auto-installed) contains common system libraries needed for building R packages. So you can directly install R packages from CRAN/Bioconductor/GitHub/source.

```R
install.packages("tidyverse")
BiocManager::install("DESeq2")
pak::pak("user/repo") # Or remotes::install_github("user/repo")
```

### Rprofile Setup

If `*.Rproj` file does not exist in the current working directory, **CondaTainer** will create:

- A `*.Rproj` file with default settings.
- A `.Rprofile` file to set up Posit Public Package Manager (P3M) CRAN and Bioconductor repositories.

See [Rprofile](https://github.com/Justype/condatainer/blob/main/helpers/slurm/.Rprofile) for more details. At the end, it will call `.set_repository_options()` for you. So you can directly download binary packages from P3M.

If `*.Rproj` file already exists, `.Rprofile` will not be created or modified. You need to set up the repositories yourself if needed.

So if you do not want to use Posit repositories, you can remove or edit the `.Rprofile` file accordingly.

If you only use CRAN packages, you can change the function in the `.Rprofile` to:

```R
.set_repository_options(repo = "cran", latest_cran = TRUE)
```

But if you plan to use Bioconductor packages, please keep the default settings to ensure compatibility.

### Install R Packages from Source

If a package is only available from source (e.g., GitHub), you need `remotes`, `pak` or other package managers. (Common system dependencies are included in `build-essential` overlay.)

I recommend using `pak`, it can tell you which system libraries are required.

```R
pak::pkg_sysreqs("user/repo@commit_hash") # or @tag
```

If system libraries are missing, you can create your `additional-deps` overlay with the required system libraries. See [Custom OS Overlays](../advanced_usage/custom_os.md) for more details.

After getting the `additional-deps.sqf` overlay, run `rstudio-server` with the `-o` option:

```bash
condatainer helper rstudio-server -o additional-deps.sqf
```

If you want to share your overlay with others, you should also provide the `def` file used to create it.

## Run without RStudio Server

You can run R scripts directly without RStudio Server:

```bash
# Basic execution
condatainer exec -o r4.4.3 -o build-essential -o env.img Rscript your_script.R
```

```bash
# If you created additional overlay `r-deps.sqf`:
condatainer exec -o r4.4.3 -o build-essential -o r-deps.sqf -o env.img Rscript your_script.R
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
