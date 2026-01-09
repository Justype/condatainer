# How to run rstudio-server on remote server with CondaTainer

Before getting started, ensure you have [CondaTainer installed](../../README.md).

## Checklist

- [ ] [Have SSH port forwarding set up to the remote server](#ssh-port-forwarding)
- [ ] [Have CondaTainer installed](#install-condatainer)
- [ ] [Have required overlay images created](#install-required-overlays)
- [ ] [Have R installed in a writable overlay image](#create-r-writable-overlay)
- [ ] [Have the `rstudio-server-helper-local` in your PATH](#rstudio-server-helper-script)

Then you can run:

```bash
rstudio-server-helper-local \
  -p <port_number> \
  -e <overlay_image>

# Default 
#   port_number=8787
#   overlay_image=env.img
```

Please change the port number to a unique one to avoid conflicts with other users.

If you have preferred settings, you can modify the script directly. See [Change the default setting](#change-the-default-setting) section below.


## SSH Port Forwarding

```
Your machine -----> Remote Server
        (port forwarding)
```

`rstudio-server` is a web-based application, so you need to access it through a port on the remote server. 

```bash
ssh -L 8787:localhost:8787 your_username@remote_server_address
```

After running this command, when your local machine accesses `localhost:8787`, it will be forwarded to `localhost:8787` on the remote server.

> [!WARNING]
> Remote server is a shared system! 
> 
> Please do not use `8787` port. This is a common port and may be used by other users.
> 
> Instead, choose a strange and unique port number like `13182`.

You can also edit the SSH config file (`~/.ssh/config`) to add a section for your HPC system:

```
Host server
    HostName remote_server_address
    User your_username
    LocalForward 8787 localhost:8787
    # change 8787 to a unique port number
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
curl -sL https://raw.githubusercontent.com/Justype/condatainer/main/assets/install.sh | bash
```

## Install Required Overlays

Creating required overlay images:

```bash
condatainer install rstudio-server texlive build-essential
```

## Create R Writable Overlay

Creating an ext3 overlay image with R

```bash
# create 30G ext3 image named as env.img under current directory
condatainer overlay -s 30720

# go into the overlay
condatainer e

# INSIDE THE OVERLAY
# install R and required packages
mm-install r-base=4.4 r-tidyverse
# pin the R version to avoid accidental updates
mm-pin r-base
```

## RStudio Server Helper Script

You can use this script: [rstudio-server-helper-local](./rstudio-server-helper-local)

It will do the following steps for you:

- Check if the port is available
- Check the integrity of the writable overlay image
- Stop jobs using the overlay image
- Start `rstudio-server` using overlay images on that port

```
Usage: rstudio-server-helper-local [options]

Options:
  -p <port>       Port for rstudio-server (default: 8787)
  -e <overlay>    Environment overlay image file (default: env.img)
  -o <overlay>    Additional overlay files (can have multiple -o options)
  -y              Accept all warnings and proceed
  -h              Show this help message
```

Let's set up and run `rstudio-server` on HPC:

```bash
# Download the helper script
# Please make sure $HOME/bin is in your PATH
mkdir -p $HOME/bin
wget https://raw.githubusercontent.com/Justype/condatainer/main/docs/tutorials/rstudio-server-helper-local -O $HOME/bin/rstudio-server-helper-local
chmod +x $HOME/bin/rstudio-server-helper-local
```

Then you can run the script: 

```bash
rstudio-server-helper-local
```

After running the script, you will see output like this:

```
rstudio-server at http://localhost:8787
You can run the following command in R to open the project directly
  rstudioapi::openProject("/scratch/your_username/current_working_directory")
```

Then you can go to your local web browser and access `http://localhost:8787`.

And run the provided R command in R console to open project directly.

## Change the default setting

You can modify the default settings in the script, such as port number and overlay image path.

They are at the beginning of the script:

```bash
#!/bin/bash

PORT=8787
```

You can use editors like `vim` or `nano` to change these values.

or use `sed`

```bash
sed -i 's/PORT=8787/PORT=13182/' rstudio-server-helper
```
