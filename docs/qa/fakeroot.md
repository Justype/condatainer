# Fakeroot

The `--fakeroot` flag allows an unprivileged user to appear as the **root user (UID 0)** inside the container. This relies on user namespace mapping (specifically `/etc/subuid` and `/etc/subgid` mappings on the host).

Fakeroot is useful when you:

- [Build the image](#build-the-image-with-fakeroot)
- [Run the container/overlay](#run-the-containeroverlay-with-fakeroot)

```{warning}
Fakeroot makes you appear as root, but it **cannot make read-only filesystems writable**. See [SIF vs SQF](#sif-vs-sqf).
```

## Build the image with fakeroot

On HPC systems, you don't have root privileges. Using `--fakeroot` allows you to create an overlay image that behaves like root inside the container.

By default, CondaTainer uses `--fakeroot` when creating overlays using `.def` files or pulling from remote sources.

Also, when running `exec`, `e`, `run`, if the uid of `upper` inside the overlay is `0` (root), CondaTainer will automatically add the `--fakeroot` flag.

## Run the container/overlay with fakeroot

When executing commands inside the container, you can use the `-f` or `--fakeroot` flag to gain root-like privileges inside the container.

```bash
condatainer exec --fakeroot -o /path/to/my_env.img <command>

# Using the 'e' shortcut
condatainer e -f my_env.img
```

This is useful when you need modification on existing images/overlays (like installing missing system libraries using `apt`)

## SIF vs SQF

**Fakeroot makes you appear as root, but it cannot make read-only filesystems writable.**

In CondaTainer, the container filesystem is assembled in layers:

| Layer | Format | Writable? | Paths |
|-------|--------|-----------|-------|
| Apptainer image | `.sif` (SIF) | No (SquashFS inside) | `/usr`, `/bin`, `/lib`, `/var`, … |
| Module/Bundle overlays | `.sqf` (SquashFS) | No | `/cnt/<name>/<version>` |
| Workspace overlay | `.img` (ext3) | Yes (when `-w`) | `/ext3/env` |

`apt` will modify `/usr/var` directories, but `sqf` overlays are read-only, so you cannot modify them even with `--fakeroot`.

The common ways to install missing system libraries:

- [Custom OS Overlays](../advanced_usage/custom_os.md) to create a new `.sqf` with missing libraries.
- Only mount the workspace overlay with `fakeroot` and install missing libraries.

In the example below, we try the second approach:

```bash
condatainer o  # create a workspace overlay named `env.img`
condatainer e r4.4.3 build-essential -- R
```

We created an overlay first, but found out that `tesseract-ocr-eng` is missing (which is not part of the `build-essential.sqf`). We can chown to root and install it:

```bash
condatainer overlay chown --root env.img
condatainer e # load env.img with fakeroot (no other overlays)
```

Within the container, run `apt` to install missing libraries:

```bash
apt update && apt install -y tesseract-ocr-eng
```

Then you need to chown back to you.

```bash
condatainer overlay chown env.img # chown back to you
condatainer e r4.4.3 build-essential -- R
```

Now, you can install R packages that require the package:

```r
install.packages("orderanalyzer") # required tesseract-ocr-eng
```

```{warning}
I don't recommend this path. IO performance of workspace overlays is bad. Try [custom OS overlays](../advanced_usage/custom_os.md) instead.
```

## Changing Overlay Ownership

You can use `overlay chown` to change the `upper` and `work` ownership inside the overlayFS.

```bash
condatainer overlay chown --root /path/to/my_env.img
```

You can also change the ownership to back to you:

```bash
condatainer overlay chown /path/to/my_env.img
```

See [condatainer manuals](../manuals/condatainer.md#overlay-chown) for more details.
