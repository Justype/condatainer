# OverlayFS and Fakeroot

## OverlayFS

Apptainer uses OverlayFS to stack multiple filesystems on top of the base `.sif` image. Each additional layer is passed via `--overlay <path>[:ro|:rw]`.

OverlayFS has three components:

| Component | Where | Role |
|-----------|-------|------|
| `lower` | the `.sif` itself  | Read-only base — never modified |
| `upper` | directory inside `.img` | Writable layer — all modifications land here (copy-on-write) |
| `work` | directory inside `.img` | Kernel-internal scratch directory |

A writable ext3 `.img` overlay stores both `upper/` and `work/` internally. When you write inside the container, the file is copied from `lower` into `upper` — the original layers are never modified.

- In `ro` (read-only) mode, all layers go to `lower` and no modifications are allowed.
- In `rw` (read-write) mode, the container writes to `upper` and modifications are persisted.

## Fakeroot

The `--fakeroot` flag allows an unprivileged user to appear as the **root user (UID 0)** inside the container. This relies on user namespace mapping (specifically `/etc/subuid` and `/etc/subgid` mappings on the host).

```{note}
Fakeroot makes you appear as root, but it **cannot make read-only filesystems writable**. A writable `.img` overlay is required for any writes inside the container.
```

### Build the image with fakeroot

On HPC systems, you don't have root privileges. Using `--fakeroot` allows you to create an overlay image that behaves like root inside the container.

By default, CondaTainer uses `--fakeroot` when creating overlays using `.def` files or pulling from remote sources.

Also, when running `exec`, `e`, `run`, if the uid of `upper` inside the overlay is `0` (root), CondaTainer will automatically add the `--fakeroot` flag.

### Run the container/overlay with fakeroot

When executing commands inside the container, you can use the `-f` or `--fakeroot` flag to gain root-like privileges inside the container.

```bash
condatainer exec --fakeroot -w -o /path/to/my_env.img <command>

# Using the 'e' shortcut
condatainer e -f my_env.img
```

This is useful when you need to modify existing images/overlays (like installing missing system libraries using `apt`).

## How to Use

### Installing missing system libraries

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
# load env.img with fakeroot (autodetects the upper ownership)
condatainer e r4.4.3 build-essential

# If apt install failed, try only loading the workspace overlay without sqfs:
condatainer e
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

### Changing Overlay Ownership

You can use `overlay chown` to change the `upper` and `work` ownership inside the overlayFS.

```bash
condatainer overlay chown --root /path/to/my_env.img
```

You can also change the ownership back to you:

```bash
condatainer overlay chown /path/to/my_env.img
```

See [condatainer manuals](../manuals/condatainer.md#overlay-chown) for more details.
