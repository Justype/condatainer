# Fakeroot

The `--fakeroot` flag allows an unprivileged user to appear as the **root user (UID 0)** inside the container. This relies on user namespace mapping (specifically `/etc/subuid` and `/etc/subgid` mappings on the host).

Fakeroot is useful when you:

- [Build the image](#build-the-image-with-fakeroot)
- [Run the container/overlay](#run-the-containeroverlay-with-fakeroot)

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

Then you can run commands that require root privileges, such as installing packages or changing file ownership inside the overlay.

```bash
# In side the container with fakeroot
apt update
apt install -y build-essential
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
