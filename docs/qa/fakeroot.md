# Fakeroot

The `--fakeroot` flag allows an unprivileged user to appear as the **root user (UID 0)** inside the container. This relies on user namespace mapping (specifically `/etc/subuid` and `/etc/subgid` mappings on the host).

Fakeroot is useful when you:

- [Build the image with fakeroot](#build-the-image-with-fakeroot)
- [Run the container/overlay with fakeroot](#run-the-containeroverlay-with-fakeroot)

## Build the image with fakeroot

On HPC systems, you don't have root privileges. Using `--fakeroot` allows you to create an overlay image that behaves like root inside the container.

By default, CondaTainer uses `--fakeroot` when creating overlays using `.def` files or pulling from remote sources.

## Run the container/overlay with fakeroot

When executing commands inside the container, you can use the `--fakeroot` flag to gain root-like privileges inside the container.

```bash
condatainer exec --fakeroot -o /path/to/my_env.img <command>
```

Then you can run commands that require root privileges, such as installing packages or changing file ownership inside the overlay.

```bash
# In side the container with fakeroot
apt update
apt install -y build-essential
```

```{note}
It works on personal machine. But on some HPC systems, this may not work due to apptainer or system configuration. In that case, just stick to normal user mode and use `condatainer overlay chown` to fix permission issues inside the overlay.

Also, RStudio Sever may not work with UID less than 1000.
```
