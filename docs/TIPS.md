# Tips

## Create writable images

Now `condatainer` can create and init writable images.

```bash
# create a 5GB image and create a conda env inside
condatainer o -s 5120 test.img
```

Then you will in the writable container shell. You will get a prompt like:

```
Overlay env:
  CNT_CONDA_PREFIX: /ext3/test

CNT /someplace> 
```

Then you can use my help commands to create/install/update/remove conda packages inside the writable container.

```bash
# install more packages
mm-install numpy pandas

# update packages
mm-update

# remove packages
mm-remove pandas

# list installed packages
mm-list

# clean all conda cache
mm-clean -ay

# export the env
mm-export --no-builds > my_env.yaml
```

## Create read-only images with env yaml

You can create read-only images with a conda env yaml file.

```bash
condatainer create -p env_prefix -f my_env.yaml
```

Then you will get an image with the conda env installed at `env_prefix.sqf`.

```bash
condatainer exec -o env_prefix.sqf bash
```

## CondaTainer is compatible with module systems

CondaTainer will scan your script for `module load` or `ml` commands and mount the corresponding overlays automatically.

**Example:**

```bash
#!/bin/bash
module load bcftools/1.22
bcftools --version
```

When you run the script with `condatainer run`, it will automatically mount the `bcftools/1.22` overlay.

```bash
condatainer run my_script.sh
```
