# Tips

## Create writable images

1. create ext3 image file using `dd` or `apptainer`/`singularity`
2. load the image using `condatainer exec -w -o <image_file>`

```bash
apptainer overlay create -s 5120 -S test.img # create 5GB sparse ext3 image

condatainer exec -w test.img
# equivalent to: condatainer exec --writable -o test.img bash
```

Then you will in the writable container shell. You will get a prompt like:

```
Overlay env:
  CNT_CONDA_PREFIX: /ext3/test

CNT /someplace> 
```

Then you can use my help commands to create/install/update/remove conda packages inside the writable container.

> [!NOTE]
> For empty writable images, run `mm-create <any_package>` first to create a conda env under `$CNT_CONDA_PREFIX`. Then you can use other commands to manage packages inside the writable image.

```bash
# initialize writable conda env with "star" package
mm-create star

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
