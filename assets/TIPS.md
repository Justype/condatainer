# Tips

## Create writable images

Now `condatainer` can create and init images.

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
