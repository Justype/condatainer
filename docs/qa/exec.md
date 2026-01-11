# Exec Troubleshooting

If you encounter errors when trying to run or mount a CondaTainer overlay, it is usually due to one of two reasons: the image is locked by another process, or the file permissions (internal or external) are incorrect.

## If an overlay image cannot be mounted

### 1. "Device or resource busy" (Used by another process)

Apptainer requires exclusive write access to mount an overlay. If another process is already using the image (even in read-only mode), you cannot mount it as writable.

#### On a Headless Server
If you are running on a standard Linux server, identifying the locking process is straightforward:

```bash
# Option 1: List processes using the file
lsof /path/to/my_env.img
# Then use `kill <PID>` to terminate the process(es) listed.

# Option 2: Kill all processes accessing the file (Use with caution!)
fuser -k -9 /path/to/my_env.img
```

#### On an HPC Cluster (Slurm)

On a cluster, the image might be locked by a running job on a different node. lsof will not show processes running on other nodes.

Check your running jobs:

```bash
squeue -u $USER
```

If you find a job that might be using the overlay, you can cancel it:

```bash
scancel <job_id>
```

#### After killing

If a job crashed or was killed forcibly while writing to the image, the overlay filesystem might be marked as "dirty." You should run a filesystem check to repair it.

```bash
e2fsck -p /path/to/my_env.img
```

If `-p` does not fix all issues, you can run it without `-p` to enter interactive mode.

```bash
e2fsck /path/to/my_env.img
```

### 2. Permission Denied Errors

#### Level 1: The Image File (Host)

Ensure you have read and write permissions on the overlay image file itself.

```bash
ls -l /path/to/my_env.img
```

If you do not have the necessary permissions, and is not the owner, you may need to contact the system administrator or the owner of the file to adjust the permissions.

```bash
# only run if you are the owner
chmod u+rw /path/to/my_env.img
```

#### Level 2: Inside the Overlay (Container)

**CondaTainer** can help you manage internal permissions. If you encounter permission errors when accessing files inside the overlay, you can change the ownership of all files to your user ID (UID) and group ID (GID):

```bash
condatainer overlay chown /path/to/my_env.img
```

You can also use `--fakeroot` option when running `condatainer exec` to avoid permission issues:

```bash
condatainer exec --fakeroot -o /path/to/my_env.img <command>
```
