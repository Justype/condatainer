# Exec Troubleshooting

Two possible reasons for errors when mounting writable overlays `.img`:

- The image is locked by another process
- File permissions (internal or external) are incorrect

If you cannot write to an overlay image:

- Overlay is full
- Permission denied errors

## Table of Contents

- [Cannot mount overlay](#if-an-overlay-image-cannot-be-mounted)
- [Cannot write to overlay](#if-cannot-write-to-an-overlay-image)

## If an overlay image cannot be mounted

### 1. "Device or resource busy" (Used by another process)

If another process is already using the image, you cannot mount it as writable.

#### On a Headless Server

If you are running on a standard Linux server, identifying the locking process is straightforward:

```bash
# Option 1: List processes using the file
lsof /path/to/my_env.img
# Then use `kill <PID>` to terminate the process(es) listed.

# Option 2: Kill all processes accessing the file (Use with caution!)
fuser -k -9 /path/to/my_env.img
```

#### On an HPC Cluster

On a cluster, the image might be locked by a running job on a different node. `lsof` will not show processes running on other nodes.

Check your running jobs and cancel if needed:

```bash
# SLURM
squeue -u $USER
scancel <job_id>

# PBS/Torque
qstat -u $USER
qdel <job_id>

# LSF
bjobs
bkill <job_id>

# HTCondor
condor_q $USER
condor_rm <job_id>
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

#### Level 2: Inside the Overlay (Container)

**CondaTainer** can help you manage internal permissions. If you encounter permission errors when accessing files inside the overlay, you can change the ownership of all files to your user ID (UID) and group ID (GID):

```bash
condatainer overlay chown /path/to/my_env.img
```

```bash
# If you want fakeroot overlay
condatainer overlay chown --root /path/to/my_env.img
```

## If cannot write to an overlay image

### 1. Overlay is Full

If often happens when installing new packages inside the overlay. Increase the size of the overlay image:

```bash
condatainer overlay resize -s <new_size> /path/to/my_env.img
```

Example:

```bash
condatainer overlay resize -s 15G /path/to/my_env.img
```

### 2. Permission Denied Errors

If you encounter permission denied errors when writing to files inside the overlay, follow the steps in [Level 2: Inside the Overlay (Container)](#level-2-inside-the-overlay-container) above to change ownership of files inside the overlay.

```bash
# change inner files ownership
condatainer overlay chown /path/to/my_env.img
```
