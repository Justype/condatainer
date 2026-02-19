# Scheduler Resources

Understanding the difference between **nodes**, **tasks**, and **CPUs** helps you allocate resources correctly and avoid common mistakes on HPC clusters.

## Table of Contents

- [Concepts: Nodes, Tasks, and CPUs](#concepts)
- [Default: Use CPUs, Not Tasks](#default-use-cpus-not-tasks)
- [MPI Jobs: When to Use Multiple Tasks](#mpi-jobs-when-to-use-multiple-tasks)
- [Cross-Scheduler Translation](#cross-scheduler-translation)

## Concepts

| Term | SLURM directive | Meaning |
|------|----------------|---------|
| **Node** | `--nodes` | A physical compute machine |
| **Task** | `--ntasks`, `--ntasks-per-node` | An MPI process (one independent program rank) |
| **CPU** | `--cpus-per-task` | Threads available to a single task |

A job with `--nodes=2 --ntasks-per-node=4 --cpus-per-task=2` runs:
- 2 machines
- 4 MPI processes per machine (8 total)
- Each process may use up to 2 threads

For PBS/LSF the terminology differs slightly but the concept is the same:

| Concept | SLURM | PBS | LSF |
|---------|-------|-----|-----|
| Nodes | `--nodes=N` | `select=N` | `-m host1+host2` |
| Tasks per node | `--ntasks-per-node=N` | `mpiprocs=N` | `-R "span[ptile=N]"` |
| CPUs per task | `--cpus-per-task=N` | `ncpus=N` | `-n N` |

## Default: Use CPUs, Not Tasks

**Most single-program jobs — Python, R, bash scripts — should request multiple CPUs (`--cpus-per-task`), not multiple tasks.**

Multiple tasks launch multiple independent copies of your program (MPI ranks). Unless your code explicitly calls `MPI_Init` / `from mpi4py import MPI`, extra tasks are wasted allocation.

**Correct: threaded job**

```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2:00:00

#DEP: myenv.img

python train_model.py
```

Inside the container, thread-parallel libraries will use all 8 CPUs.

**Wrong: accidental multi-task**

```bash
#SBATCH --ntasks=8   # ❌ Launches 8 copies of python — each does the full job
```

This starts 8 separate Python processes, each running `train_model.py` from scratch wasting 7× the resources.

## MPI Jobs

Use multiple tasks **only** when your program uses MPI for inter-process communication.

```bash
#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4   # 8 MPI ranks total across 2 nodes
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1:00:00

#DEP: mpi.img

python mpi_simulation.py
```

When `condatainer run` detects `ntasks > 1`, it **automatically wraps the command with `mpiexec`**:

```bash
# Generated job command (via module system):
module purge && module load openmpi/4.1.5 && mpiexec condatainer run mpi_job.sh
```

Each MPI rank launches its own container, all sharing the same MPI communicator via SLURM's process management interface.

````{important}
You need to have the same major and minor version of OpenMPI installed inside the container as on the host.

```bash
ml av openmpi
# openmpi/4.1.5
condatainer e mpi.img -- mm-install mpi4py openmpi=4.1 -y
```
````

### Hybrid MPI + threading

Some programs combine MPI ranks with per-rank threads:

```bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2   # 2 MPI ranks per node
#SBATCH --cpus-per-task=4     # 4 OpenMP threads per rank

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
```

## CondaTainer MPI Limitations

**Tasks must distribute evenly across nodes.**

This constraint exists to enable **cross-scheduler translation**: schedulers like PBS use `ntasks-per-node` (`mpiprocs=N`) rather than a global task count, so CondaTainer must be able to convert `--ntasks` into a per-node value.

If the division is uneven, that translation is ambiguous and CondaTainer falls back to **passthrough mode**, and `mpiexec` wrapping is disabled.

**Will work** (8 tasks across 2 nodes = 4 per node):

```bash
#SBATCH --nodes=2
#SBATCH --ntasks=8
```

**Passthrough mode triggered** (3 tasks across 2 nodes — uneven):

```bash
#SBATCH --nodes=2
#SBATCH --ntasks=3   # ⚠ 3 % 2 != 0 → mpiexec wrapping disabled
```

Use `--ntasks-per-node` instead of `--ntasks` to guarantee even distribution and avoid passthrough mode:

```bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4   # Always evenly distributed, translates cleanly to PBS/LSF
```

See [Passthrough Mode Workaround](#passthrough-mode-workaround) for a workaround if you need to use an unsupported topology directive.

## Cross-Scheduler Translation

CondaTainer can run a script written for one scheduler on a cluster running a different scheduler. When `condatainer run` detects that the script's directives (e.g., `#SBATCH`) do not match the host scheduler, it **automatically translates** the resource spec into the host's format.

Supported schedulers: **SLURM**, **PBS/Torque**, **LSF/Spectrum LSF**, **HTCondor**.

### How it works

CondaTainer parses all scheduler directives into a normalized internal representation (nodes, tasks-per-node, CPUs-per-task, memory, GPU, time, job control).

When generating the submission script for the host scheduler, it writes those fields out in the host's syntax.

**Example: SLURM script submitted on a PBS cluster**

```bash
#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
```

CondaTainer translates this into the PBS equivalent at submission time:

```bash
#PBS -l select=2:ncpus=2:mpiprocs=4:mem=32768mb
#PBS -l walltime=02:00:00
```

### What is translated

| Field | SLURM | PBS | LSF |
|-------|-------|-----|-----|
| Nodes | `--nodes` | `select=N` | `span[hosts=N]` |
| Tasks per node | `--ntasks-per-node` | `mpiprocs=N` | `span[ptile=N]` |
| CPUs per task | `--cpus-per-task` | `ncpus=N` | `-n N` |
| Memory per node | `--mem` | `mem=NMB` | `-M NMB` |
| GPU | `--gpus-per-node:N` | `ngpus=N` | `-gpu "num=N"` |
| Wall time | `--time` | `walltime=HH:MM:SS` | `-W HH:MM` |
| Job name | `-J` / `--job-name` | `-N` | `-J` |
| Stdout / Stderr | `-o` / `-e` | `-o` / `-e` | `-o` / `-e` |
| Email | `--mail-type/user` | `-m` / `-M` | `-B` / `-N` / `-u` |

### What is NOT translated

Two fields are intentionally dropped during cross-scheduler translation:

- **Partition / queue** — queue names are site-specific and have no equivalent on the target scheduler. Set the target partition via `condatainer scheduler -p <partition>` or `scheduler.partition` in the config.
- **Scheduler-specific flags** — unrecognized directives (e.g., `#SBATCH --constraint=...`) are dropped with a warning. They cannot be expressed in the target scheduler's syntax.

### Passthrough mode

If CondaTainer cannot safely translate the resource spec — due to topology directives that have no cross-scheduler equivalent — it enters **passthrough mode**: the job is submitted but no resource directives are written. The original unrecognized flags are forwarded as-is.

Directives that trigger passthrough mode on SLURM scripts:

```
--gpus-per-socket   --sockets-per-node   --cores-per-socket
--threads-per-core  --ntasks-per-socket  --ntasks-per-core
--distribution
```

### HTCondor

HTCondor uses native `.sub` submit files rather than in-script directives. CondaTainer only parses HTCondor specs from files with a `.sub` extension.

HTCondor is inherently single-node, so `nodes` and `tasks-per-node` are not translated. (`universe = vanilla` is assumed for all jobs.)

Job dependencies (after-OK chains) are not supported by HTCondor's simple submission interface (`DAGMan` is not supported in current version).

### Passthrough Mode Workaround

Just don't use `condatainer run` for the initial submission. Submit the job with `sbatch` as normal, and CondaTainer will detect the `IN_CONDATAINER` environment variable and re-run the script in container mode on the compute node.

```bash
#!/bin/bash
#SBATCH --nodes=3
#SBATCH --ntasks=8   # Triggers passthrough mode
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00

#DEP: mpi.img

if [ -z "$IN_CONDATAINER" ] && command -v condatainer >/dev/null 2>&1; then
    if [ -n "$SLURM_JOB_ID" ]; then
        FULL_COMMAND=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}' | head -n 1)
        ORIGINAL_SCRIPT_PATH=$(echo "$FULL_COMMAND" | awk '{print $1}')
    else
        ORIGINAL_SCRIPT_PATH=$(realpath "$0")
    fi
    module purge && module load openmpi/4.1.5 && \
        mpiexec condatainer run "$ORIGINAL_SCRIPT_PATH"
    exit $?
fi

python my_script.py
```

Then use `sbatch my_script.sh` to submit as normal.
