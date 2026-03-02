# Scheduler Resources

Understanding the difference between **nodes**, **tasks**, and **CPUs** helps you allocate resources correctly and avoid common mistakes on HPC clusters.

## Table of Contents

- [Concepts: Nodes, Tasks, and CPUs](#concepts)
  - [Single-Task vs Multi-Task (MPI) Jobs](#single-task-vs-multi-task-mpi-jobs)
  - [HTCondor: Single-Node Only](#htcondor-single-node-only)
- [How Schedulers Define Resources](#how-schedulers-define-resources)
- [How CondaTainer Handles Directives](#how-condatainer-handles-directives)
  - [Normalizing Directives and Passthrough Mode](#normalizing-directives-and-passthrough-mode)
  - [Cross-Scheduler Translation](#cross-scheduler-translation)
  - [MPI Auto-Detection](#mpi-auto-detection)
- [Passthrough Mode Workaround](#passthrough-mode-workaround)

## Concepts

| Term | Meaning |
|------|---------|
| **Node** | A physical compute machine |
| **Task** | An MPI process (one independent program rank) |
| **CPU** | Threads available to a single task |

A job with 2 nodes, 4 tasks per node, and 2 CPUs per task runs:
- 2 machines
- 4 MPI processes per machine (8 total)
- Each process may use up to 2 threads

### Single-Task vs Multi-Task (MPI) Jobs

**Most jobs — Python, R, bash scripts — should request multiple CPUs, not multiple tasks.**

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

This starts 8 separate Python processes, each running `train_model.py` from scratch, wasting 7× the resources.

**Correct: MPI job**

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

### HTCondor: Single-Node Only

HTCondor is inherently single-node (`universe = vanilla` is assumed for all jobs), so `nodes` and `tasks-per-node` do not apply.

When CondaTainer reads `.sub` files, it will get the script path by `executable = ...` and parse other resource requests. Then CondaTainer will read the overlay requirements from the script `#DEP:` tags and launch the container with the appropriate resources.

### Array Jobs

Use `--array` with `condatainer run` to submit the same script over a list of inputs. Each line in the input file becomes one subjob; its space-separated tokens arrive as positional arguments (`$1`, `$2`, …) inside the script.

```
sample1 condition_A
sample2 condition_B
```

The file should have no blank lines, and all non-empty lines must have the same number of shell-split tokens (quoted strings count as one). Blank lines are flagged on `--dry-run` and block submission — remove them before submitting.

```bash
condatainer run --array samples.txt quant.sh
condatainer run --array samples.txt --array-limit 4 quant.sh  # max 4 concurrent
condatainer run --dry-run --array samples.txt quant.sh        # preview without submitting
```

### Job Chaining

CondaTainer supports job dependencies via `--afterok` (SLURM-style). You can submit a job and capture its ID for downstream submission.

```bash
TRIM=$(condatainer run trim.sh sample1)
ALIGN=$(condatainer run --afterok "$TRIM" align.sh sample1)
condatainer run --afterok "$ALIGN" quant.sh sample1
```

Multiple job IDs can be joined with colons: `--afterok 123:456:789`.

> **Note**: `DAGMan` is not supported. Use `DAGMan` directly for complex workflows on HTCondor clusters.

### Array Jobs + Chaining

`--array` and `--afterok` can be combined: each stage submits an array job and waits for the previous stage to finish before starting. A final single job can collect results after all subjobs complete.

```bash
# Stage 1: trim all samples (no dependency)
JOB=$(condatainer run --array samples.txt --array-limit 10 trim.sh)
# Stage 2: align — waits for ALL trim subjobs to finish
JOB=$(condatainer run --array samples.txt --array-limit 10 --afterok "$JOB" align.sh)
# Stage 3: quant — waits for ALL align subjobs to finish
JOB=$(condatainer run --array samples.txt --array-limit 10 --afterok "$JOB" quant.sh)
# Final: single job collecting results — waits for ALL quant subjobs
condatainer run --afterok "$JOB" collect_results.sh samples.txt
```

## How Schedulers Define Resources

Each scheduler has its own syntax for the same underlying resource concepts:

| Concept | SLURM | PBS | LSF | HTCondor |
|---------|-------|-----|-----|----------|
| Nodes | `--nodes=N` | `select=N` | `span[hosts=N]` | *(single-node only)* |
| Tasks per node | `--ntasks-per-node=N` | `mpiprocs=N` | `span[ptile=N]` | *(not applicable)* |
| CPUs per task | `--cpus-per-task=N` | `ncpus=N` | `-n N` | `request_cpus = N` |
| Memory per node | `--mem=N` | `mem=NMB` | `-M NMB` | `request_memory = N` |
| GPU | `--gpus-per-node=N` | `ngpus=N` | `-gpu "num=N"` | `request_gpus = N` |
| Wall time | `--time=HH:MM:SS` | `walltime=HH:MM:SS` | `-W HH:MM` | `+MaxRunningTime = N` |
| Job name | `-J` / `--job-name` | `-N` | `-J` | `MyDescription = ...` |
| Stdout / Stderr | `-o` / `-e` | `-o` / `-e` | `-o` / `-e` | `output` / `error` |
| Email | `--mail-type/user` | `-m` / `-M` | `-B` / `-N` / `-u` | `notify_user` |

For SLURM/PBS/LSF, directives appear as in-script comments (`#SBATCH`, `#PBS`, `#BSUB`). For HTCondor, they appear in a separate `.sub` file.

### GPU support

GPU resources are normalized across schedulers. SLURM and LSF support specifying a GPU model; PBS and HTCondor translate to count only.

| Scheduler | Num | Model | Directive |
|-----------|-----|-------|-----------|
| SLURM | Yes | Yes | `--gpus-per-node=<model>:<num>` |
| PBS Pro | Yes | No | `ngpus=<num>` (in `-l select=...`) |
| LSF | Yes | Yes | `-gpu "num=<num>:gmodel=<model>"` |
| HTCondor | Yes | No | `request_gpus = <num>` |

## How CondaTainer Handles Directives

When `condatainer run` is invoked, it reads the script's scheduler directives and processes them in three steps.

### Normalizing Directives and Passthrough Mode

CondaTainer first parses all directives into a **normalized internal representation** (nodes, tasks-per-node, CPUs-per-task, memory, GPU, time, job control).

If a directive cannot be normalized — because it has no cross-scheduler equivalent — CondaTainer enters **passthrough mode**: no resource directives are written to the submission script and the unrecognized flags are forwarded as-is.

SLURM directives that trigger passthrough mode:

```
--gpus-per-socket   --sockets-per-node   --cores-per-socket
--threads-per-core  --ntasks-per-socket  --ntasks-per-core
--distribution
```

**MPI task distribution must be even.** Schedulers like PBS use `mpiprocs=N` (tasks *per node*) rather than a global task count, so CondaTainer must convert `--ntasks` into a per-node value. If the division is uneven, it cannot be represented and passthrough mode is triggered.

```bash
# Will work — 8 tasks across 2 nodes = 4 per node
#SBATCH --nodes=2
#SBATCH --ntasks=8

# Passthrough mode triggered — 3 % 2 != 0
#SBATCH --nodes=2
#SBATCH --ntasks=3
```

Use `--ntasks-per-node` to guarantee even distribution:

```bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4   # Always evenly distributed, translates cleanly to PBS/LSF
```

### Cross-Scheduler Translation

Once directives are normalized, CondaTainer can emit them in any supported scheduler's syntax. This allows a script written for one scheduler to run on a cluster running a different one.

**Example: SLURM script submitted on a PBS cluster**

```bash
# Original SLURM directives:
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

Two fields are intentionally dropped during translation:

- **Partition / queue** — queue names are site-specific. Set the target partition via `condatainer scheduler -p <partition>` or `scheduler.partition` in the config.
- **Unrecognized flags** — directives with no target-scheduler equivalent (e.g., `#SBATCH --constraint=...`) are dropped with a warning.

### MPI Auto-Detection

When CondaTainer detects `ntasks > 1`, it **automatically wraps the command with `mpiexec`**:

```bash
# Generated job command:
mpiexec condatainer run mpi_job.sh
```

`mpiexec` must be available in `PATH` at submission time (e.g. load your MPI module before calling `condatainer run`). If it is not found, submission fails with an error.

Each MPI rank launches its own container, all sharing the same MPI communicator via the scheduler's process management interface.

````{important}
The OpenMPI version inside the container must match the major and minor version on the host.

```bash
ml av openmpi
# openmpi/4.1.5
condatainer e mpi.img -- mm-install mpi4py openmpi=4.1 -y
```
````

**Hybrid MPI + threading** combines MPI ranks with per-rank threads:

```bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2   # 2 MPI ranks per node
#SBATCH --cpus-per-task=4     # 4 OpenMP threads per rank

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
```

## Passthrough Mode Workaround

If you need to use a directive that triggers passthrough mode, bypass `condatainer run` for the initial submission. Submit with `sbatch` directly — CondaTainer detects the `IN_CONDATAINER` environment variable on the compute node and re-runs the script in container mode automatically.

```bash
#!/bin/bash
#SBATCH --nodes=3
#SBATCH --ntasks=8   # Triggers passthrough mode — uneven across 3 nodes
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
