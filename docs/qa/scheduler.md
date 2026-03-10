# Scheduler Resources

Understanding the difference between **nodes**, **tasks**, and **CPUs** helps you allocate resources correctly and avoid common mistakes on HPC clusters.

## Table of Contents

- [Concepts: Nodes, Tasks, and CPUs](#concepts)
  - [Single-Task vs Multi-Task (MPI) Jobs](#single-task-vs-multi-task-mpi-jobs)
  - [HTCondor: Single-Task Only](#htcondor-single-task-only)
- [How Schedulers Define Resources](#how-schedulers-define-resources)
- [How CondaTainer Handles Directives](#how-condatainer-handles-directives)
  - [Normalizing Directives and Passthrough Mode](#normalizing-directives-and-passthrough-mode)
  - [Cross-Scheduler Translation](#cross-scheduler-translation)
  - [MPI Auto-Detection](#mpi-auto-detection)
- [Manual Submission](#manual-submission)

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

**Most jobs — Python, R, bash scripts — should request multiple treads, not multiple tasks.**

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

### HTCondor: Single-Task Only

HTCondor is inherently single-task (`universe = vanilla` is assumed for all jobs), so `nodes` and `tasks-per-node` do not apply.

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

CondaTainer supports three job dependency flags, each mapping to a native scheduler dependency type:

| Flag | Runs when upstream job… | SLURM / PBS | LSF |
|---|---|---|---|
| `--afterok IDS` | Succeeds | `afterok` | `done()` |
| `--afternotok IDS` | Fails | `afternotok` | `exit()` |
| `--afterany IDS` | Finishes (any outcome) | `afterany` | `ended()` |

```bash
TRIM=$(condatainer run trim.sh sample1)
ALIGN=$(condatainer run --afterok "$TRIM" align.sh sample1)
condatainer run --afterok "$ALIGN" quant.sh sample1
```

Multiple job IDs can be joined with colons: `--afterok 123:456:789`. All three flags can be combined in one submission.

> **Note**: `DAGMan` is not supported. Please use `DAGMan` directly for complex workflows on HTCondor clusters.

### Array Jobs + Chaining

`--array` and `--afterok` can be combined: each stage submits an array job and waits for the previous stage to **succeed** before starting. Use `--afterany` to proceed even if some subjobs failed. A final single job can collect results after all subjobs complete.

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

Each scheduler has its own syntax for the same underlying resource concepts. All directives below appear as in-script comments (`#SBATCH`, `#PBS`, `#BSUB`) for SLURM/PBS/LSF, or in a separate `.sub` file for HTCondor.

| Concept | SLURM | PBS (`#PBS -l`) | LSF | HTCondor |
|---------|-------|-----------------|-----|----------|
| Nodes | `--nodes=N` | `select=N:...` | *(derived from `-n` ÷ `ptile`)* | *(single-task only)* |
| Total tasks | `--ntasks=N` | *(Chunks × mpiprocs)* | `-n N` | *(not applicable)* |
| Tasks per node | `--ntasks-per-node=N` | `select=N:mpiprocs=P:...` | `-R "span[ptile=N]"` | *(not applicable)* |
| CPUs per task | `--cpus-per-task=N` | `select=1:ncpus=N:mpiprocs=1:...` | `-R "affinity[cores(N)]"` | `request_cpus = N` |
| Memory | `--mem=N` / `--mem-per-cpu=N` | `select=1:mem=N:...` (per node) | `-R "rusage[mem=N]"` (per slot) | `request_memory = N` |
| GPU | `--gpus-per-node=N` | `select=1:ngpus=N:...` | `-gpu "num=N"` | `request_gpus = N` |
| Wall time | `--time=HH:MM:SS` | `walltime=HH:MM:SS` | `-W HH:MM` | `+MaxRuntime = N` (sec) |
| Job name | `--job-name=NAME` | `-N NAME` | `-J NAME` | *(not supported)* |
| Stdout / Stderr | `-o` / `-e` | `-o` / `-e` | `-o` / `-e` | `output` / `error` |
| Email type | `--mail-type=TYPE` | `-m <aben>` | `-B` / `-N` | `notification = ...` |
| Email address | `--mail-user=ADDR` | `-M ADDR` | `-u ADDR` | `notify_user = ADDR` |

**PBS** packs chunks, tasks-per-node, CPUs, and memory into a single `-l select=N:ncpus=M:mpiprocs=P:mem=X` command. Multiple chunks separated by `+` allow different resource shapes on different nodes.

**LSF** uses a single `-R` string to express node constraint (`span[hosts=1]` for single-node), task distribution (`span[ptile=N]`), CPU affinity (`affinity[cores(N)]`), and memory allocation (`rusage[mem=N]`). `-n N` sets the total slot count. `-M N` is a per-slot memory ulimit and passes through separately.

### GPU support

GPU resources are normalized across schedulers. SLURM and LSF support specifying a GPU model; PBS and HTCondor translate to count only.

| Scheduler | Num | Model | Directive |
|-----------|-----|-------|-----------|
| SLURM | Yes | Yes | `--gpus-per-node=<model>:<num>` |
| PBS Pro | Yes | No | `ngpus=<num>` (in `-l select=...`) |
| LSF | Yes | Yes | `-gpu "num=<num>:gmodel=<model>"` |
| HTCondor | Yes | No | `request_gpus = <num>` |

## How CondaTainer Handles Directives

When `condatainer run` is invoked, it reads the script's scheduler directives and processes them.

### Normalizing Directives and Passthrough Mode

CondaTainer parses all directives into a **normalized internal representation**: 

```go
type ResourceSpec struct {
	Nodes        int           // Number of nodes
	Ntasks       int           // Total MPI tasks (0 = not set; use Nodes*TasksPerNode)
	TasksPerNode int           // MPI ranks per node (0 = non-uniform distribution, e.g. multi-chunk PBS)
	CpusPerTask  int           // CPU threads per task (OpenMP)
	MemPerCpuMB  int64         // RAM per logical CPU in MB (SLURM --mem-per-cpu, LSF rusage[mem])
	MemPerNodeMB int64         // Total RAM per node in MB (SLURM --mem, PBS mem=)
	Gpu          *GpuSpec      // GPU requirements (nil = no GPU; Count = per node)
	Time         time.Duration // Job walltime limit
	Exclusive    bool          // Request exclusive node access (no other jobs on the same node)
}
```

If a directive cannot be represented in this schema — for example, low-level CPU topology constraints — CondaTainer enters **passthrough mode**: the script's directives cannot be regenerated and `condatainer run` will **reject the submission with an error**, directing you to submit manually with the native scheduler command. See [Manual Submission](#manual-submission) below.

**For MPI jobs, all tasks must share the same CPUs-per-task and memory-per-CPU.** CondaTainer stores a single CPUs-per-task and a memory-per-CPU for the entire job. PBS can use multi chunks (`select=...+...`) for uneven task distributions.

### Cross-Scheduler Translation

Once normalized, CondaTainer can emit directives in any target scheduler's syntax, enabling scripts written for one scheduler to run on a cluster running a different one.

**Example: SLURM script submitted on a PBS cluster**

```bash
# Original SLURM directives:
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
```

```bash
# PBS equivalent generated at submission time:
#PBS -l select=2:ncpus=8:mpiprocs=4:mem=32768mb
#PBS -l walltime=02:00:00
```

Uneven MPI distributions are reconstructed as multi-chunk PBS select:

```bash
# SLURM: 8 tasks across 3 nodes (3+3+2)
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=3
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100M
```

```bash
# PBS: two chunks — 2 full nodes + 1 remainder
#PBS -l select=2:ncpus=3:mpiprocs=3:mem=300mb+1:ncpus=2:mpiprocs=2:mem=200mb
```

Two fields are intentionally dropped during translation:

- **Partition / queue** — queue names are site-specific.
- **Unrecognized flags** — directives with no cross-scheduler equivalent are dropped with a warning.

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

## Manual Submission

When a script triggers passthrough mode, `condatainer run` rejects it with an error. You must submit directly with the native scheduler command (`sbatch`, `qsub`, `bsub`).

To still run your payload inside a CondaTainer container in this case, use the **self-bootstrap pattern**: the script detects whether it is already running inside a container (`$IN_CONDATAINER`) and, if not, re-invokes itself via `condatainer run`. CondaTainer then resolves the `#DEP:` tags and starts the container.

```bash
#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=200M
#SBATCH --time=02:00:00

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

Submit with the native scheduler — **not** `condatainer run`:

```bash
sbatch my_script.sh
```
