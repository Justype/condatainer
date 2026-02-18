# Scheduler Resources: Nodes, Tasks, and CPUs

Understanding the difference between **nodes**, **tasks**, and **CPUs** helps you allocate resources correctly and avoid common mistakes on HPC clusters.

## Table of Contents

- [Concepts: Nodes, Tasks, and CPUs](#concepts)
- [Default: Use CPUs, Not Tasks](#default-use-cpus-not-tasks)
- [MPI Jobs: When to Use Multiple Tasks](#mpi-jobs-when-to-use-multiple-tasks)

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

## MPI Jobs: When to Use Multiple Tasks

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

