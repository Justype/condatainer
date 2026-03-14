# scheduler

HPC scheduler abstraction layer for job submission, script parsing, and resource management across SLURM, PBS, LSF, and HTCondor.

## Supported Schedulers

| Scheduler | Binary | Directive | Status |
|-----------|--------|-----------|--------|
| SLURM | `sbatch` | `#SBATCH` | Stable |
| PBS | `qsub` | `#PBS` | Experimental |
| LSF | `bsub` | `#BSUB` | Experimental |
| HTCondor | `condor_submit` | Native `.sub` format | Experimental |

## Architecture

```
scheduler.go        Core interface, types, resource resolution
registry.go         Thread-safe active scheduler singleton + debugMode
error.go            Structured error types
script_helpers.go   Shared parsing/formatting helpers (memory, time, flags, job header/footer)
slurm.go / pbs.go / lsf.go / htcondor.go   Scheduler implementations
```

## Interface

```go
type Scheduler interface {
    IsAvailable() bool
    IsInsideJob() bool
    ReadScriptSpecs(scriptPath string) (*ScriptSpecs, error)
    CreateScriptWithSpec(spec *JobSpec, outputDir string) (string, error)
    Submit(scriptPath string, dependencyJobIDs []string) (string, error)
    GetClusterInfo() (*ClusterInfo, error)
    GetInfo() *SchedulerInfo
    GetJobResources() *ResourceSpec // from scheduler env vars; zero fields = not set
}
```

## Key Types

**`ScriptSpecs`** — parsed job specification:
- `HasDirectives bool` — `false` when no `#SBATCH`/`#PBS`/`#BSUB` lines exist
- `Spec *ResourceSpec` — resource geometry; `nil` in passthrough mode
- `Control RuntimeConfig` — job name, stdout/stderr, email, partition

**`ResourceSpec`** — compute geometry: `Nodes`, `Ntasks`, `TasksPerNode` (0 = not set / freely distributed), `CpusPerTask`, `MemPerCpuMB`, `MemPerNodeMB`, `Time`, `Gpu`, `Exclusive`

| Memory field | Source | Meaning |
|---|---|---|
| `MemPerCpuMB` | SLURM `--mem-per-cpu`; LSF `rusage[mem=X]` ÷ CpusPerTask; PBS multi-chunk `+` (uniform per-cpu mem) | Per-CPU memory in MB |
| `MemPerNodeMB` | SLURM `--mem`; PBS `select=…:mem=X` (per-node) | Per-node total memory in MB |

`GetMemPerNodeMB()` checks `MemPerCpuMB × CpusPerTask × TasksPerNode` first (requires `TasksPerNode > 0`); otherwise returns `MemPerNodeMB` (0 if not set).

`GetMemPerTaskMB()` checks `MemPerCpuMB × CpusPerTask` first; otherwise returns `MemPerNodeMB ÷ TasksPerNode` (treats 0 as 1). Returns 0 if no memory is specified.

**`JobSpec`** — submission input: `Name`, `Command`, `Specs`, `DepJobIDs`, `Metadata`, `Array *ArraySpec`

**`ArraySpec`**: `InputFile`, `Count` (non-empty lines), `Limit` (max concurrent), `ArgCount`, `SampleArgs`, `BlankLines`

## Resource Spec Priority Chain

```
defaults  →  script directives (HasDirectives=true only)  →  job resources (non-zero fields)
```

```go
ResolveResourceSpec(jobResources, scriptSpecs)              // uses GetSpecDefaults() as baseline
ResolveResourceSpecFrom(defaults, jobResources, scriptSpecs) // custom baseline
(*ResourceSpec).Override(other *ResourceSpec)               // in-place merge
```

## Script Parsing

`ReadScriptSpecsFromPath` always returns non-nil `*ScriptSpecs`:

| State | `HasDirectives` | `Spec` | Meaning |
|---|---|---|---|
| No directives | `false` | non-nil | No scheduler lines found |
| Normal | `true` | non-nil | Parsed successfully |
| Passthrough | `true` | `nil` | Directives found, resource parse failed |

## Key Functions

```go
scheduler.Init("")                  // auto-detect and set active scheduler
scheduler.ActiveScheduler()         // get current scheduler

ReadScriptSpecsFromPath(path)       // always non-nil
HasSchedulerSpecs(specs)            // HasDirectives=true
IsPassthrough(specs)                // HasDirectives=true && Spec=nil
```

## Normalized Environment Variables

Normalized variables are **exported in the generated batch script** for LSF and HTCondor (which lack sufficient native environment variables). For SLURM and PBS, they are **computed on-demand** when entering containers via `condatainer run/exec`.

### Always Available

| Variable | Meaning | Always Set |
|----------|---------|------------|
| `NNODES` | Number of nodes | yes (default 1) |
| `NTASKS` | Total MPI tasks | yes (default 1) |
| `NCPUS` | CPUs per task (= `OMP_NUM_THREADS`) | yes (default 1) |
| `OMP_NUM_THREADS` | OpenMP thread count | yes (= `NCPUS`) |

### Conditionally Available

| Variable | OpenMP | MPI | Hybrid | Free-Dist (no Nodes) | Availability Logic |
|----------|--------|-----|--------|----------------------|-------------------|
| `NTASKS_PER_NODE` | ✓ | ✓ | ✓ | ✗ | `TasksPerNode > 0` (known layout) |
| `MEM` / `MEM_GB` | ✓ | ✓ | ✓ | ✓ | Always if memory specified |

**Memory Variable Logic:**
- `MEM` is memory **per task** in MB — always available when memory is specified
  - From `MemPerCpuMB`: `MEM = MemPerCpuMB × CpusPerTask`
  - From `MemPerNodeMB`: `MEM = MemPerNodeMB ÷ TasksPerNode` (treats 0 as 1)
- `MEM_GB` is `MEM ÷ 1024` (integer)

**Job Types:**
- **OpenMP**: 1 task, multiple CPUs per task (e.g., `NNODES=1, NTASKS=1, NCPUS=8`)
- **MPI**: Multiple tasks, 1 CPU per task, fixed layout (e.g., `NNODES=2, NTASKS=16, NTASKS_PER_NODE=8, NCPUS=1`)
- **Hybrid**: Multiple tasks, multiple CPUs per task, fixed layout (e.g., `NNODES=2, NTASKS=8, NTASKS_PER_NODE=4, NCPUS=4`)
- **Free-Dist (no Nodes)**: Only tasks, no node info (e.g., `NTASKS=7, NCPUS=1, MEM=8192`)

**Ceiling Division for Non-Uniform Distribution:**

When only get `Nodes` and `Ntasks` without `TasksPerNode`, CondaTainer will try to distribute the tasks by nodes as evenly as possible. (e.g. ntasks=7 and nodes=2 => tasks_per_node=4)

### Implementation Notes

- `NTASKS` uses `rs.Ntasks` directly when set; otherwise derives from `Nodes × TasksPerNode`
- `NTASKS_PER_NODE` is **absent** only when node count is not specified:
  - SLURM: `--ntasks=N` without `--nodes` (free distribution across any nodes)
  - LSF: `-n N` without `span[]` directive (free distribution)
  - PBS: requires explicit node count; single-task jobs default to 1 node

> [PBS 2024 User Manual Page 118: Hybrid MPI-OpenMP Jobs](https://help.altair.com/2024.1.0/PBS%20Professional/PBSUserGuide2024.1.pdf#page=118)

## Per-Scheduler Environment Variables

| Field | SLURM | PBS | LSF |
|-------|-------|-----|-----|
| Job ID | `SLURM_JOB_ID` | `PBS_JOBID` | `LSB_JOBID` |
| Nodes | `SLURM_JOB_NUM_NODES` | `PBS_NUM_NODES` | *(injected `NNODES`)* |
| Total tasks | `SLURM_NTASKS` | `PBS_NP` / `PBS_TASKNUM` | `LSB_DJOB_NUMPROC` / `LSB_MAX_NUM_PROCESSORS` (total slots) |
| Tasks per node | `SLURM_NTASKS_PER_NODE` | *(derived)* | *(injected `NTASKS_PER_NODE`)* |
| CPUs per task | `SLURM_CPUS_PER_TASK` | `NCPUS` (cpus per task) | *(injected `NCPUS`)* |
| Memory | `SLURM_MEM_PER_NODE` (MB) → `MemPerNodeMB` | `PBS_VMEM` (bytes) ÷ Nodes → `MemPerNodeMB` | `LSB_MAX_MEM_RUSAGE` (KB) → `MemPerCpuMB` |
| GPU | `CUDA_VISIBLE_DEVICES` count → `Gpu` | `CUDA_VISIBLE_DEVICES` count → `Gpu` | `CUDA_VISIBLE_DEVICES` count → `Gpu` |
| Array task index | `SLURM_ARRAY_TASK_ID` | `PBS_ARRAY_INDEX` | `LSB_JOBINDEX` |

HTCondor does not set standard resource environment variables.

## Error Types

- `SubmissionError` — job submission failed; includes scheduler name, job name, and output
- `ClusterError` — cluster info query failed; includes scheduler name and operation
- `ScriptCreationError` — batch script creation failed
- `TimeoutError` — scheduler command timed out

## Schedulers

Parallelism model fields: `Nodes`, `Ntasks` (0 = derive from `Nodes × TasksPerNode`), `TasksPerNode` (0 = free distribution), `CpusPerTask`.

### Job Type Directives

| Job type | SLURM | PBS | LSF | HTCondor |
|----------|-------|-----|-----|----------|
| **OpenMP** (`Ntasks=1, CpusPerTask=T`) | `--ntasks=1 --cpus-per-task=T` | `select=1:ncpus=T` | `-n T -R "span[hosts=1]"` | `request_cpus=T` |
| **MPI** (`CpusPerTask=1`) | `--nodes=N --ntasks-per-node=M --cpus-per-task=1` | `select=N:ncpus=M:mpiprocs=M` | `-n N*M -R "span[ptile=M]"` | *(single-task)* |
| **Hybrid** | `--nodes=N --ntasks-per-node=M --cpus-per-task=T` | `select=N:ncpus=M*T:mpiprocs=M:ompthreads=T` | `-n N*M -R "span[ptile=M] affinity[cores(T)]"` | *(single-task)* |

Free-distribution (non-divisible `Ntasks`): SLURM emits `--ntasks=Total`; LSF emits `-n Total`; PBS uses two-chunk `select=`.

### Slurm

Phase 2 resolution order (after collecting all directives):

1. **Nodes** — `--nodes`
2. **GPU** — `--gres=gpu:type:count` > `--gpus-per-node` > `--gpus`÷nodes > `--gpus-per-task`×tasks
3. **Tasks** — `--ntasks`; `--ntasks-per-node`; `--ntasks-per-gpu` × gpuCount × nodes (requires GPU)
4. **CPU** — `--cpus-per-task`; `--cpus-per-gpu` × gpuCount ÷ tasksPerNode
5. **Memory** — `--mem` → `MemPerNodeMB`; `--mem-per-cpu` → `MemPerCpuMB`; `--mem-per-gpu` × gpu count → `MemPerNodeMB`
6. **Time** — `--time` (independent)

Blacklisted (passthrough): `--gpus-per-socket`, `--sockets-per-node`, `--cores-per-socket`, `--threads-per-core`, `--ntasks-per-socket`, `--ntasks-per-core`, `--distribution`

### PBS

PBS Pro/OpenPBS only; Torque not supported.

- Resources: `select=N:ncpus=M:mpiprocs=P:mem=Xgb:ngpus=G`
- Old-style standalone resources (`-l nodes=`, `-l ncpus=`, etc.) → **passthrough**
- Multi-chunk `+`: all chunks must share the same `CpusPerTask = ncpus/mpiprocs`; `ngpus=` in any chunk → passthrough

### LSF

- No node-count directive; nodes derived from `-n / span[ptile=M]`
- If `span[hosts=1]`, `-n` = total CPUs; otherwise `-n` = total tasks
- Memory: `rusage[mem=X]` is per-slot → `MemPerCpuMB = rusage_mem / affinity.cores` (default 1); without `-n` → `MemPerNodeMB`
- `-M` (KB, ulimit) used as fallback; same division rule applies
- GPU: `-gpu "num=N:gmodel=T"` or `-R "rusage[ngpus_physical=N]"`

**LSF `-n` semantics:**

| `-R` directive | `-n` meaning | `Nodes` |
|----------------|-------------|---------|
| `span[hosts=1]` | Total CPUs (= `CpusPerTask`) | 1 |
| `span[ptile=M]` | Total tasks (= `Ntasks`) | `Ntasks / ptile` |
| `affinity[cores(T)]` only | Total tasks | 0 |
| *(none)* | Total tasks | 0 |

### HTCondor

- Parses native `.sub` files; `executable` → `ScriptSpecs.ScriptPath`
- Resource keys: `request_cpus`, `request_memory`, `request_gpus`, `+MaxRuntime` (seconds)
- Always single-task; no dependency support

## Array Jobs

Pass `JobSpec.Array = &ArraySpec{...}` to `CreateScriptWithSpec`.

- Scheduler directive silences per-subjob stdout (`--output=/dev/null`).
- `writeArrayBlock` extracts the current line via `sed -n "${IDX}p"` → `ARRAY_ARGS`; redirects output to `{logDir}/{name}_{idx}_{tag}.log`.
- HTCondor uses `queue array_args from {inputFile}`.

| Scheduler | Directive |
|-----------|-----------|
| SLURM | `#SBATCH --array=1-N%limit` |
| PBS | `#PBS -J 1-N%limit` |
| LSF | `#BSUB -J name[1-N]%limit` |
| HTCondor | `queue array_args from file` |

## Dependency Formats

| Scheduler | Format |
|-----------|--------|
| SLURM | `--dependency=afterok:ID1:ID2` |
| PBS | `-W depend=afterok:ID1:ID2` |
| LSF | `-w "done(ID1) && done(ID2)"` |
| HTCondor | not supported (requires DAGMan) |

## Time Formats

| Scheduler | Format | Example |
|-----------|--------|---------|
| SLURM | `D-HH:MM:SS` or `HH:MM:SS` | `1-12:00:00` |
| PBS | `HH:MM:SS` | `04:00:00` |
| LSF | `HH:MM` | `04:00` |
| HTCondor | integer seconds (`+MaxRuntime`) | `14400` |

## Scheduler Documentation

- [Slurm Sbatch](https://slurm.schedmd.com/sbatch.html)
- [LSF 10.1.0](https://www.ibm.com/docs/en/spectrum-lsf/10.1.0)
- [PBS Pro 2024.1 User Guide](https://help.altair.com/2024.1.0/PBS%20Professional/PBSUserGuide2024.1.pdf)
- [HTCondor Latest](https://htcondor.readthedocs.io/en/latest/)
