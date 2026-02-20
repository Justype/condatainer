# scheduler

HPC scheduler abstraction layer for job submission, script parsing, and resource management across SLURM, PBS/Torque, LSF, and HTCondor.

## Supported Schedulers

| Scheduler | Binary | Directive | Status |
|-----------|--------|-----------|--------|
| SLURM | `sbatch` | `#SBATCH` | Stable |
| PBS/Torque | `qsub` | `#PBS` | Experimental |
| LSF | `bsub` | `#BSUB` | Experimental |
| HTCondor | `condor_submit` | Native `.sub` format | Experimental |

## Architecture

```
scheduler.go      Core interface, types, resource resolution, validation
registry.go       Thread-safe active scheduler singleton + debugMode
error.go          Structured error types
gpu.go            GPU compatibility database, tier matching, MIG profiles
slurm.go          SLURM implementation
pbs.go            PBS/Torque implementation
lsf.go            LSF implementation
htcondor.go       HTCondor implementation
```

## Interface

```go
type Scheduler interface {
    IsAvailable() bool
    ReadScriptSpecs(scriptPath string) (*ScriptSpecs, error)
    CreateScriptWithSpec(spec *JobSpec, outputDir string) (string, error)
    Submit(scriptPath string, dependencyJobIDs []string) (string, error)
    GetClusterInfo() (*ClusterInfo, error)
    GetInfo() *SchedulerInfo
    GetJobResources() *ResourceSpec // from scheduler env vars; zero fields = not set
}
```

## Resource Spec Priority Chain

Resources are resolved in three layers (lowest → highest priority):

```
defaults  →  script directives (HasDirectives=true only)  →  job resources (non-zero fields)
```

- **defaults**: caller-supplied baseline (`GetSpecDefaults()` for run; `buildDefaults` for build)
- **script directives**: parsed `#SBATCH`/`#PBS`/`#BSUB` values — only applied when `HasDirectives=true`
- **job resources**: live allocation from scheduler env vars (read by `GetJobResources()`) — only non-zero fields override

### Why `HasDirectives`, not `Spec != nil`

When a script has no directives, the parser fills `Spec` with `GetSpecDefaults()`. Checking only `Spec != nil` would silently replace the caller's baseline with scheduler defaults. Checking `HasDirectives` ensures no-directive scripts fall back to the caller's baseline correctly.

### Functions

```go
// Run system: uses GetSpecDefaults() as baseline
ResolveResourceSpec(jobResources, scriptSpecs) *ResourceSpec

// Build system (or any caller with custom baseline)
ResolveResourceSpecFrom(defaults, jobResources, scriptSpecs) *ResourceSpec

// In-place merge: overwrites rs fields with non-zero/non-nil fields from other
(*ResourceSpec).Override(other *ResourceSpec)

// Scheduler env vars → *ResourceSpec; starts from zero (no defaults)
sched.GetJobResources() *ResourceSpec
```

## Key Types

**`ScriptSpecs`** — parsed job specification:
- `HasDirectives bool` — `false` when no `#SBATCH`/`#PBS`/`#BSUB` lines exist
- `Spec *ResourceSpec` — resource geometry; `nil` in passthrough mode
- `Control RuntimeConfig` — job name, stdout/stderr, email, partition
- `RawFlags []string` — all directives verbatim
- `RemainingFlags []string` — directives not consumed by Spec or Control

**`ResourceSpec`** — compute geometry: `Nodes`, `TasksPerNode`, `CpusPerTask`, `MemPerNodeMB`, `Time`, `Gpu`, `Exclusive`

**`JobSpec`** — submission input: `Name`, `Command`, `Specs`, `DepJobIDs`, `Metadata`, `OverrideOutput`

## Script Parsing

Two-stage pipeline:
1. **Stage 1 (critical):** `parseRuntimeConfig` — job control (name, I/O, email). Failure stops parsing.
2. **Stage 2 (best-effort):** `parseResourceSpec` — compute geometry. Failure → **passthrough mode** (`Spec=nil`).

`ReadScriptSpecsFromPath` always returns non-nil `*ScriptSpecs`:

| State | `HasDirectives` | `Spec` | Meaning |
|---|---|---|---|
| No directives | `false` | non-nil (caller defaults) | No scheduler lines found |
| Normal | `true` | non-nil (parsed values) | Parsed successfully |
| Passthrough | `true` | `nil` | Directives found, resource parse failed |

Passthrough is triggered by: invalid directive values, unresolvable relationships (e.g. `--ntasks` not divisible by `--nodes`), or blacklisted topology flags (SLURM only).

## Configurable Scheduler Defaults

Applied when a script omits resource directives. Set via `SetSpecDefaults()` at CLI init (wired from `config.Global.Scheduler`).

| Key | Default | Description |
|-----|---------|-------------|
| `scheduler.nodes` | 1 | Node count |
| `scheduler.tasks_per_node` | 1 | Tasks per node |
| `scheduler.ncpus_per_task` | 2 | CPUs per task |
| `scheduler.mem_mb_per_node` | 8192 | Memory per node (MB) |
| `scheduler.time` | `4h` | Time limit |

```yaml
# ~/.config/condatainer/config.yaml
scheduler:
  ncpus_per_task: 8
  mem_mb_per_node: 16384
  time: "8h"
```

## Convenience Functions

```go
// Detection
scheduler.Init("")                  // auto-detect and set active scheduler
scheduler.InitIfAvailable("")       // non-failing
scheduler.ActiveScheduler()         // get current scheduler

// Script helpers
ReadScriptSpecsFromPath(path)       // always non-nil
ParseScriptAny(path)                // cross-scheduler parsing
HasSchedulerSpecs(specs)            // HasDirectives=true
IsPassthrough(specs)                // HasDirectives=true && Spec=nil

// Submission
SubmitJob(sched, scriptPath, deps)  // single job
SubmitJobs(sched, jobs, outputDir)  // multi-job with dependencies

// Validation
ValidateAndConvertSpecs(specs)      // validates against cluster limits (no auto-adjust)
ValidateSpecs(specs, limits)        // validate against single partition
ValidateGpuAvailability(gpu, info)  // returns GpuValidationError with suggestions
```

## Normalized Environment Variables

Generated job scripts export these variables regardless of scheduler:

| Variable | Description |
|----------|-------------|
| `NNODES` | Number of nodes |
| `NTASKS` | Total tasks |
| `NTASKS_PER_NODE` | Tasks per node |
| `NCPUS` | CPUs per node (`CpusPerTask × TasksPerNode`) |
| `NCPUS_PER_TASK` | CPUs per task |
| `MEM` / `MEM_MB` | Memory in MB |
| `MEM_GB` | Memory in GB |

### Scheduler-Native Variables (Read-Only)

Each scheduler also sets its own runtime environment variables when inside a job:

| Variable | SLURM | PBS | LSF | HTCondor |
|----------|-------|-----|-----|----------|
| Job ID | `SLURM_JOB_ID` | `PBS_JOBID` | `LSB_JOBID` | `_CONDOR_JOB_AD` |
| CPUs/Task | `SLURM_CPUS_PER_TASK` | `PBS_NCPUS` | `LSB_DJOB_NUMPROC` | `_CONDOR_REQUEST_CPUS` |
| Memory | `SLURM_MEM_PER_NODE` | `PBS_VMEM` | `LSB_MAX_MEM_RUSAGE` | `_CONDOR_REQUEST_MEMORY` |
| GPUs | `CUDA_VISIBLE_DEVICES` | `CUDA_VISIBLE_DEVICES` | `CUDA_VISIBLE_DEVICES` | `CUDA_VISIBLE_DEVICES` |
| Nodes | `SLURM_JOB_NUM_NODES` | `PBS_NUM_NODES` | N/A | N/A |
| Tasks/Node | `SLURM_NTASKS_PER_NODE` | `PBS_NP` ÷ nodes | N/A | N/A |

## Error Types

| Type | Returned by | Key fields |
|------|-------------|------------|
| `ValidationError` | `ValidateSpecs`, `ValidateAndConvertSpecs` | `Field`, `Requested`, `Limit`, `Partition` |
| `GpuValidationError` | `ValidateGpuAvailability` | `Suggestions` |
| `ParseError` | script parsing | — |
| `SubmissionError` | `Submit` | — |
| `ClusterError` | `GetClusterInfo` | — |

Check with: `scheduler.IsValidationError(err)`, `scheduler.IsGpuValidationError(err)`, etc.

## Scheduler-Specific Notes

### SLURM
- GPU: `--gres=gpu:type:count` > `--gpus-per-node` > `--gpus`÷nodes > `--gpus-per-task`×tasks
- Memory: `--mem` (per-node), `--mem-per-cpu` (×cpus), `--mem-per-gpu` (×gpu count)
- Blacklisted flags (trigger passthrough): `--gpus-per-socket`, `--sockets-per-node`, `--cores-per-socket`, `--threads-per-core`, `--ntasks-per-socket`, `--ntasks-per-core`, `--distribution`
- Cluster info from `sinfo` + `scontrol`

### PBS
- Resources: `select=N:ncpus=M:mpiprocs=P:mem=Xgb:ngpus=G` or `nodes=N:ppn=M`
- Cluster info from `pbsnodes -a`

### LSF
- CPU: `-n`, memory: `-M` (KB), nodes: `-R "span[hosts=N]"`
- GPU: `-gpu "num=N:type=T"` or `-R "rusage[ngpus_physical=N]"`
- Custom rusage params beyond mem/ngpus are lost during parsing
- Cluster info from `bhosts -w`

### HTCondor
- Parses native `.sub` files; `executable` field → `ScriptSpecs.ScriptPath`
- Resource keys: `request_cpus`, `request_memory`, `request_gpus`, `+MaxRuntime` (seconds)
- Generates `.sub` + `.sh` wrapper; always single-node; no dependency support

## Single-Node Enforcement

`CreateScriptWithSpec` normalizes to single-node by default:

| Scheduler | Mechanism |
|-----------|-----------|
| SLURM | `--nodes=1 --ntasks=1` |
| PBS | `select=1:ncpus=N` |
| LSF | `span[hosts=1]` + `-n N` |
| HTCondor | inherently single-node |

Multi-node: set `Nodes > 1` in `ScriptSpecs`.

## Dependency Formats

| Scheduler | Format |
|-----------|--------|
| SLURM | `--dependency=afterok:ID1,ID2` |
| PBS | `-W depend=afterok:ID1:ID2` |
| LSF | `-w "done(ID1) && done(ID2)"` |
| HTCondor | not supported (requires DAGMan) |
