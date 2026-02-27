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
scheduler.go        Core interface, types, resource resolution, validation
registry.go         Thread-safe active scheduler singleton + debugMode
error.go            Structured error types
gpu.go              GPU compatibility database, tier matching, MIG profiles
script_helpers.go   Shared parsing/formatting helpers (memory, time, flags, job header/footer)
slurm.go            SLURM implementation
pbs.go              PBS/Torque implementation
lsf.go              LSF implementation
htcondor.go         HTCondor implementation
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

## Key Types

**`ScriptSpecs`** — parsed job specification:
- `HasDirectives bool` — `false` when no `#SBATCH`/`#PBS`/`#BSUB` lines exist
- `Spec *ResourceSpec` — resource geometry; `nil` in passthrough mode
- `Control RuntimeConfig` — job name, stdout/stderr, email, partition

**`ResourceSpec`** — compute geometry: `Nodes`, `TasksPerNode`, `CpusPerTask`, `MemPerNodeMB`, `Time`, `Gpu`, `Exclusive`

**`JobSpec`** — submission input: `Name`, `Command`, `Specs`, `DepJobIDs`, `Metadata`, `Array *ArraySpec`

**`ArraySpec`** — array job descriptor:
- `InputFile string` — absolute path to the input list (one entry per line)
- `Count int` — number of non-empty lines (= number of subjobs)
- `Limit int` — max concurrently running subjobs (0 = unlimited)
- `ArgCount int` — number of shell-split tokens per line (validated uniform across all lines)
- `SampleArgs []string` — tokens from the first line, used for dry-run display
- `BlankLines []int` — 1-based line numbers of all blank lines; warned on dry-run, blocks submission

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
// Detection
scheduler.Init("")                  // auto-detect and set active scheduler
scheduler.InitIfAvailable("")       // non-failing
scheduler.ActiveScheduler()         // get current scheduler

// Script parsing
ReadScriptSpecsFromPath(path)       // always non-nil
HasSchedulerSpecs(specs)            // HasDirectives=true
IsPassthrough(specs)                // HasDirectives=true && Spec=nil

// Submission
SubmitJob(sched, scriptPath, deps)  // single job
SubmitJobs(sched, jobs, outputDir)  // multi-job with dependencies

// Validation
ValidateAndConvertSpecs(specs)      // validates against cluster limits
ValidateGpuAvailability(gpu, nodes, info) // returns GpuValidationError with suggestions
```

## Configurable Scheduler Defaults

Applied when a script omits resource directives. Set via `SetSpecDefaults()` at CLI init.

```yaml
# ~/.config/condatainer/config.yaml
scheduler:
  nodes: 1
  tasks_per_node: 1
  ncpus_per_task: 2
  mem_per_node_mb: 8192
  time: "4h"
```

## Normalized Environment Variables

Generated job scripts export these regardless of scheduler:
`NNODES`, `NTASKS`, `NTASKS_PER_NODE`, `NCPUS`, `NCPUS_PER_TASK`, `MEM_MB`, `MEM_GB`

## Error Types

- `ValidationError` - `IsValidationError(err)` — fields: `Field`, `Requested`, `Limit`, `Partition`
- `GpuValidationError` - `IsGpuValidationError(err)` — has `Suggestions`
- `ParseError`, `SubmissionError`, `ClusterError` — checked with `Is*` helpers

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

## Array Jobs

Pass `JobSpec.Array = &ArraySpec{...}` to `CreateScriptWithSpec` to generate an array job.

**How it works (SLURM / PBS / LSF):**

1. Scheduler directive silences per-subjob output (`--output=/dev/null`).
2. `writeArrayBlock` is injected after directives. It:
   - Extracts the current subjob's input line via `sed -n "${IDX}p"` into `ARRAY_ARGS`.
   - Exports `ARRAY_ARGS` (word-splits into positional args when used unquoted in the command).
   - Redirects per-subjob output to `{logDir}/{jobName}_{paddedIdx}_{tag}.log` (combined) or `.out`/`.err` (separated, when the original stderr differs from stdout). `{tag}` is `ARRAY_ARGS` with all non-alphanumeric characters replaced by `_`, consecutive underscores squeezed, capped at 20 characters.
3. The user command receives `$ARRAY_ARGS` prepended as positional arguments.

**HTCondor:**

Uses `queue array_args from {inputFile}` — HTCondor sets `ARRAY_ARGS` via an `environment` directive for each subjob. Per-subjob output goes to `{logDir}/{jobName}_$(Process)_$(array_args).log`.

**Scheduler directives:**

| Scheduler | Directive | Range format |
|-----------|-----------|--------------|
| SLURM | `#SBATCH --array=1-N[%limit]` | 1-based, optional concurrency cap |
| PBS | `#PBS -J 1-N` + `qdel` limit | 1-based |
| LSF | `#BSUB -J name[1-N]` | 1-based |
| HTCondor | `queue array_args from file` | item-based |

## Dependency Formats

| Scheduler | Format |
|-----------|--------|
| SLURM | `--dependency=afterok:ID1,ID2` |
| PBS | `-W depend=afterok:ID1:ID2` |
| LSF | `-w "done(ID1) && done(ID2)"` |
| HTCondor | not supported (requires DAGMan) |
