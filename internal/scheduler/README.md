# scheduler

HPC scheduler abstraction layer providing a unified interface for job submission, script parsing, and resource management across SLURM, PBS/Torque, LSF, and HTCondor.

## Supported Schedulers

| Scheduler | Binary | Directive | Status |
|-----------|--------|-----------|--------|
| SLURM | `sbatch` | `#SBATCH` | Stable |
| PBS/Torque | `qsub` | `#PBS` | Experimental |
| LSF | `bsub` | `#BSUB` | Experimental |
| HTCondor | `condor_submit` | Native `.sub` format | Experimental |

## Architecture

```
scheduler.go      Core interface, types, detection, cross-scheduler parsing,
                  ValidateAndConvertSpecs (resource auto-adjustment)
registry.go       Thread-safe global active scheduler singleton
error.go          Structured error types (validation, parse, submission, etc.)
gpu.go            GPU compatibility database, tier matching, MIG profiles
slurm.go          SLURM implementation
pbs.go            PBS/Torque implementation
lsf.go            LSF implementation
htcondor.go       HTCondor implementation
```

## Interface

All schedulers implement the `Scheduler` interface:

```go
type Scheduler interface {
    IsAvailable() bool
    ReadScriptSpecs(scriptPath string) (*ScriptSpecs, error)
    CreateScriptWithSpec(spec *JobSpec, outputDir string) (string, error)
    Submit(scriptPath string, dependencyJobIDs []string) (string, error)
    GetClusterInfo() (*ClusterInfo, error)
    GetInfo() *SchedulerInfo
    GetJobResources() *JobResources
}
```

## Key Types

**ScriptSpecs** - Parsed job specifications (scheduler-agnostic):
- `ScriptPath` - Absolute path of the parsed script (for HTCondor: the `executable` from `.sub` file)
- `Spec *ResourceSpec` - Compute geometry: `Nodes`, `TasksPerNode`, `CpusPerTask`, `MemPerNodeMB`, `Time`, `Gpu`, `Exclusive`
- `Control RuntimeConfig` - Job control: `JobName`, `Stdout`, `Stderr`, email settings, `Partition`
- `RawFlags []string` - Immutable audit log of all directives
- `RemainingFlags []string` - Directives not absorbed by Spec or Control

**JobSpec** - Job submission input:
- `Name`, `Command`, `Specs *ScriptSpecs`, `DepJobIDs`, `Metadata`

**JobResources** - Runtime-allocated resources (from env vars, pointer fields):
- `Ncpus`, `Ntasks`, `Nodes`, `MemMB`, `Ngpus`

## Script Parsing Behavior

Parsing uses a two-stage pipeline:

1. **Stage 1 (critical):** `parseRuntimeConfig` — extracts job control settings (name, I/O, email). Failures stop parsing.
2. **Stage 2 (best-effort):** `parseResourceSpec` — extracts compute geometry (CPU, memory, GPU, time). Returns nil on any parse failure → **passthrough mode** (the script is forwarded to the scheduler unchanged).

**Passthrough mode** is triggered when:
- A known directive has an invalid value (e.g., non-integer CPU count)
- A directive value cannot be safely resolved (e.g., `--ntasks` not evenly divisible by `--nodes`)
- A blacklisted topology flag is encountered (SLURM only — see below)

`RawFlags` is an immutable audit log of ALL directives. `RemainingFlags` contains only directives not absorbed by either stage. This ensures that when generating scripts for a different scheduler, standard resources are regenerated using the target scheduler's native syntax.

## Key Functions

**Resource Validation & Adjustment:**
- `ValidateAndConvertSpecs(specs) → (err, cpuAdjusted, cpuMsg, gpuConverted, gpuMsg)` - Auto-adjusts CPUs/time, converts GPUs
- `ValidateSpecs(specs, limits) → error` - Validates specs against single partition limits
- `ValidateGpuAvailability(gpuSpec, clusterInfo) → error` - Checks GPU availability

**GPU Conversion:**
- `ConvertGpuSpec(requestedSpec, clusterInfo) → (*GpuSpec, error)` - Finds best compatible GPU
- `FindCompatibleGpu(requestedSpec, clusterInfo) → ([]*GpuConversionOption, error)` - Returns all compatible options

**Script Parsing:**
- `ReadScriptSpecsFromPath(path) → (*ScriptSpecs, error)` - Parses with active scheduler
- `ParseScriptAny(path) → (*ParsedScript, error)` - Cross-scheduler parsing

**Convenience Wrappers:**
- `SubmitJob(scheduler, scriptPath, deps) → (jobID, error)` - Single job submission
- `SubmitJobs(scheduler, jobs, outputDir) → (jobIDs, error)` - Multi-job with dependencies
- `ReadScript(scheduler, path)`, `ValidateScript(scheduler, path)`, `WriteScript(...)` - Helper aliases

## Configurable Defaults

`ReadScriptSpecs` applies configurable defaults when a script omits resource directives. Defaults are set via config (`scheduler.*` keys) and wired through `SetSpecDefaults()` at CLI init.

| Key | Default | Description |
|-----|---------|-------------|
| `scheduler.nodes` | 1 | Node count |
| `scheduler.tasks_per_node` | 1 | Tasks per node |
| `scheduler.ncpus_per_task` | 2 | CPUs per task |
| `scheduler.mem_mb_per_node` | 8192 | Memory per node in MB |
| `scheduler.time` | `4h` | Time limit |

```yaml
# ~/.config/condatainer/config.yaml
scheduler:
  nodes: 1
  tasks_per_node: 1
  ncpus_per_task: 8
  mem_mb_per_node: 16384
  time: "8h"
```

These are independent from build defaults (`build.default_cpus`), which control `$NCPUS` during image builds.

## Usage

### Detection and Initialization

```go
// Auto-detect and set active scheduler
schedType, err := scheduler.Init("")

// Or detect with a preferred binary
schedType, err := scheduler.Init("/usr/bin/sbatch")

// Non-failing detection
schedType := scheduler.InitIfAvailable("")

// Access the active scheduler
sched := scheduler.ActiveScheduler()
```

### Parsing Scripts

```go
// Parse using detected scheduler
specs, err := sched.ReadScriptSpecs("job.sh")

// Cross-scheduler parsing (tries all parsers)
parsed := scheduler.ParseScriptAny("job.sh")
```

### Generating Scripts

```go
jobSpec := &scheduler.JobSpec{
    Name:    "my_job",
    Command: "python train.py",
    Specs: &scheduler.ScriptSpecs{
        Spec: &scheduler.ResourceSpec{
            CpusPerTask:  8,
            MemPerNodeMB: 16384,
            Time:         time.Hour,
            Gpu:          &scheduler.GpuSpec{Count: 1, Type: "a100"},
        },
        Control: scheduler.RuntimeConfig{JobName: "my_job"},
    },
}
scriptPath, err := sched.CreateScriptWithSpec(jobSpec, "/output/dir")
```

### Job Submission

```go
// Single job
jobID, err := scheduler.SubmitJob(sched, scriptPath, nil)

// Multiple jobs with dependencies
jobs := []*scheduler.JobSpec{job1, job2, job3}
jobIDs, err := scheduler.SubmitJobs(sched, jobs, outputDir)
```

### Validation and Resource Auto-Adjustment

`ValidateAndConvertSpecs` validates job specs against cluster limits and auto-adjusts resources when possible. Called automatically by `condatainer run` and build submissions.

| Resource | Auto-Adjust | Notes |
|----------|-------------|-------|
| **CPUs** | ✓ Reduced | Time scaled proportionally; fails if adjusted time > limit |
| **Memory** | ✗ Critical | Always fails if exceeded |
| **Time** | ✗ Critical | Always fails if exceeded |
| **GPUs** | ✓ Converted | Upgrades allowed; downgrades blocked if time limit set |

```go
// Auto-adjusts specs in-place
validationErr, cpuAdjusted, cpuMsg, gpuConverted, gpuMsg := scheduler.ValidateAndConvertSpecs(specs)

// Example warnings:
// "CPUs: 32 → 16 (reduced to fit limit); Time: 4h → 8h (adjusted proportionally)"
// "GPU: v100 → a100 (upgrade)"
```

**Key behaviors:**
- CPU reduction: Assumes linear scaling, adjusts time proportionally
- GPU upgrades: Always safe (faster = better runtime)
- GPU downgrades: Blocked when time limit exists (unpredictable performance)
- Validation skipped if cluster info unavailable

## Single-Node Enforcement

By default, all schedulers generate scripts that run on a **single node** (`Nodes=1`, `Ntasks=1`). `CreateScriptWithSpec` normalizes these defaults and regenerates resource directives explicitly, skipping any raw flags that would conflict:

| Scheduler | Mechanism |
|-----------|-----------|
| SLURM | `--nodes=1 --ntasks=1` |
| PBS | `select=1:ncpus=N` |
| LSF | `span[hosts=1]` + `-n N` |
| HTCondor | Inherently single-node (vanilla universe) |

Multi-node is supported by setting `Nodes > 1` in `ScriptSpecs`.

## Dependency Formats

| Scheduler | Format |
|-----------|--------|
| SLURM | `--dependency=afterok:ID1,ID2` |
| PBS | `-W depend=afterok:ID1:ID2` |
| LSF | `-w "done(ID1) && done(ID2)"` |
| HTCondor | Not supported (requires DAGMan) |

## Scheduler-Specific Notes

### SLURM
- Parses `--nodes`, `--ntasks`, `--ntasks-per-node`, `--cpus-per-task`, `--cpus-per-gpu`
- GPU via `--gres=gpu:type:count`, `--gpus-per-node`, `--gpus` (total), or `--gpus-per-task`; supports MIG profiles
  - Resolution priority: `--gres=gpu:` > `--gpus-per-node` > `--gpus` (÷ nodes) > `--gpus-per-task` (× tasks/node)
  - `--gpus` total must divide evenly by `--nodes`; otherwise triggers passthrough
- Memory via `--mem` (per-node), `--mem-per-cpu` (× CpusPerTask × TasksPerNode), or `--mem-per-gpu` (× Gpu.Count)
- Resource spec uses a **two-phase approach**: Phase 1 scans all directives into temp vars; Phase 2 resolves in dependency order (Nodes/Tasks → GPU → CPU → Memory → Time), ensuring derived values are computed correctly
- `--ntasks` must divide evenly by `--nodes`; otherwise triggers passthrough
- **Blacklisted topology flags** (immediately trigger passthrough — cannot be reliably translated):
  `--gpus-per-socket`, `--sockets-per-node`, `--cores-per-socket`, `--threads-per-core`,
  `--ntasks-per-socket`, `--ntasks-per-core`, `--distribution`
- Cluster info from `sinfo` and `scontrol`

### PBS
- Resource formats: `select=N:ncpus=M:mpiprocs=P:mem=Xgb:ngpus=G` or `nodes=N:ppn=M`
- `Ntasks` derived from `nodes * mpiprocs`
- Invalid values for `ncpus`, `mem`, `walltime`, `ngpus`, `gpus`, `select`, `nodes`, `ppn`, `mpiprocs` trigger passthrough
- Cluster info from `pbsnodes -a`

### LSF
- CPU count via `-n`, memory via `-M` (KB), node count via `-R "span[hosts=N]"`
- GPU via `-gpu "num=N:type=T"` or `-R "rusage[ngpus_physical=N]"`
- `-R` resource flags are parsed: `span[hosts=N]` sets `Nodes`, `rusage[mem=X]` sets `MemMB`, `rusage[ngpus*=N]` sets GPU count
- Script generation always regenerates `-R "span[hosts=N]"` from `specs.Nodes`
- **Note**: Custom rusage parameters in `-R` flags (beyond mem/ngpus) are currently lost during parsing; use separate flags for custom resources
- Cluster info from `bhosts -w`

### HTCondor
- Parses native `.sub` submit files directly
- `executable = script.sh` in the `.sub` file is extracted as `ScriptSpecs.ScriptPath`
- Comments (`#`), empty lines, and structural keywords (`universe`, `executable`, `transfer_executable`, `queue`, `arguments`) are skipped during parsing
- Resource keys: `request_cpus`, `request_memory`, `request_gpus`, `+MaxRuntime` (seconds)
- Generates two files: `.sub` (submit description) + `.sh` (wrapper script)
- Always single-node (vanilla universe)
- No native job dependency support

## Error Types

- `ValidationError` - Resource limits exceeded (CPU, memory, time, GPU count)
  - Returned by `ValidateSpecs()` and `ValidateAndConvertSpecs()`
  - Contains `Field`, `Requested`, `Limit`, `Partition` details
  - Examples: Memory 128GB > limit 64GB, Time 48h > limit 24h
- `GpuValidationError` - GPU type unavailable or incompatible
  - Returned by `ValidateGpuAvailability()` when no conversion possible
  - Includes `Suggestions` list with available GPU options
- `ParseError` - Directive parsing failure
- `SubmissionError` - Job submission failure
- `ClusterError` - Cluster info query failure
- `DependencyError` - Missing job dependencies
- `ScriptCreationError` - Script file creation failure

Check errors with: `scheduler.IsValidationError(err)`, `scheduler.IsGpuValidationError(err)`, etc.

## Environment Variables

### Normalized Variables (Exported by CondaTainer)

All generated job scripts export **normalized environment variables** for use in build scripts, regardless of scheduler:

| Variable | Description | Example |
|----------|-------------|---------|
| `NNODES` | Number of nodes | `1` |
| `NTASKS` | Total number of tasks | `1` |
| `NCPUS` | CPUs per task | `8` |
| `MEM` | Memory in MB | `16384` |
| `MEM_MB` | Memory in MB (alias) | `16384` |
| `MEM_GB` | Memory in GB | `16` |

These variables are set by CondaTainer's script generation and are available inside all build scripts for resource-aware compilation (e.g., `make -j $NCPUS`).

### Scheduler-Native Variables (Read-Only)

Each scheduler also sets its own runtime environment variables when inside a job:

| Variable | SLURM | PBS | LSF | HTCondor |
|----------|-------|-----|-----|----------|
| Job ID | `SLURM_JOB_ID` | `PBS_JOBID` | `LSB_JOBID` | `_CONDOR_JOB_AD` |
| CPUs | `SLURM_CPUS_PER_TASK` | `PBS_NCPUS` | `LSB_DJOB_NUMPROC` | `_CONDOR_REQUEST_CPUS` |
| Memory | `SLURM_MEM_PER_NODE` | `PBS_VMEM` | `LSB_MAX_MEM_RUSAGE` | `_CONDOR_REQUEST_MEMORY` |
| GPUs | `CUDA_VISIBLE_DEVICES` | `PBS_NGPUS` | `CUDA_VISIBLE_DEVICES` | `CUDA_VISIBLE_DEVICES` |
| Nodes | `SLURM_JOB_NUM_NODES` | `PBS_NUM_NODES` | `LSB_HOSTS` | N/A |
| Tasks | `SLURM_NTASKS` | `PBS_NP` | N/A | N/A |
