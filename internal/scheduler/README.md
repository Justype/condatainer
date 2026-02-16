# scheduler

HPC scheduler abstraction layer providing a unified interface for job submission, script parsing, and resource management across SLURM, PBS/Torque, LSF, and HTCondor.

## Supported Schedulers

| Scheduler | Binary | Directive | Status |
|-----------|--------|-----------|--------|
| SLURM | `sbatch` | `#SBATCH` | Stable |
| PBS/Torque | `qsub` | `#PBS` | Experimental |
| LSF | `bsub` | `#BSUB` | Experimental |
| HTCondor | `condor_submit` | `#CONDOR` | Experimental |

## Architecture

```
scheduler.go      Core interface, types, detection, cross-scheduler parsing
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
- `JobName`, `Ncpus`, `Ntasks`, `Nodes`, `MemMB`, `Time`
- `Gpu *GpuSpec`, email settings, `RawFlags []string`

**JobSpec** - Job submission input:
- `Name`, `Command`, `Specs *ScriptSpecs`, `DepJobIDs`, `Metadata`

**JobResources** - Runtime-allocated resources (from env vars, pointer fields):
- `Ncpus`, `Ntasks`, `Nodes`, `MemMB`, `Ngpus`

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
        Ncpus: 8,
        MemMB: 16384,
        Time:  3600,
        Gpu:   &scheduler.GpuSpec{Count: 1, Type: "a100"},
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
- Parses `--nodes`, `--ntasks`, `--ntasks-per-node`, `--cpus-per-task`
- GPU via `--gres=gpu:type:count` or `--gpus=type:count`; supports MIG profiles
- Cluster info from `sinfo` and `scontrol`

### PBS
- Resource formats: `select=N:ncpus=M:mpiprocs=P:mem=Xgb:ngpus=G` or `nodes=N:ppn=M`
- `Ntasks` derived from `nodes * mpiprocs`
- Cluster info from `pbsnodes -a`

### LSF
- CPU count via `-n`, memory via `-M` (KB), node count via `-R "span[hosts=N]"`
- GPU via `-gpu "num=N:type=T"`
- Mixed `-R` flags (e.g., `span[...] rusage[...]`) are split: span is regenerated, rusage is preserved
- Cluster info from `bhosts -w`

### HTCondor
- Uses `key = value` directive format (e.g., `request_cpus = 8`)
- Generates two files: `.sub` (submit description) + `.sh` (wrapper script)
- Time via `+MaxRuntime` (seconds)
- Always single-node (vanilla universe)
- No native job dependency support

## Error Types

- `ValidationError` - Resource limits exceeded
- `ParseError` - Directive parsing failure
- `SubmissionError` - Job submission failure
- `ClusterError` - Cluster info query failure
- `DependencyError` - Missing job dependencies
- `ScriptCreationError` - Script file creation failure
- `GpuValidationError` - GPU unavailable (includes suggestions)

Check errors with: `scheduler.IsValidationError(err)`, `scheduler.IsParseError(err)`, etc.

## Environment Variables

Each scheduler reads runtime resources from environment variables set by the scheduler when inside a job:

| Variable | SLURM | PBS | LSF | HTCondor |
|----------|-------|-----|-----|----------|
| Job ID | `SLURM_JOB_ID` | `PBS_JOBID` | `LSB_JOBID` | `_CONDOR_JOB_AD` |
| CPUs | `SLURM_CPUS_PER_TASK` | `PBS_NCPUS` | `LSB_DJOB_NUMPROC` | `_CONDOR_REQUEST_CPUS` |
| Memory | `SLURM_MEM_PER_NODE` | `PBS_VMEM` | `LSB_MAX_MEM_RUSAGE` | `_CONDOR_REQUEST_MEMORY` |
| GPUs | `CUDA_VISIBLE_DEVICES` | `PBS_NGPUS` | `CUDA_VISIBLE_DEVICES` | `CUDA_VISIBLE_DEVICES` |
| Nodes | `SLURM_JOB_NUM_NODES` | `PBS_NUM_NODES` | `LSB_HOSTS` | N/A |
| Tasks | `SLURM_NTASKS` | `PBS_NP` | N/A | N/A |
