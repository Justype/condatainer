# CLAUDE.md

## Project Overview

CondaTainer is an HPC-oriented CLI tool that packs Conda environments into single image files (SquashFS/OverlayFS) via Apptainer and Micromamba to avoid inode quota issues on HPC systems.

## Commands

```bash
make build          # Build binary to bin/condatainer_go
make test           # Run all tests (or: go test ./...)
go test -v ./internal/build -run TestParseScriptMetadata  # Single test
go test -v ./internal/scheduler/...                       # Package tests
```

## Architecture

- **CLI**: spf13/cobra + spf13/viper. Entry: `main.go` → `cmd.Execute()`
- **cmd/**: One file per subcommand, delegates to internal packages

**internal/** packages:
- `apptainer/` - Apptainer binary wrapper (exec, build, instance management)
- `build/` - Build system: resolves name/version → Conda, Script, or Def build type; dependency graphs; remote script fetching
- `config/` - Multi-level config (flags > env > user > portable > system > defaults), data directory search
- `container/` - Container setup pipeline: overlay resolution, bind dedup, env collection, GPU detection
- `exec/` - Ephemeral container execution
- `instance/` - Persistent instance management (start/exec/stop) with state persistence
- `overlay/` - Overlay image CRUD (ext3/SquashFS), resize, chown, locking
- `scheduler/` - HPC scheduler abstraction (SLURM, PBS, LSF, HTCondor); auto-detection, directive parsing, cross-scheduler translation
- `utils/` - Console output (`Print*`), file ops, downloads, script parsing

## Build Scripts

Located in `build-scripts/` with naming:
- Apps: `name/version` (e.g., `cellranger/9.0.1`)
- References: `assembly/data-type/version` (e.g., `grch38/star/2.7.11b/gencode47-101`)

Must define an `install()` function. Available vars: `$NCPUS`, `$target_dir`, `$tmp_dir`, `$app_name`, `$version`.

Metadata headers: `#DEP:name/version` (deps), `#SBATCH`/`#PBS`/`#BSUB`/`#CONDOR` (scheduler job params), `#ENV:VAR=$app_root` (env vars), `#INTERACTIVE:prompt` (user input).

Overlays are stored as `.sqf` (SquashFS, read-only) or `.img` (ext3, writable).

## Data Directory Search Order

1. `CONDATAINER_EXTRA_BASE_DIRS` / config `extra_base_dirs`
2. `<install-dir>/` (portable, auto-detected)
3. `$SCRATCH/condatainer/`
4. `~/.local/share/condatainer/`

Each contains `images/`, `build-scripts/`, `helper-scripts/`. Writes go to first writable dir.

## Helper Scripts

Bash scripts in `helpers/` launch interactive services inside CondaTainer on HPC. Modes: `headless/` (direct) and `<scheduler>/` (submit + SSH tunnel). See `helpers/README.md` for details.

Shared library `.common.sh` provides: config management (`config_init/load/require`), port helpers (`choose_port`, `validate_port`), overlay checks (`check_overlay_integrity`, `check_and_install_overlays`), job state (`read_job_state`, `wait_for_job`), reuse mode (`handle_reuse_mode`), and display (`spec_line`, `print_specs`, `countdown`).

Scheduler script flow: config_init → config_load → read_job_state → getopts (sets `_ARG_*` flags) → handle_reuse_mode → port resolution → print_specs → sanity checks → submit → wait_for_job → connect.

## Key Patterns

- Global config singleton: `config.Global`
- Error types: `ApptainerError`, `ValidationError` with structured fields
- Console output: `utils.PrintMessage`, `PrintWarning`, `PrintError`, `PrintDebug`
- Always use absolute paths; `config.Get*Dir()` for standard locations
- File permissions: `utils.PermFile` (0664), `utils.PermDir` (0775), `utils.PermExec` (0775)
