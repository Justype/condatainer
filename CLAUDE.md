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
- `config/` - Multi-level config (flags > env > user > extra-root > app-root > system > defaults), data directory search
- `container/` - Container setup pipeline: overlay resolution, bind dedup, env collection, GPU detection
- `exec/` - Ephemeral container execution
- `instance/` - Persistent instance management (start/exec/stop) with state persistence
- `overlay/` - Overlay image CRUD (ext3/SquashFS), resize, chown, locking
- `scheduler/` - HPC scheduler abstraction (SLURM, PBS, LSF, HTCondor); auto-detection, directive parsing, cross-scheduler translation
- `utils/` - Console output (`Print*`), file ops, downloads, script parsing

## Build Scripts

Build scripts live in the [`cnt-scripts`](https://github.com/Justype/cnt-scripts) repo (fetched remotely via `scripts_link` config, or auto-detected from a local `cnt-scripts/` clone next to the binary). Three categories:

- **OS**: `<distro>/<name>` (e.g., `ubuntu24/igv`) — Apptainer definition files for distro system tools
- **Apps**: `<name>/<version>` (e.g., `cellranger/9.0.1`) — Apps that not available as conda packages, or specific versions not in conda
- **Data**: `<assembly|project>/<datatype>/<version>` (e.g., `grch38/star/2.7.11b/gencode47-101`) — any data, including genome reference indexes

Must define an `install()` function. Available vars: `$NCPUS`, `$target_dir`, `$tmp_dir`, `$app_name`, `$version`.

Metadata headers: `#DEP:name/version` or `#DEP:name/version>=min` (deps; preferred version is implicit upper bound, so valid range is `[min, version]`), `#SBATCH`/`#PBS`/`#BSUB` (scheduler job params), `#ENV:VAR=$app_root` (env vars), `#INTERACTIVE:prompt` (user input).

Overlays are stored as `.sqf` (SquashFS, read-only) or `.img` (ext3, writable).

## Data Directory Search Order

1. `extra_image_dirs` / `extra_build_dirs` / `extra_helper_dirs` (config keys)
2. `CNT_EXTRA_ROOT` (group/lab root, env only)
3. `CNT_ROOT` / `<install-dir>/` (app-root, auto-detected)
4. `$SCRATCH/condatainer/`
5. `~/.local/share/condatainer/`

Each contains `images/`, `build-scripts/`, `helper-scripts/`. Writes go to first writable dir.

## Helper Scripts

Bash scripts in [`cnt-scripts/helpers/`](https://github.com/Justype/cnt-scripts) launch interactive services inside CondaTainer on HPC. Modes: `headless/` (direct) and `<scheduler>/` (submit + SSH tunnel). See `helpers/README.md` for details.

## File Locking

`exec`/`run` hold `LOCK_SH` on `.sqf`/`.sif` files during execution (`.img` skipped — Apptainer flocks those itself); `remove` and `build --update` probe `LOCK_EX` before modifying. See `internal/overlay/lock.go`.

## Coding Rules

- When editing a function's behaviour, also update its doc comment to match.
- When changing UX (flags, output format, command behaviour), also update `docs/manuals/condatainer.md` and any relevant `docs/` pages.

## Key Patterns

- Global config singleton: `config.Global`
- Error types: `ApptainerError`, `ValidationError` with structured fields
- Console output: `utils.PrintMessage`, `PrintWarning`, `PrintError`, `PrintDebug`
- Always use absolute paths; `config.Get*Dir()` for standard locations
- File permissions: `utils.PermFile` (0664), `utils.PermDir` (0775), `utils.PermExec` (0775)
- Prefer `utils.CreateFileWritable(path)` for runtime file creation so new files consistently use `utils.PermFile`.
