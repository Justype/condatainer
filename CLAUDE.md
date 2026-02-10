# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

CondaTainer is an HPC-oriented CLI tool that streamlines environment and data management using Apptainer, SquashFS, OverlayFS, and Micromamba. It packs Conda environments into single image files to avoid inode quota issues on HPC systems.

## Build and Development Commands

```bash
make build          # Build binary to bin/condatainer_go
make test           # Run all tests (or: go test ./...)
make clean          # Remove build artifacts
make deps           # Download Go dependencies

# Run a single test
go test -v ./internal/build -run TestParseScriptMetadata

# Run tests for a specific package
go test -v ./internal/scheduler/...
```

## Architecture

### CLI Framework
- Uses spf13/cobra for command parsing, spf13/viper for configuration
- Entry point: `main.go` → `cmd.Execute()`
- Root command in `cmd/root.go` handles global flags and initialization sequence

### Package Structure

**cmd/** - CLI command implementations
- Each file defines a subcommand (exec, create, overlay, config, etc.)
- Commands call into internal packages for business logic

**internal/** - Core implementation packages:
- `apptainer/` - Apptainer binary wrapper, exec/build operations, GPU detection
- `build/` - Build system for overlays (BuildObject interface with Conda, Def, Script types)
- `config/` - Configuration loading (viper), path resolution, data directories
- `exec/` - Container execution with overlay mounting, environment setup, bind paths
- `overlay/` - Overlay image creation (ext3/ext4), resize, chown, inspection
- `scheduler/` - HPC scheduler abstraction (SLURM, PBS support), job submission
- `utils/` - Console output, file operations, script parsing

### Build System Flow

1. `NewBuildObject()` resolves build type from name/version:
   - Conda: no build script found → install via micromamba
   - Script: shell script in `build-scripts/` → runs install() function
   - Def: Apptainer definition file → builds container

2. Build scripts use headers for metadata:
   - `#DEP:name/version` - dependency declaration
   - `#SBATCH` - SLURM job parameters (also cross-translated for PBS)
   - `#ENV:VAR=$app_root` - environment variables for overlay
   - `#INTERACTIVE:prompt` - requires user input during build

3. Overlays are stored as `.sqf` (SquashFS, read-only) or `.img` (ext3, writable)

### Configuration

Configuration uses a multi-level priority system (highest to lowest):
1. Command-line flags
2. Environment variables (`CONDATAINER_*` prefix)
3. User config: `~/.config/condatainer/config.yaml`
4. Portable config: `<install-dir>/config.yaml` (auto-detected when binary lives in `<dir>/bin/`)
5. System config: `/etc/condatainer/config.yaml`
6. Hardcoded defaults

Key settings: `apptainer_bin`, `scheduler_bin`, `submit_job`, `extra_base_dirs`, `build.*` (ncpus, mem_mb, time, tmp_size_mb, compress_args, overlay_type)

### Data Directory Search

Images, build scripts, and helper scripts are searched in this order:
1. Extra base dirs: user-specified via `CONDATAINER_EXTRA_BASE_DIRS` (colon-separated) or config `extra_base_dirs`
2. Portable dir: `<install-dir>/` (auto-detected from executable location)
3. Scratch dir: `$SCRATCH/condatainer/` (HPC systems with `$SCRATCH` set)
4. User XDG dir: `~/.local/share/condatainer/`

Each base dir contains subdirectories: `images/`, `build-scripts/`, `helper-scripts/`. Write operations use the first writable directory in search order.

### Scheduler Integration

The scheduler package provides a unified interface for SLURM and PBS:
- Auto-detects scheduler from environment
- Parses script directives (#SBATCH, #PBS) into normalized ScriptSpecs
- Cross-scheduler translation: SLURM scripts can run on PBS systems

## Build Scripts Convention

Build scripts live in `build-scripts/` with naming:
- Apps: `name/version` (e.g., `cellranger/9.0.1`)
- References: `assembly/data-type/version` (e.g., `grch38/star/2.7.11b/gencode47-101`)

Shell scripts must define an `install()` function. Available variables: `$NCPUS`, `$target_dir`, `$tmp_dir`, `$app_name`, `$version`.

## Key Patterns

- Global config singleton: `config.Global`
- Error types: `ApptainerError`, `ValidationError` with structured fields
- Console output: Use `utils.Print*` functions (PrintMessage, PrintWarning, PrintError, PrintDebug)
- Path handling: Always use absolute paths; `config.Get*Dir()` functions for standard locations
- File permissions: Use `utils.PermFile` (0664), `utils.PermDir` (0775), and `utils.PermExec` (0775) constants for consistent group access on HPC systems
