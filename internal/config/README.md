# config

Multi-level configuration management and data directory search following XDG Base Directory specification.

## Architecture

```
config.go       Global config singleton, defaults, structure
datapaths.go    Data directory search paths and resolution
persist.go      Configuration file loading (Viper)
```

## Configuration Hierarchy

Priority order (highest to lowest):
1. Command-line flags
2. Environment variables
3. User config (`~/.config/condatainer/config.yaml`)
4. Portable config (`<install>/config.yaml`)
5. System config (`/etc/condatainer/config.yaml`)
6. Hardcoded defaults

## Global Config

```go
config.Global  // Singleton instance
```

**Fields:**
- `Debug`, `SubmitJob`, `Version`
- `ProgramDir`, `LogsDir`
- `ApptainerBin`, `SchedulerBin`
- `Branch`, `PreferRemote` (remote script fetching)
- `Scheduler` - Default resource specs (`Ncpus`, `MemMB`, `Time`, `Nodes`, `Ntasks`)
- `Build` - Build settings (`DefaultCPUs`, `DefaultMemMB`, `TmpSizeMB`, `CompressArgs`, `OverlayType`)

## Data Directory Search

Follows XDG + HPC conventions with multiple search paths:

**Search Order:**
1. `CONDATAINER_EXTRA_BASE_DIRS` / `extra_base_dirs` config
2. Portable dir (auto-detected from install location)
3. Scratch dir (`$SCRATCH/condatainer/`)
4. User dir (`$XDG_DATA_HOME/condatainer/` or `~/.local/share/condatainer/`)

Each contains: `images/`, `build-scripts/`, `helper-scripts/`

**Write Operations:** First writable directory in search order

## Usage

### Initialization

```go
// Load defaults
config.LoadDefaults("/path/to/condatainer")

// Load config files (Viper)
config.InitConfig()
```

### Directory Resolution

```go
// Get base image path (search all paths)
baseImage := config.GetBaseImage()
findPath := config.FindBaseImage()

// Get writable path for new images
writePath, err := config.GetWritableImagesDir()
baseImagePath, err := config.GetBaseImageWritePath()

// Build scripts
scriptPath, err := config.FindBuildScriptFile("cellranger/9.0.1")

// Helper scripts
helperPath := config.FindHelperScript("jupyter")
```

### Data Paths

```go
config.InitDataPaths()  // Initialize GlobalDataPaths

// Search paths (ordered)
config.GlobalDataPaths.ImagesDirs
config.GlobalDataPaths.BuildScriptsDirs
config.GlobalDataPaths.HelperScriptsDirs

// Individual directories
config.GetExtraBaseDirs()
config.GetPortableDataDir()
config.GetScratchDataDir()
config.GetUserDataDir()
config.GetSystemDataDir()
config.GetUserStateDir()  // For instance state
```

## Environment Variables

- `CONDATAINER_EXTRA_BASE_DIRS` - Colon-separated extra data dirs
- `XDG_DATA_HOME` - User data directory
- `XDG_CONFIG_HOME` - User config directory
- `XDG_STATE_HOME` - User state directory
- `XDG_DATA_DIRS` - System data directories
- `SCRATCH` - HPC scratch directory

## Configuration File

Location: `~/.config/condatainer/config.yaml`

```yaml
debug: false
submit_job: true
apptainer_bin: "apptainer"
branch: "main"
prefer_remote: false

extra_base_dirs:
  - "/custom/path/condatainer"

scheduler:
  ncpus: 4
  mem_mb: 16384
  time: "4h"
  nodes: 1
  ntasks: 1

build:
  default_cpus: 4
  default_mem_mb: 16384
  default_time: "4h"
  tmp_size_mb: 3000
  compress_args: "-comp zstd -Xcompression-level 3"
  overlay_type: "squashfs"
```

## Portable Installation

When CondaTainer is installed in a portable location:
- Auto-detects `<install>/` as portable data dir
- Searches `<install>/images/`, `<install>/build-scripts/`, etc.
- Allows self-contained deployments without user home directory

## Constants

- `VERSION` - Current version
- `GitHubRepo` - GitHub repository (`Justype/condatainer`)
- `PrebuiltBaseURL` - Prebuilt image download URL
