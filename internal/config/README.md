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
- `DefaultDistro` - Base OS slug (e.g. `"ubuntu24"`)
- `Branch`, `PreferRemote` (remote script fetching)
- `ParseModuleLoad` - Treat `module load` lines as `#DEP` dependencies
- `Scheduler scheduler.ResourceSpec` - Default scheduler specs (`Nodes`, `TasksPerNode`, `CpusPerTask`, `MemPerNodeMB`, `Time`)
- `Build BuildConfig` - Build settings (`Defaults scheduler.ResourceSpec`, `TmpSizeMB`, `CompressArgs`, `OverlayType`)

## Data Directory Search

**Search Order:**
1. `CNT_EXTRA_BASE_DIRS` / `extra_base_dirs` config
2. Portable dir (auto-detected from install location)
3. Scratch dir (`$SCRATCH/condatainer/`)
4. User dir (`$XDG_DATA_HOME/condatainer/` or `~/.local/share/condatainer/`)

Each contains: `images/`, `build-scripts/`, `helper-scripts/`

**Write Operations:** First writable directory in search order

## Usage

```go
config.LoadDefaults("/path/to/condatainer")  // load defaults
config.InitViper()                           // load config files

config.GetBaseImage()                        // base image path (search all)
config.GetWritableImagesDir()                // writable images directory
config.FindImage("cellranger--9.0.1.sqf")   // search all image paths
config.FindHelperScript("jupyter")           // search helper script paths

config.GlobalDataPaths.ImagesDirs           // ordered image search paths
config.GlobalDataPaths.BuildScriptsDirs     // ordered script search paths
config.GetUserStateDir()                    // instance state directory
```

## Environment Variables

- `CNT_EXTRA_BASE_DIRS` - Colon-separated extra data dirs
- `XDG_DATA_HOME` / `XDG_CONFIG_HOME` / `XDG_STATE_HOME` - XDG dirs
- `SCRATCH` - HPC scratch directory

## Configuration File

Location: `~/.config/condatainer/config.yaml`

```yaml
debug: false
submit_job: true
apptainer_bin: "apptainer"
branch: "main"
prefer_remote: false
default_distro: "ubuntu24"
parse_module_load: false

extra_base_dirs:
  - "/custom/path/condatainer"

scheduler:
  nodes: 1
  tasks_per_node: 1
  ncpus_per_task: 2
  mem_per_node_mb: 8192
  time: "4h"

build:
  ncpus: 4
  mem_mb: 8192
  time: "2h"
  tmp_size_mb: 20480
  compress_args: ""   # empty = auto-detect based on apptainer version
  overlay_type: "ext3"
```

## Constants

- `VERSION` - Current version
- `GitHubRepo` - GitHub repository (`Justype/condatainer`)
- `PrebuiltBaseURL` - Prebuilt image download URL
