# config

Multi-level configuration management and data directory search following XDG Base Directory specification.

## Architecture

```
config.go       Global config singleton, defaults, structure
datapaths.go    Data directory search paths and resolution
persist.go      Configuration file loading (Viper)
```

## Configuration Hierarchy

Two independent mechanisms with different semantics:

**1. Environment variables (`CNT_*`)** — admin-level control (e.g. set in a module file).
Always **replaces** the config file value for that key entirely.

**2. Config files** — layered user/group/system config. All existing files are loaded:
- Scalar keys (`apptainer_bin`, `default_distro`, etc.): highest-priority file that sets the key wins.
- Array keys (`extra_*_dirs`, `extra_scripts_links`): **merged** across all layers (user ++ extra-root ++ root ++ system), deduplicated, user entries first.
- `channels`: **overwrite** — highest-priority config file that sets it wins (not merged).

Priority order (highest to lowest):
1. Command-line flags
2. Environment variables (`CNT_*`) — replaces, not merged
3. User config (`~/.config/condatainer/config.yaml`)
4. Extra-root config (`$CNT_EXTRA_ROOT/config.yaml`, group/lab layer)
5. App-root config (`$CNT_ROOT/config.yaml` or `<install>/config.yaml`)
6. System config (`/etc/condatainer/config.yaml`)
7. Hardcoded defaults

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
- `Notification` - Notification method when a helper job starts (default: `""` = none). Values: `"bell"` (terminal bell), `"email"` (scheduler email directive), ≥5-char string (ntfy.sh topic, fires from compute node), `""` or `"none"` (silent).
- `Scheduler scheduler.ResourceSpec` - Default scheduler specs (`Nodes`, `TasksPerNode`, `CpusPerTask`, `MemPerNodeMB`, `Time`)
- `Build BuildConfig` - Build settings (`Defaults scheduler.ResourceSpec`, `TmpSizeMB`, `CompressArgs`, `OverlayType`)

## Data Directory Search

**Search order for images:**
1. `extra_image_dirs` — explicit image directories (direct paths, support `:ro`/`:rw`)
2. `CNT_EXTRA_ROOT` → `<extra-root>/images/` (group/lab layer)
3. Root dir → `$CNT_ROOT/images/` or `<install>/images/`
4. Scratch dir → `$SCRATCH/condatainer/images/`
5. User dir → `$XDG_DATA_HOME/condatainer/images/` or `~/.local/share/condatainer/images/`

**Search order for scripts:**
**Build scripts:**
1. `extra_build_dirs` — explicit build-scripts directories (direct paths)
2. `CNT_EXTRA_ROOT`, Root dir, Scratch dir, User dir (same pattern)

**Helper scripts:**
1. `extra_helper_dirs` — explicit helper-scripts directories (direct paths)
2. `CNT_EXTRA_ROOT`, Root dir, Scratch dir, User dir (same pattern)

**Write operations:**
- **Images / helpers**: first writable directory in search order. `:ro` entries and `extra_build_dirs` are always skipped. Personal dirs (scratch, user) are always created on first use. Shared dirs (extra-root, root): subdirs (`images/`, `build-scripts/`, etc.) are auto-created if the parent directory already exists — the parent itself is never auto-created.
- **Cache**: always written to a personal directory (scratch → user cache) to avoid cross-user pollution. Shared dirs are never written to.

**`:ro` / `:rw` markers** (image and helper dirs, config file only):
- `:ro` — search-only; condatainer never writes here even if filesystem allows it
- `:rw` — explicit writable annotation (same as no marker; for documentation clarity)
- Only applies to `extra_image_dirs` and `extra_helper_dirs`; `extra_build_dirs` entries are plain paths

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
config.GetWritableCacheDir()                // personal writable cache directory (scratch → user)
config.GetCacheSearchPaths()                // personal cache search paths
```

## Environment Variables

All multi-value env vars use `|` as separator. `CNT_EXTRA_BUILD_DIRS` also accepts `:`.

| Variable | Separator | Description |
|---|---|---|
| `CNT_ROOT` | — | Cluster/system root dir (loads `config.yaml` + data dirs; replaces bin/ heuristic) |
| `CNT_EXTRA_ROOT` | — | Group/lab root dir (single path; loads `config.yaml` + data dirs) |
| `CNT_EXTRA_IMAGE_DIRS` | `\|` | Extra image directories; entries support `:ro`/`:rw` |
| `CNT_EXTRA_BUILD_DIRS` | `\|` or `:` | Extra build-scripts directories |
| `CNT_EXTRA_HELPER_DIRS` | `\|` | Extra helper-scripts directories; entries support `:ro`/`:rw` |
| `CNT_EXTRA_SCRIPTS_LINKS` | `\|` | Extra remote build script source URLs |
| `CNT_CHANNELS` | `\|` or `:`  | Conda channels |
| `CNT_NOTIFICATION` | — | Override `notification` for the current session (e.g. `bell`, `email`, ntfy.sh topic) |
| `CNT_TMPDIR` | — | Override build temp directory |
| `SCRATCH` | — | HPC scratch directory (`$SCRATCH/condatainer/`) |
| `XDG_DATA_HOME` / `XDG_CONFIG_HOME` / `XDG_STATE_HOME` | — | XDG base dirs |

## Configuration File

Location: `~/.config/condatainer/config.yaml`

```yaml
apptainer_bin: "apptainer"
scheduler_bin: ""         # auto-detect if empty

prefer_remote: false
default_distro: "ubuntu24"
parse_module_load: false

# Extra directories — team/lab use via module file or config
extra_image_dirs:
  - "/shared/lab/images:ro"   # search-only shared store
  - "/fast/scratch/images"    # writable personal store
extra_build_dirs:
  - "/shared/lab/scripts"
extra_helper_dirs:
  - "/shared/lab/helpers"
# For a group/lab root with standard layout (images/, build-scripts/, etc.),
# set CNT_EXTRA_ROOT=/proj/condatainer in the module file instead.

channels:
  - conda-forge
  - bioconda

build:
  ncpus: 4
  mem: 8192   # MB
  time: "2h"
  compress_args: ""   # auto-detect: zstd-medium (apptainer>=1.4), lz4 (older), gzip (singularity)
```

## Constants

- `VERSION` - Current version
- `GitHubRepo` - GitHub repository (`Justype/condatainer`)
- `PrebuiltBaseURL` - Prebuilt image download URL
