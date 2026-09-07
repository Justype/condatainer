# config

Multi-level configuration management and data directory search following XDG Base Directory specification.

## Architecture

```
config.go       Global config singleton, defaults, structure
datapaths.go    Data directory search paths and resolution
persist.go      Configuration file loading (Viper)
sources.go      Recipe collections, the catalog handle, and the default distro
```

## Recipe sources and the default distro

`sources` is an ordered list of recipe collections, read from every config layer
and concatenated strongest first, so a user entry shadows a site entry of the
same name. Each entry is a single-key mapping — `- lab: /shared/lab/recipes` —
because both the order and the handle are load-bearing and a plain map gives
neither. `CNT_SOURCES` overrides the lot: `cnt=https://…|lab=/shared/lab`.

The public `cnt` collection ships as a default *value*, not as a fallback the
resolver reaches for, and is **appended** rather than prepended, so every
configured entry outranks it. A site redefining `cnt` replaces it outright,
which is how the handle is pointed elsewhere without rewriting the `#DEP:` lines
that name it.

**Unreachable sources are reported, not fatal.** The remaining sources still
answer, and a compute node with no route out is ordinary — but silence is not an
option either, because sources are first-wins: an unreachable one promotes the
next source's recipe or falls through to conda, and the build would otherwise
look normal. `WarnUnreachableSources` warns once per process, and must be called
*after* the catalog has been consulted, since a source's `Err` is set when it is
first read rather than when it is opened.

`avail` and `create` can restrict that configured list with repeatable source
handles. Flag order becomes lookup precedence; dependencies and the default
distro use the same restricted catalog:

```text
condatainer avail -s lab
condatainer create -s lab -s cnt star/2.7.11b
```

With no `--source`, the full configured list is used. An unknown handle is an
error. A missing recipe still follows the normal Conda fallback.

**The default distro is recorded once and never revised.** `EnsureDefaultDistro`
writes it from the first source declaring a `default_distro` the first time one
is needed. Changing it rebuilds the container root and every `os` overlay
stacked on it, so following an upstream bump would invalidate a whole set of
images on an ordinary update; a later default is something the user opts into
with `config set default_distro`.

`ResolvedDefaultDistro` reads config alone and never opens the catalog. It is
the bare-name prefix for installed overlays (`build-essential` →
`ubuntu24/build-essential`), so it is called on offline paths like `list` and
`info`, where reaching for a source descriptor would mean a network fetch to
expand a local name. The `default_distro` fallback belongs where a catalog is
already open — base image resolution, in `internal/build`.

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
- `Notification` - Notification method when a helper job starts (default: `""` = none). Values: `"bell"` (terminal bell), `"email"` (scheduler email directive), ≥5-char string (ntfy.sh topic, fires from compute node), `""` or `"none"` (silent).
- `ProxyPerJob` - Auto-start a per-job SOCKS5 proxy inside submitted jobs when no active proxy is found (`proxy_perjob` config key, default: `false`)
- `Scheduler scheduler.ResourceSpec` - Default scheduler specs (`Nodes`, `TasksPerNode`, `CpusPerTask`, `MemPerNodeMB`, `Time`)
- `Build BuildConfig` - Build settings (`Defaults scheduler.ResourceSpec`, `AppTmpOverlay`, `AppTmpOverlaySizeMB`, `CompressArgs`, `BlockSize`, `DataBlockSize`)

## Data Directory Search

Reads go nearest-first and writes furthest-first, so a build lands as far out as
permissions allow (one copy for the whole group) while a personal build shadows a
shared one for the person who made it. Order within a tier is the same both ways.

**Search order for images** (personal, then shared):
1. Scratch dir → `$SCRATCH/condatainer/images/`
2. User dir → `$XDG_DATA_HOME/condatainer/images/` or `~/.local/share/condatainer/images/`
3. `CNT_EXTRA_ROOT` → `<extra-root>/images/` (group/lab layer)
4. Root dir → `$CNT_ROOT/images/` or `<install>/images/`

**Search order for scripts:**

Recipes are not searched for on disk — they come from the catalog's `sources`
(see [sources.go](sources.go)), which is an ordered list of collections, not a
data directory.

The default collection is always reachable: unless something already answers to
the handle `cnt`, `https://raw.githubusercontent.com/condatainer/recipes/main` is
**appended** to whatever the config lists. Appended, never prepended, so every
configured entry outranks it — and defining your own `cnt` replaces it, which is
how a site points the handle elsewhere without rewriting the `#DEP:` lines that
name it.

**Helper scripts:**
1. `CNT_EXTRA_ROOT`, Root dir, Scratch dir, User dir (same pattern, with `helper-scripts/`)

**Write operations:**
- **Images / helpers**: first writable directory in the *reverse* order — extra-root → root → scratch → user. Personal dirs (scratch, user) are always created on first use. Shared dirs (extra-root, root): subdirs (`images/`, `helper-scripts/`) are auto-created if the parent directory already exists — the parent itself is never auto-created.
- **Cache**: always written to a personal directory (scratch → user cache) to avoid cross-user pollution. Shared dirs are never written to.

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

All multi-value env vars use `|` as separator.

| Variable | Separator | Description |
|---|---|---|
| `CNT_ROOT` | — | Cluster/system root dir (loads `config.yaml` + data dirs; replaces bin/ heuristic) |
| `CNT_EXTRA_ROOT` | — | Group/lab root dir (single path; loads `config.yaml` + data dirs) |
| `CNT_SOURCES` | `\|` | Recipe collections as `name=base` pairs; overrides the `sources` config key |
| `CNT_CHANNELS` | `\|` or `:`  | Conda channels |
| `CNT_NOTIFICATION` | — | Override `notification` for the current session (e.g. `bell`, `email`, ntfy.sh topic) |
| `CNT_PROXY_PERJOB` | — | Override `proxy_perjob` for the current invocation (`1` = enable) |
| `CNT_TMPDIR` | — | Override build temp directory |
| `SCRATCH` | — | HPC scratch directory (`$SCRATCH/condatainer/`) |
| `XDG_DATA_HOME` / `XDG_CONFIG_HOME` / `XDG_STATE_HOME` | — | XDG base dirs |

## Configuration File

Location: `~/.config/condatainer/config.yaml`

```yaml
apptainer_bin: "apptainer"
scheduler_bin: ""         # auto-detect if empty

default_distro: "ubuntu24"

# Recipe collections, in order — first match wins, like PATH.
# `cnt` is appended automatically; list it yourself only to point it elsewhere.
sources:
  - lab: "/shared/lab/recipes"

# For a group/lab root with standard layout (images/, helper-scripts/),
# set CNT_EXTRA_ROOT=/proj/condatainer in the module file.

channels:
  - conda-forge
  - bioconda

build:
  ncpus: 4
  mem: 8192   # MB
  time: "2h"
  compress_args: "-comp zstd -Xcompression-level 8"   # zstd-medium, always — every reader is a version-checked apptainer
```

## Constants

- `VERSION` - Current version
- `GitHubRepo` - GitHub repository (`Justype/condatainer`)
- `PrebuiltBaseURL` - Prebuilt image download URL
