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
- `apptainer/` - Apptainer binary wrapper (exec, build, version detection)
- `build/` - Build system: resolves name/version → Conda, Script, or Def build type; dependency graphs; remote script fetching
- `config/` - Multi-level config (flags > env > user > extra-root > app-root > system > defaults), data directory search
- `container/` - Container setup pipeline: overlay resolution, bind dedup, env collection, GPU detection
- `exec/` - Ephemeral container execution
- `helper/` - Helper service job lifecycle: resolve params, submit to scheduler (or run headless), monitor state, JSONL run history
- `overlay/` - Overlay image CRUD (ext3/SquashFS), resize, chown, locking
- `proxy/` - SSH tunnel + dual-protocol proxy management for HPC compute nodes
- `scheduler/` - HPC scheduler abstraction (SLURM, PBS, LSF, HTCondor); auto-detection, directive parsing, cross-scheduler translation
- `server/` - Dashboard HTTP server (web UI + REST API, SSE log streaming)
- `utils/` - Console output (`Print*`), file ops, downloads, script parsing

## Recipes

Recipes live in collections listed in the ordered `sources` config (first match wins, like `PATH`);
`catalog/` resolves a name to a recipe and walks its dependency graph. A recipe is `<name>/<version>`,
where the name may carry slashes (`grch38/genome/gencode` is name `grch38/genome`, version `gencode`).
Four kinds: `base` (produces the container root), `os`, `app` (contributes to `PATH`), `data`.

A recipe runs top to bottom as `bash -euo pipefail <recipe>` — no `install()` wrapper.
Available vars: `$CNT_NAME`, `$CNT_VERSION`, `$CNT_KIND`, `$CNT_PREFIX` (where the payload goes),
`$CNT_TMP` (also `$TMPDIR`), plus the scheduler's normalized `$NCPUS`, `$MEM`, `$MEM_GB`.

`#DEP:` is a **build** dependency only — what must be mounted while the recipe runs. It is not
recorded in the artifact and never re-expanded at run time; there is no runtime dependency tree.
So an `app` is self-contained (a conda env or a prebuilt package carrying its own libraries), and
`data` is the kind that normally has deps, since producing an index needs the producing tool.

Metadata headers: `#DEP:name/version` or `#DEP:name/version>=min` (build deps; preferred version is
implicit upper bound, so valid range is `[min, version]`), `#SBATCH`/`#PBS`/`#BSUB` (scheduler job params),
`#ENV:VAR={prefix}/sub  ## note` (env vars; `{prefix}` is filled with the mount root at load time),
`#INPUT:prompt` (user input, fed on stdin in order — read with `IFS= read -r VAR`), `#PH:`/`#TARGET:` (templates),
`#WHATIS:`, `#URL:`, `#TYPE:`.

Overlays are stored as `.sqf` (SquashFS, read-only) or `.img` (ext3, writable).

## Data Directory Search Order

1. `extra_image_dirs` / `extra_helper_dirs` (config keys)
2. `CNT_EXTRA_ROOT` (group/lab root, env only)
3. `CNT_ROOT` / `<install-dir>/` (app-root, auto-detected)
4. `$SCRATCH/condatainer/`
5. `~/.local/share/condatainer/`

Each contains `images/` and `helper-scripts/`. Writes go to first writable dir.
Recipes are not searched here — they come from the ordered `sources` list.

## Helper Scripts

Bash scripts in [`cnt-scripts/helpers/`](https://github.com/Justype/cnt-scripts) launch interactive services inside CondaTainer on HPC. Modes: `headless/` (direct) and `<scheduler>/` (submit + SSH tunnel). See `helpers/README.md` for details.

## File Locking

`exec`/`run` hold `LOCK_SH` on `.sqf`/`.sif` files during execution (`.img` skipped — Apptainer flocks those itself); `remove` and `build --update` probe `LOCK_EX` before modifying. See `internal/image/lock.go`.

## Coding Rules

- When editing a function, update its doc comment to match. Keep comments concise and behavior-first — say what it does; give a reason only when the behavior is surprising.
- When editing code, check whether the change is reflected in docs and README files, and update them when needed: `docs/manuals/condatainer.md` and relevant `docs/` pages for UX changes (flags, output format, command behaviour), and the nearest `README.md` (e.g. `internal/helper/README.md`) for package-level changes.

## Key Patterns

- Global config singleton: `config.Global`
- Error types: `ApptainerError`, `ValidationError` with structured fields
- Console output: `utils.PrintMessage`, `PrintWarning`, `PrintError`, `PrintDebug`
- Always use absolute paths; `config.Get*Dir()` for standard locations
- File/dir creation: use the wired helpers — `utils.CreateFileWritable` (files), `utils.MkdirAllShared` (dirs), `utils.MakeExecutable` (make a file executable). They apply `utils.PermFile`/`PermDir` (umask-subject) and call `ShareWithParentGroup` so children inherit group-write inside a shared `2775` install while personal installs stay umask-default. Never force perms with a bare `os.Chmod(path, utils.Perm*)` — see `internal/utils/files.go`.
