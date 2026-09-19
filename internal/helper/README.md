# internal/helper

Orchestrates the full lifecycle of a helper service job: resolve params, check/create overlays, submit to scheduler (or run headless), monitor NFS state files, and clean up.

Scheduler-vs-headless is normally automatic (`config.Global.SubmitJob`, scheduler binary detection, `sched.IsInsideJob()`), but `RunOptions.NoSubmit` (CLI: `condatainer helper --no-submit`; dashboard: the "Run headless" checkbox, shown only when a scheduler is actually available) forces headless for one run without touching global config — the server needs this since it is one long-lived process serving concurrent requests. A failed `Submit` never falls back to headless automatically; the error names the override instead, since silently running a compute-node job on the login node without consent is bad etiquette on a shared system.

Every `HelperRun` records `Runner` — `"local"` for headless, or the scheduler type (`"slurm"`, `"pbs"`, `"lsf"`, `"htcondor"`) for a submitted job — set once at creation, the same convention `internal/image/producer.Info.Runner` uses for build locks. It exists because `JobID` emptiness and the initial `"pending"`/`"starting"` status both leak the same fact only until a run reaches `"running"`, after which nothing else says which path it took.

`RunOptions.Account`/`.Partition` resolve the same way on every entry point — CLI flag (`-A`/`--account`, `-p`/`--partition` on `cmd/helper.go`; the terminal `PromptSettings` table also lets either be edited interactively before launch) or the dashboard's `cfg-account`/`cfg-partition` fields, falling back to `config.Global.Scheduler.Account`/`.Partition` when left unset — applied once, in `buildHelperScriptSpecs` (`run.go`), never earlier. There is no script-header tier (no `#ACCOUNT:`/`#PARTITION:`): unlike `#NCPUS:`/`#MEM:`/`#TIME:`/`#GPU:`, which describe what the workload needs, account/partition describe the user's own cluster access, which a script has no way to know.

## Files

| File | Purpose |
|------|---------|
| `run.go` | `RunHelper()` — top-level entry point; also server auto-start, scheduler detection, wrapper generation |
| `state.go` | NFS state files (`ready`, `done`, `messages`), JSONL history, `HelperRun` struct |
| `monitor.go` | `PollOnce()` — shared non-blocking read of all three NFS state files |
| `params.go` | `ApplyParamFlags()`, `PromptMissingParams()` — resolve `#PARAM:` headers |
| `config.go` | Per-helper saved config (`~/.local/share/condatainer/helper-config/<name>.conf`) |
| `shared_history.go` | Thin `helper.RecordUsed`/`ListUsed`/`ListAll` re-exports of `internal/helperhistory` |

---

## Shared, cross-user helper-setup history (project-scoped)

Personal reuse (`state.go`'s JSONL history) is per-`$HOME`: two people running
`rstudio-server` from the same project directory each rediscover which
overlays make it work, with nothing recording that the other already solved
it. `internal/helperhistory` (a separate package — see below for why) fixes
that by recording one entry per distinct **overlay combination** per
`(helper name, project-root-relative location)`, shared through the
project's own `cnt-lock/.helper-history/`.

**Gated on an existing project — never created as a side effect.** This only
ever reads or writes when the launch's `cwd` already stands inside a project
(`project.StandingAt` is non-nil); it never runs `project init`, and outside
a project it silently falls back to exactly today's personal-JSONL-only
behavior. A project's mere presence already changes helper resolution
(`CheckRequiredOverlays` switches from ambient install to strict,
no-fallback `Standing.ResolveNames`), so auto-creating one to make a place
for this history could break the very launch it was meant to help.

**Storage is one flat directory, one JSON file per combination — no
subdirectories, no hash.** A combination's identity (`Helper`, `Location`,
`Overlays`) lives entirely inside the file; the filename
(`<helper>-<location-with-slashes-converted>-<unix-nano>.json`) exists only
to be unique on disk, never parsed back. Every read — the fresh-start lookup
for one location, or the project-wide aggregate below — reads every file in
the directory and filters in memory, never a targeted lookup; the file count
is bounded by how many genuinely distinct setups a project's helpers ever
use, not by run count.

**Dedup and recency come from content and mtime, never a hash or a written
field.** A write reads every existing file, and on a
`Helper`/`Location`/sorted-`Overlays` match bumps that file's mtime instead
of writing a new one. The bump goes through `unix.UtimesNanoAt` with
`UTIME_OMIT` (atime) and the `UTIME_NOW` sentinel (mtime) — **never
`os.Chtimes`**, which always passes an explicit timestamp that `utimensat(2)`
only lets the file's owner set. `UTIME_NOW` waives that and needs only write
access, which `ShareWithParentGroup` (already called inside
`utils.CreateFileWritable`) guarantees every member of a shared, group-writable
`cnt-lock/` has — the whole point, since this is written by whoever's UID
happens to run the helper next.

**What counts as the combination:** two lists kept apart: `required`, the
helper's `#REQUIRED_OVERLAYS:` names as declared with the launch's `{KEY}`
params filled in — a `{R_VERSION}` choice is recorded as `r/4.4.3`, the same
string that keys its pin — and `overlays`, the `-o` list the user added, in the
form `normalizeOverlayForHistory` produces (a catalog `name/version`, or a
path relative to the launch `cwd`). Both take part in dedup. Every overlay in either list counts as used, so both
feed the project's status; only `overlays` is offered back by the reuse lookup,
since the required ones follow from the helper and its params. Resolved
`#PARAM:` values and the autoloaded env overlay (`-e`) are not recorded — a
`#PARAM:` choice is a per-launch settings question already served by saved
config, and every helper autoloads its env overlay unconditionally, so it is
not a decision two setups could differ on.

**No browsing view, anywhere.** There is no `helper history <name>` command
and no dashboard tab for past combinations. The data has exactly two
consumers: the fresh-start reuse lookup (`ListUsed`, one helper, one location,
offered as a confirmed default — never applied silently) and
`internal/project`'s project-wide "used but not pinned" / "which helper uses
this manual pin" signals (`ListAll`, aggregated, never displayed as raw
combinations). See `internal/project/README.md`, "Used but not pinned, and who
uses a manual pin".

**Why a separate package (`internal/helperhistory`) instead of living here.**
`internal/project` needs to read this data too (for the aggregate signals
above), but `internal/helper` already imports `internal/project`
(`CheckRequiredOverlays` calls `project.StandingAt`) — so `internal/project`
importing `internal/helper` back would cycle. `internal/helperhistory` is a
leaf package with no dependency on either, so both can depend on it; this
package's `RecordUsed`/`ListUsed`/`ListAll` are thin re-exports so a CLI or
dashboard caller never has to import `internal/helperhistory` directly.

---

## Helper Script Headers

Helper scripts are plain bash files with metadata in comment headers. All headers are parsed before the job is submitted.

### Resource headers

Parsed by `scheduler.ReadScriptSpecs()`. Standard scheduler directives also work (`#SBATCH`, `#PBS`, `#BSUB`).

```bash
#NCPUS: 4          # default CPU count
#MEM:   16G        # default memory       — bare number defaults to GB (16 = 16G)
#TIME:  12:00:00   # default walltime     — bare number defaults to hours (12 = 12h)
                   # also accepts: 12h, 2h30m, 4d, 1-12:00:00
#GPU:   a100:1     # optional GPU spec (type:count)
```

`#GPU:` drives two independent things: it becomes a scheduler resource request
(`--gpus-per-node=...`) so the job lands on a GPU node, and `buildCondatainerCmd`
appends `--gpu` to the generated `condatainer exec` line so that node's
container still gets `--nv`/`--rocm` even if it has `autoload_gpu` disabled.

### `#DESC:` — short description

Shown in `condatainer helper --list` and the server dashboard.

```bash
#DESC: Jupyter Lab
```

### `#PARAM:` — helper-specific parameters

Defines a configurable variable. The Go runner resolves each param from (in order): CLI flag → saved config → script default.

```
#PARAM: KEY=default --long-flag,-s "Description shown in --help and prompts"
```

Three forms for the default value:

| Default | Behavior |
|---------|----------|
| `KEY=` | **Required** — user must supply a value via CLI flag, saved config, or interactive prompt. |
| `KEY=?` | **Auto** — filled from the first entry of the matching `#VALUE:` list; empty (optional) if no list. |
| `KEY=literal` | **Literal** — always uses `literal` unless overridden by flag or saved config. |

- `--long-flag,-s` — optional CLI flags. All combinations are valid: both (`--flag,-s`), long-only (`--flag`), short-only (`-s`), or omitted entirely (prompt-only param).
- `"Description"` — shown in `condatainer helper <name> --help` and the settings table.

The resolved value is exported as `$KEY` inside the job wrapper.

**Auto default (`KEY=?`) example** — picks the latest R version from `#VALUE:` automatically; user can override with `-r`:

```bash
#PARAM: POSIT_R=? --rversion,-r "R version"
#VALUE: POSIT_R=4.5.3,4.5.2,4.5.1,...
#REQUIRED_OVERLAYS: r{POSIT_R} rstudio-server
```

**First-run save pattern** — for params computed at runtime (e.g. tunnel names), leave the default empty (`KEY=`) and compute+save inside the script on first run:

```bash
if [ -z "${MACHINE_NAME:-}" ]; then
    MACHINE_NAME="$(whoami)-$(hostname -s)"
    condatainer helper "$CNT_HELPER_NAME" config set MACHINE_NAME "$MACHINE_NAME"
fi
```

Subsequent runs: Go reads the saved config and injects `$MACHINE_NAME` before the script runs.

### `#VALUE:` — allowed values / version list

Provides the ordered list of valid values for a `#PARAM:` key, used for the interactive prompt selector and the web UI dropdown. Updated by `auto.py`.

```bash
#VALUE: CONDA_PYTHON=3.13.0,3.12.7,3.11.9,...
```

Two formats:
- **Comma-separated** `KEY=val1,val2,...` — version list, **sorted descending** automatically (newest first). Order in the script doesn't matter.
- **Pipe-separated** `KEY=opt1 | opt2 | opt3` — ordered option labels, **not sorted**. Order in the script is preserved.

In both, `a-b` expands to every integer in the range and `*` marks the list open-ended (sorted last, so index 0 stays the default). Parsed by `catalog.ParseValues`, which recipe `#PH:` lists also use — the two headers are one dialect.

### `#IMG_PACKAGES:` — conda packages for guided overlay creation

Presence signals that a writable conda-env overlay is required. Value is the default package list for guided creation; `{KEY}` tokens are substituted from resolved `#PARAM:` values.

```bash
#IMG_PACKAGES: python={CONDA_PYTHON} jupyterlab   # overlay required + package check
#IMG_PACKAGES:                                     # overlay required, no package check
```

`#IMG_PACKAGES:` is really two requirements: the packages must be installed, and there
must be a *writable* form to install into if they aren't. `ResolveEnvOverlayInDir`
resolves `.img` and `.sqf` as two forms of one environment overlay, `.img` first, falling
back to the read-only `.sqf` when that's all that exists — so the first requirement can
be satisfied by a bare snapshot, but the second cannot: `PlanRun` still refuses a resolved
`.sqf` (`!utils.IsImg(opts.EnvImg)`) and names it in the message, and the CLI still runs
the guided-creation wizard in that case, because there is nowhere to write to yet.

- **With a writable `.img` resolved**: if packages are listed, Go runs a pre-submission
  install check (fatal if packages not found). If the value is empty, the overlay is
  still required but no check is performed.
- **Without one** (nothing found, or only a read-only `.sqf`): Go prompts to create a
  writable overlay (guided overlay creation flow). This is not wasteful even when a
  snapshot already exists: `exec.CreateCondaOverlay` looks up the paired snapshot against
  the *final* destination path and mounts it during the scratch install, so packages the
  snapshot already provides are not reinstalled — only the incremental diff lands in the
  new `.img`'s own layer.

The install check itself (`checkPackages` → `container.PairedPackages`) reads the
overlay's paired snapshot as well as the `.img` itself, unioning both package sets before
checking — or, when `EnvImg` resolved to a bare `.sqf`, checks that directly. A thin
`.img` sitting on a frozen snapshot holds only the newest delta in its own `conda-meta` —
checking the `.img` alone would report everything the snapshot carries as "not installed"
and fatally block every `#IMG_PACKAGES:` helper against a perfectly good, just-frozen
environment. `CheckEnv`'s `SizeMB` and `Snapshot` fields (consumed by the dashboard) go
through `container.PairedSize`, the sibling function for size instead of packages, for
the same reason: a thin `.img`'s own size says nothing about the environment it actually
provides.

### `#POST_INSTALL_CMD:` — post-install hook

A single shell command run inside the container after guided overlay creation installs the packages. Used to pin package versions or run setup steps.

```bash
#POST_INSTALL_CMD: mm pin add r-base
```

### `#CHECK_PATH:` — pre-submission binary check

One path per line. Checked inside the overlay (fatal) when `#IMG_PACKAGES:` is set, or in the base image (non-fatal warning) otherwise. Optional quoted message shown when the path is missing.

```bash
#CHECK_PATH: /cnt_env/bin/jupyter-lab "Jupyter Lab not found — install with: mm install jupyterlab"
#CHECK_PATH: /cnt_env/bin/R "Conda R not found — install with: mm install r-base=<version>"
```

Multiple lines accumulate independently (one path checked per line).

### `#REQUIRED_OVERLAYS:` — named SquashFS overlays

Space-separated overlay names. `{KEY}` tokens are substituted from resolved params. **Declaration order is preserved** in the Apptainer overlay stack — the last name is the topmost layer (wins on file conflicts such as `/var/lib/dpkg/status`).

Standing in a project (the launch's working directory holds `cnt-lock/`), each name resolves through that project's lock instead of by installed name — the same substitution `exec -o` makes — and a name the project has not pinned refuses rather than falling back to whatever currently answers to it. Outside a project, each name — bare, carrying a partial version, or missing its distro prefix — is resolved against what's already installed first, no network call, then built via `condatainer create <name>` (passed exactly as declared) only once nothing local satisfies it. A name that resolved to something other than itself is printed once, e.g. `rstudio-server -> ubuntu24/rstudio-server/2026.09.0-174`.

```bash
#REQUIRED_OVERLAYS: r{POSIT_R} rstudio-server build-essential
# resolves to: r4.4.3 rstudio-server build-essential
# → build-essential is topmost; its dpkg database is visible to tools like pak
```

### `#BIND:` — extra bind mounts

One bind spec per line, in `src:dest` or `src:dest:opts` form (same as Apptainer `--bind`).
`{KEY}` tokens are substituted from resolved `#PARAM:` values. `$VAR` references are left
verbatim and expanded by the shell when the wrapper runs on the compute node — any of the
[injected env vars](#environment-variables-injected-by-the-wrapper) can be used directly.

```bash
#BIND: $CNT_HELPER_STATE_DIR/vnc:$HOME/.vnc
#BIND: /shared/data:/data:ro
#BIND: /path/to/{KEY_DIR}:/mnt/custom
```

### `#SINGLETON:` — at most one running instance

When `true`, the "start new instance" option is suppressed if an instance is already running. `--new`/`--force` prints an error instead.

```bash
#SINGLETON: true
```

---

## Environment Variables Injected by the Wrapper

The Go runner injects these into every job before the helper script body runs:

| Variable | Value |
|----------|-------|
| `CNT_HELPER_NAME` | Helper script name (e.g. `jupyterlab`) |
| `CNT_HELPER_ID` | Unique run ID: `{name}-{job_id}` or `{name}-{pid}` |
| `CNT_HELPER_PORT` | Free TCP port on the compute node (helper must bind here) |
| `CNT_HELPER_BIND_ADDR` | `127.0.0.1` normally; `0.0.0.0` when `helper_bind_all` is set |
| `CNT_HELPER_BIND_ALL` | Set to `1` when `helper_bind_all` is enabled; absent otherwise |
| `CNT_HELPER_CWD` | Working directory |
| `CNT_HELPER_STATE_DIR` | NFS state directory for this run — write runtime files here |
| `CNT_HELPER_WALLTIME_SECS` | Walltime in seconds |
| `CNT_HELPER_JOB_ID` | Scheduler job ID (empty for headless) |
| `CNT_HELPER_SCRIPT_DIR` | Directory containing the helper script |
| `CNT_JOB_TMPDIR` | Node-local scratch dir (Unix sockets, etc.); cleaned up on exit |
| `$KEY` | One var per resolved `#PARAM:` key |

When writing helper scripts, use `${CNT_HELPER_BIND_ADDR:-127.0.0.1}` as the bind address so the script works correctly in both normal and `helper_bind_all` modes.

---

## Stopping a helper

A stop is a SIGTERM to the wrapper: the scheduler's cancel signals the batch shell only, and headless
runs signal the process group. The wrapper runs the container command in the background and `wait`s on
it, because bash defers a trap while a foreground command runs. Its trap forwards TERM to
`condatainer exec` and the wrapper keeps waiting, so the container unmounts its images before the
wrapper records `done` (exit code 130) and exits.

`condatainer exec` is launched with `--stop-grace 30` (`stopGraceSecs`): after TERM, Apptainer gets that
long to exit before it is killed. The value stays inside the scheduler's own TERM-to-KILL delay, which
is what ends the job when the container ignores TERM. Killing a writable image's `fuse2fs` mount leaves
the image needing recovery, so nothing in the wrapper kills the container before that delay expires.

Until the container has started there is nothing to forward to, so an earlier trap records `done` and
cleans up directly.

## NFS State Files

Written by the script on the compute node; read by the CLI monitor and server watcher on the login node.

```
~/.local/state/condatainer/helper/{id}/
    ready       — JSON written by _server_ready when service is accessible
    done        — JSON written by the wrapper on exit (always fires, even on crash)
    messages    — append-only JSONL written by _server_message
    job.log     — scheduler stdout/stderr (streamed by server dashboard)
```

### `ready` JSON

```json
{"port": 8888, "node": "cn01", "label": "Jupyter Lab", "timestamp": "...",
 "walltime_secs": 43200, "job_id": "12345",
 "url_path": "?token=abc123", "external_url": ""}
```

`port: 0` + non-empty `external_url` = no proxy tunnel needed (e.g. vscode-tunnel).

### `done` JSON

```json
{"exit_code": 0, "ts": "2024-01-01T12:00:00Z"}
```

### `messages` JSONL

One JSON object per line:

```json
{"level": "info", "text": "Open project: rstudioapi::openProject(...)", "ts": "..."}
{"level": "warn", "text": "Low memory warning", "ts": "..."}
{"level": "error", "text": "Port already in use", "ts": "..."}
```

---

## Hidden Commands Used by Scripts

Scripts call these condatainer subcommands from inside the container:

### `condatainer _server_ready`

Writes the `ready` file and appends to JSONL history.

```bash
condatainer _server_ready \
    --port "$CNT_HELPER_PORT" \
    --label "Jupyter Lab" \
    --url-path "?token=$TOKEN"

# For external-URL services (no port forwarding):
condatainer _server_ready \
    --port 0 \
    --label "VS Code Tunnel: $MACHINE_NAME" \
    --external-url "https://vscode.dev/tunnel/$MACHINE_NAME$CNT_HELPER_CWD"
```

### `condatainer _server_message`

Appends one line to the `messages` file. Shown in the terminal and the server dashboard.

```bash
condatainer _server_message "Open project: rstudioapi::openProject(\"$CNT_HELPER_CWD\")"
condatainer _server_message --level warn "Low disk space"
condatainer _server_message --level error "Failed to start"
```

### `condatainer helper <name> config`

Reads and writes the per-helper saved config. Callable from inside the container (NFS-accessible).

```bash
condatainer helper "$CNT_HELPER_NAME" config set MACHINE_NAME "$MACHINE_NAME"
condatainer helper "$CNT_HELPER_NAME" config get MACHINE_NAME
condatainer helper "$CNT_HELPER_NAME" config show   # all keys
condatainer helper "$CNT_HELPER_NAME" config path   # config file path
```

---

## Minimal Example Script

```bash
#!/bin/bash
#DESC: My Custom Service
#NCPUS: 2
#MEM: 8G
#TIME: 4:00:00
#PARAM: MODE=light --mode,-m "Interface mode: light or full"

my-service \
    --mode "$MODE" \
    --bind 127.0.0.1:"$CNT_HELPER_PORT" &
PID=$!

condatainer _server_ready \
    --port "$CNT_HELPER_PORT" \
    --label "My Custom Service ($MODE)"

wait $PID
```

## Example: conda-env service (with overlay)

```bash
#!/bin/bash
#DESC: My App (conda)
#IMG_PACKAGES: my-app={VERSION}
#CHECK_PATH: /cnt_env/bin/my-app "my-app not found — install with: mm install my-app"
#NCPUS: 4
#MEM: 16G
#TIME: 12:00:00
#VALUE: VERSION=2.0,1.9,1.8

my-app --port "$CNT_HELPER_PORT" &
PID=$!

condatainer _server_ready --port "$CNT_HELPER_PORT" --label "My App $VERSION"
wait $PID
```

## Example: external-URL service (no port forwarding)

```bash
#!/bin/bash
#DESC: My Tunnel Service
#SINGLETON: true
#NCPUS: 2
#MEM: 4G
#TIME: 24:00:00
#PARAM: NAME= --name,-n "Tunnel name (auto-set on first run)"

if [ -z "${NAME:-}" ]; then
    NAME="$(whoami)-$(hostname -s | tr '.' '-')"
    condatainer helper "$CNT_HELPER_NAME" config set NAME "$NAME"
fi

my-tunnel --name "$NAME" 2>&1 | while IFS= read -r line; do
    echo "$line"
    case "$line" in
        *"device code"*) condatainer _server_message "$line" ;;
        *"tunnel ready"*)
            condatainer _server_ready \
                --port 0 \
                --label "My Tunnel: $NAME" \
                --external-url "https://my-service.example.com/tunnel/$NAME" ;;
    esac
done
```
