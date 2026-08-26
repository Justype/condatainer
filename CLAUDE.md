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
- `registry/` - Authenticated OCI transport: upstream digest resolution plus artifact publish,
  platform-aware resolve, verification, pull, tag listing, and credential-store login/logout
- `scheduler/` - HPC scheduler abstraction (SLURM, PBS, LSF, HTCondor); auto-detection, directive parsing, cross-scheduler translation
- `server/` - Dashboard HTTP server (web UI + REST API, SSE log streaming)
- `utils/` - Console output (`Print*`), file ops, downloads, script parsing

## Recipes

Recipes live in collections listed in the ordered `sources` config (first match wins, like `PATH`);
`catalog/` resolves a name to a recipe and walks its dependency graph. A recipe is `<name>/<version>`,
where the name may carry slashes (`grch38/genome/gencode` is name `grch38/genome`, version `gencode`).
Four types: `base` (produces the container root), `os`, `app` (contributes to `PATH`), `data`.

A recipe runs top to bottom as `bash -euo pipefail <recipe>` — no `install()` wrapper.
Available vars: `$CNT_NAME` (the complete name, e.g. `samtools/1.23.1`), `$CNT_TYPE`,
`$CNT_PREFIX` (where the payload goes), `$CNT_TMP` (also `$TMPDIR`, both `/cnt_tmp`), plus the
scheduler's normalized `$NCPUS`, `$MEM`, `$MEM_GB`. There is no `$CNT_VERSION`: a recipe that
varies by version uses a `#PH:` placeholder, and one pinned to a version writes it literally.

`#DEP:` is a **build** dependency only — what must be mounted while the recipe runs. It is not
recorded in the image and never re-expanded at run time; there is no runtime dependency tree.
**Only a `data` recipe may declare one**, and validation rejects it on any other type: an `app` is
self-contained (a conda env or a prebuilt package carrying its own libraries), an `os` is
self-contained by definition, and a `base` *is* the build environment. Producing an index needs the
producing tool, which is why data is the type with deps. A `#DEP:` may name an app, data, or os, but
never a base — resolution refuses that edge.

**A *build's* `#DEP:` is a `name/version`, never an overlay path** — but a *running* script's may be
either. The asymmetry is the point: a running script mounts what it names and records nothing, while
a build's declaration becomes an edge in an artifact that has to mean the same thing on another
machine, and a path is neither resolvable there nor a key anything can regenerate. So
`#DEP: ./overlays/x.sqf` and `#DEP: env.img  ## unpinned` stay valid in an analysis script and are
rejected in a recipe.

Both build rules live in `catalog.ValidateDeps`, whose only callers are `Recipe.Validate` and
`build.FromExternalSource` — never `run` or `check`. An **external build** (`create -p -f <script>`)
answers to them exactly as a catalog recipe does, which is why they live there rather than on
`Recipe`: an external build has no recipe path to be parsed from. `catalog.IsPathDep` is the single
answer to "is this dep a path", since `catalog` owns the `Normalize`/`ParseDep` grammar that applies
only to names; its extension set must stay `utils.IsOverlay`'s, `.ext3` included.

A **version constraint is a recipe-only feature.** `#DEP:name/version>=min` lets a build reuse a
satisfying version that is already installed instead of producing another, and for reference data
the tool version often does not matter — `samtools faidx` writes the same `.fai` whichever recent
samtools ran. An **analysis script must name an exact `name/version`**: "any version in this range"
is not a claim to attach to a result, and a project lock exists to say which one. A constrained
declaration in a project script is a scan finding, not a request.

Metadata headers: `#DEP:name/version`, or `#DEP:name/version>=min` in a recipe (preferred version is
implicit upper bound, so valid range is `[min, version]`), `#SBATCH`/`#PBS`/`#BSUB` (scheduler job params),
`#ENV:VAR={prefix}/sub  ## note` (env vars; `{prefix}` is filled with the install prefix at load time),
`#INPUT:prompt` (user input, fed on stdin in order — read with `IFS= read -r VAR`), `#PH:`/`#TARGET:` (templates),
`#ARCH:noarch` (app and data script recipes only; default `native`), `#DESC:`, `#URL:`, `#TYPE:`,
`#LICENSE:` (an SPDX expression, verbatim), `#REDISTRIBUTE:` (`yes` or `no` — see *Publishing*).

**`#TARGET:` without `#PH:` is not a template — it names the artifact**, and that is how an *external*
build (`create -p <path> -f <script>`) gets a name at all. Without it the name is the `-p` basename,
which is a single component, so `key.Role`'s component match never fires and every `#DEP:` is silently
downgraded to build history. An external script declaring `#DEP:` must therefore declare `#TARGET:`,
and `FromExternalSource` refuses it otherwise. The name and the file path are separate namespaces:
`#TARGET:` fixes `/cnt/<name>` and the role classification, `-p` fixes where the `.sqf` lands. That is
why a path-addressed artifact's filename carries no naming claim — `project.LookupAt` matches it by
manifest name alone, while a flat or store scan still requires the filename to encode the name,
because there the filename *is* the address.

An annotation is read wherever it is written; position carries no meaning, so there is no header
block and no boundary. `catalog.ScanAnnotations` is the one tokenizer — recipes, user scripts and
project scanning all go through it, and a key means the same thing everywhere. A line qualifies when
it starts with `#`, an upper-case key, and `:`, which is what keeps prose like `# note: rerun weekly`
out. Scheduler directives are matched by prefix instead, because `#SBATCH --time=01:00:00` has no
`#KEY: value` shape to cut on.

The cost is deliberate: a recipe that writes a job script in a heredoc also declares whatever that
script declares. Both keys hash `catalog.StripComments`, which drops every whole-line comment —
annotations with them — so never compute that preimage a second way. Neither affects what is stored
or what runs: the recipe is embedded and executed byte for byte.

Overlays are stored as `.sqf` (SquashFS, read-only) or `.img` (ext3, writable).

## Base and Layering

A `base` provides the container root and the tooling a build and a run need — Apptainer and
Micromamba — and nothing an overlay's payload links against. Overlays mount `--overlay <path>:ro`
with payloads under their own `/cnt/<name>` prefix, so they are disjoint subtrees, not stacked
diffs. The only ordering rule is that the single writable `.img` goes last
(`putImgToLast`, `internal/runtime/container/setup.go`).

Because they are disjoint rather than stacked, **two images claiming one prefix do not combine** —
the later mount takes the whole subtree and the earlier contributes nothing, while both still reach
PATH and the environment. `ensureDistinctPrefixes` refuses that mount rather than warning: it is a
wrong container that looks like a working one. Two builds of one name are the case it catches. A
base, an OS image and anything without readable metadata record no prefix and are exempt.

**A base is never a compatibility gate.** It carries no identity or equivalence key, and any base
information an artifact records is build diagnostics only — never compared, never warned on, never
refused. A base is rebuilt whenever its OS ships patches, so a base digest differing from the one an
overlay was built against says nothing about whether they work together. Never add a check that
treats it as if it did.

## Publishing

What may be pushed depends on **who can pull it**, declared per endpoint as `audience` (`public` by
default, or `restricted`), and on what the artifact says about itself. Everything is read from the
embedded manifest, never guessed from the filename, and `push` refuses rather than warns.

A `restricted` endpoint takes anything. At a `public` one:

| the recipe declared | may publish publicly |
|---|---|
| `#REDISTRIBUTE: no` | **never**, whatever the type, and no flag overrides it |
| `#REDISTRIBUTE: yes` | yes, whatever the type |
| nothing, and it is `base` or `os` | yes — ours to publish: a container root, and packages from a public distribution |
| nothing, and it is `data` | yes — the type asserts public reference data and the indexes built from it |
| nothing, and it is an `app` | **no** — someone else's software with unstated terms, and unknown is not permission |
| nothing, and it is a Conda build | yes — see below |

**`#REDISTRIBUTE:` is the answer; the type is only the default for an unanswered question.** The
declaration lives in the recipe rather than on the command line, so it is authored once, reviewed in a
commit, and travels with every artifact built from it — the same standard `audience` is held to, and
the reason a `yes` is trusted here while a `--force` would not be. Nothing verifies it.

**A Conda build is never asked.** It embeds no recipe (`build.SourceSpec.RecipeFile`), so there is
nowhere for the declaration to be written that travels with the artifact, and gating it would be a
permanent refusal wearing a default's clothes. Its packages were also vetted for redistributable
licensing as a condition of being in the channel. What the channels cannot answer — a private or
vendor channel — is reported from `Build.Channels` at push time and decided by a person; do not add a
check that adjudicates it, since the only ways to try are a config allowlist or interpreting a few
hundred licence strings.

`#LICENSE:` is an SPDX expression published verbatim as `org.opencontainers.image.licenses`. It is
never parsed and never gates anything: deriving redistribution permission from a licence expression is
a judgement a tool gets wrong in the permissive direction, and that direction cannot be taken back.

`.img` is never published to any endpoint. That is structural, not policy — a writable overlay has no
identity, so there is nothing to publish it *as*.

`audience` states a fact about a registry and CondaTainer derives the permitted set from it. Never
let config enumerate types directly: `types: [app]` beside a public endpoint would erase the rule
with no error, whereas a wrong `audience` is a claim someone has to write down and defend. It is a
declaration, not enforcement — nothing verifies the registry is actually restricted, and it is
unrelated to a GitHub package's own visibility setting, which CondaTainer never reads or changes.

## Data Directory Order

Reads go nearest-first, writes furthest-first — opposite directions, same four tiers.

| | read (first match wins) | write (first writable) |
|---|---|---|
| 1 | `$SCRATCH/condatainer/` | `CNT_EXTRA_ROOT` (group/lab root, env only) |
| 2 | `~/.local/share/condatainer/` | `CNT_ROOT` / `<install-dir>/` (app-root, auto-detected) |
| 3 | `CNT_EXTRA_ROOT` | `$SCRATCH/condatainer/` |
| 4 | `CNT_ROOT` / `<install-dir>/` | `~/.local/share/condatainer/` |

A build lands as far out as permissions allow, so one copy serves the whole group;
anyone who wants their own version of a name builds it into their own directory and
has it win for them. Order *within* a tier is the same both ways.

Each contains `images/` and `helper-scripts/`.
Recipes are not searched here — they come from the ordered `sources` list.

## Helper Scripts

Bash scripts in [`cnt-scripts/helpers/`](https://github.com/Justype/cnt-scripts) launch interactive services inside CondaTainer on HPC. Modes: `headless/` (direct) and `<scheduler>/` (submit + SSH tunnel). See `helpers/README.md` for details.

## File Locking

`exec`/`run` hold `LOCK_SH` on `.sqf`/`.sif` files during execution (`.img` skipped — Apptainer flocks those itself); `remove` and `build --update` probe `LOCK_EX` before modifying. See `internal/image/lock.go`.

**An unwritable image is protected and is never modified or removed** — not even for its owner, who
can unlink it through the directory anyway and can restore the bit. Clearing the write bit
(`chmod a-w`) is how an artifact is pinned. This is why a write lock opens `O_RDWR`: `flock` does not
need it, so never "simplify" that to `O_RDONLY` — the open mode *is* the protection check.

A failed lock has three distinct causes and they must stay distinct: `ErrProtected` (write bit
clear), `ErrInUse` (a conflicting flock), and a missing file. Reporting a protected image as "in
use" sends the reader hunting for a container that is not running.

## Coding Rules

- When editing a function, update its doc comment to match. Keep comments concise and behavior-first — say what it does; give a reason only when the behavior is surprising.
- When editing code, check whether the change is reflected in docs and README files, and update them when needed: `docs/manuals/condatainer.md` and relevant `docs/` pages for UX changes (flags, output format, command behaviour), and the nearest `README.md` (e.g. `internal/helper/README.md`) for package-level changes.

## Key Patterns

- Global config singleton: `config.Global`
- Errors: **sentinels by default** — a package-level `ErrX`, wrapped with `%w`, matched by callers with
  `errors.Is`. That is what nearly every package does: `image.ErrProtected`/`ErrInUse`,
  `meta.ErrNoManifest`/`ErrUnsupportedSchema`/`ErrInvalid`, `tool.ErrFileNotFound`,
  `catalog.ErrNotProvided`, `registry.ErrIncompatible`/`ErrMismatch`, the `scheduler.Err*` set.
  Wrapping keeps the detail in the message and the category matchable in one value.
  Reach for a **struct** only when a caller must read a field rather than recognize a case —
  `apptainer.ApptainerError` (`ExitCode()`), `scheduler.SubmissionError`, `tool.Error`.
- Console output: `utils.PrintMessage`, `PrintWarning`, `PrintError`, `PrintDebug`
- Always use absolute paths; `config.Get*Dir()` for standard locations
- File/dir creation: use the wired helpers — `utils.CreateFileWritable` (files), `utils.MkdirAllShared` (dirs), `utils.MakeExecutable` (make a file executable). They apply `utils.PermFile`/`PermDir` (umask-subject) and call `ShareWithParentGroup` so children inherit group-write inside a shared `2775` install while personal installs stay umask-default. Never force perms with a bare `os.Chmod(path, utils.Perm*)` — see `internal/utils/files.go`.
