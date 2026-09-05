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
Four declarable types: `base` (produces the container root), `os`, `app` (contributes to `PATH`),
`data`. A fifth, `env`, is what `overlay freeze` captures from a writable overlay; no recipe may
declare it and `DeriveType` never returns it. It is written "environment" in anything a user reads,
since `env` elsewhere in the tool means environment variables.

A recipe runs top to bottom as `bash -euo pipefail <recipe>` — no `install()` wrapper.
Available vars: `$CNT_NAME` (the complete name, e.g. `samtools/1.23.1`), `$CNT_TYPE`,
`$CNT_PREFIX` (where the payload goes), `$CNT_TMP` (also `$TMPDIR`, both `/cnt_tmp`), plus the
scheduler's normalized `$NCPUS`, `$MEM`, `$MEM_GB`. There is no `$CNT_VERSION`: a recipe that
varies by version uses a `#PH:` placeholder, and one pinned to a version writes it literally.

`#DEP:` is a **build** dependency only — what must be mounted while the recipe runs. It is not
recorded in the image and never re-expanded at run time; there is no runtime dependency tree.
**Only a `data` recipe may declare one**, and **a build's `#DEP:` is a `name/version`, never a
path** — where a running script's may be either. Both rules live in `catalog.ValidateDeps`, and the
reasoning is in [`catalog/README.md`](catalog/README.md) with the `#TARGET:` naming rule an external
build depends on.

A **version constraint is a recipe-only feature**: it lets a build reuse a satisfying version already
installed. An **analysis script must name an exact `name/version`**, and a constrained declaration in
a project script is a scan finding, not a request.

Metadata headers: `#DEP:name/version`, or `#DEP:name/version>=min` in a recipe (preferred version is
implicit upper bound, so valid range is `[min, version]`), `#SBATCH`/`#PBS`/`#BSUB` (scheduler job params),
`#ENV:VAR={prefix}/sub  ## note` (env vars; `{prefix}` is filled with the install prefix at load time),
`#INPUT:prompt` (user input, fed on stdin in order — read with `IFS= read -r VAR`), `#PH:`/`#TARGET:` (templates),
`#ARCH:noarch` (app and data script recipes only; default `native`), `#DESC:`, `#URL:`, `#TYPE:`,
`#LICENSE:` (an SPDX expression, verbatim), `#REDISTRIBUTE:` (`yes` or `no` — see *Publishing*).

**`#TARGET:` without `#PH:` is not a template — it names the artifact**, which is how an external
build (`create -p <path> -f <script>`) gets a name at all; one declaring `#DEP:` must declare it.

`catalog.ScanAnnotations` is the one tokenizer for every `#KEY: value` line, in recipes, user scripts
and project scanning alike; position carries no meaning, so there is no header block. Both keys hash
`catalog.StripComments`, so never compute that preimage a second way. See
[`catalog/README.md`](catalog/README.md), *The header boundary*.

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

What may be pushed depends on **who can pull it** — declared per endpoint as `audience` (`public` by
default, or `restricted`) — and on what the artifact says about itself. Everything is read from the
embedded manifest, never guessed from the filename; `push` refuses rather than warns, and no flag
overrides a refusal.

What a public endpoint takes is the table in [`docs/manuals/condatainer.md`](docs/manuals/condatainer.md),
*Publishing rules*. Why each row falls that way, and the three checks that must never be added, are in
[`internal/registry/README.md`](internal/registry/README.md).

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

`exec`/`run` hold `LOCK_SH` on `.sqf`/`.sif` files during execution (`.img` skipped — Apptainer flocks
those itself); `remove` and `build --update` probe `LOCK_EX` before modifying.

**An unwritable image is protected and is never modified or removed.** Clearing the write bit
(`chmod a-w`) is how an artifact is pinned, a write lock opens `O_RDWR` because that open mode *is*
the check — never "simplify" it to `O_RDONLY` — and a failed lock's three causes stay distinct:
`ErrProtected`, `ErrInUse`, and a missing file. See
[`internal/image/README.md`](internal/image/README.md).

## Coding Rules

- When editing a function, update its doc comment to match.
- **Before changing behavior, read the package's `README.md`.** It records the decisions the code
  cannot show and the orderings that are load bearing. A change that contradicts one is either wrong,
  or right and the README is updated in the same change — never left to disagree with the code.
- When editing code, check whether the change is reflected in docs and README files, and update them when needed: `docs/manuals/condatainer.md` and relevant `docs/` pages for UX changes (flags, output format, command behaviour), and the nearest `README.md` (e.g. `internal/helper/README.md`) for package-level changes.

### Where writing goes

Three places, three jobs. Writing put in the wrong one is not a style problem: a
design argument in the manual is noise to the person reading it, and a paragraph
defending a decision beside the code is read as documentation of behavior.

- **Comments — what the code does.** Concise and behavior-first. Give a reason
  only where the behavior would look wrong without it, and then in a clause, not
  a paragraph. **A function comment over 6 lines is a mistake** unless the
  function branches and each branch decides something; the rest belongs in a
  README. Never argue with an objection nobody made at the call site.
- **`docs/` and CLI help — what a user does and sees.** End-user facing:
  commands, flags, output, what is refused and how to proceed. This covers every
  `Short`, `Long` and `Example` string as much as it covers `docs/`. **Never a
  design argument, never why a decision went one way, and never how it works
  inside** — no "reads the lock rather than scanning", no "one extraction rather
  than one read per file". State the rule and the consequence the user can act
  on, not its defence.
- **A package `README.md` — the design.** The decisions, the constraints they
  answer, the orderings that are load bearing, and anything true of the package
  that no single file states. **Not a retelling of what the code already says.**
  Reasoning cut from a comment or from the manual lands here.

**Every document describes the present, in all three places.** When behavior
changes, **edit the part that is now wrong** — never add a paragraph beside it
saying what it used to do or why that was wrong. No "previously", no "note that
this no longer", no correction layered on a correction. The reader has not seen
the old version and git holds it anyway. A section that has grown by accretion is
rewritten, not appended to: **documents must not grow like a tumor.**

### Tests describe the present too

The same rule, applied to test code. A test names behavior that exists now, so
when behavior changes **delete the tests for what it used to do** rather than
adapting them or adding new ones beside them.

**A removed special case takes its tests with it.** Do not keep a test alive to
prove the old rule no longer applies — nothing is left to apply it. When a
special case is replaced by a general rule, the replacement gets **one** test,
not one per case the old rule enumerated.

**Test code must not grow like a tumor either.** A change that removes behavior
should leave fewer test lines than it found.

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
