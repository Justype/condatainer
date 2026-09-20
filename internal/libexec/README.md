# internal/libexec

Resolves and provisions CondaTainer's self-provisioned toolchain: one micromamba prefix in one of
the four data-directory tiers. It always holds `micromamba`, and holds `squashfs-tools`
(`mksquashfs`, `unsquashfs`), `squashfuse` and an ordinary (non-fakeroot) `apptainer` only when
asked for by name. It never contributes to any artifact's identity: the toolchain is host-local
infrastructure, the same status as the Apptainer binary itself. A Conda or script build creates it with `micromamba` alone when no tier has one
(`EnsureMicromamba`). A script build can call
`micromamba` from it, appended to `PATH` so anything the base provides wins; that is a convenience,
not an input, so the version is not recorded and a recipe that solves with it is identified by its
text.

## Placement

`libexec/` is a sibling of `images/` and `helper-scripts/` in each of the four data-directory
tiers, using the exact same nearest-read / furthest-write order (`internal/config`'s
`GetLibexecSearchPaths`/`GetWritableLibexecDir`). "One copy serves the whole group" applies at
least as much to a toolchain this size as it does to a single image.

Named `libexec`, not `tools`: this codebase already uses "tool" for other things (`internal/image/
tool`'s archive-read error shaping, `internal/toolpath`'s resolver, `Build.Tools`'s build-provenance
record). `libexec` is the FHS term for exactly this — executables invoked by other programs, never
typed by a user.

## The package table

`packages` in `provision.go` is the one list of what a prefix may hold: each package's conda name,
the binaries it installs, which of them are checked for running, which one is reported and
floor-checked, and the floor. `Provides`, `Installed`, `Versions`, `verifyToolchain` and every
message that names a package derive from it. What a prefix holds is read from its `bin/`, not from
code, so a prefix with only `micromamba` is a complete, valid one.

`micromamba` is the base: its binary is what `Dir()` treats as "this tier is provisioned", and it is
what installs and updates everything else. Nothing else is installed unless named.

## Create and update

Each `update --libexec` call takes an exclusive lock on the tier (a sibling `.libexec.lock`, so it
exists before any prefix does) and, when a prefix is live, on that prefix's own `.lock` as well.
Then one of three things happens.

- **No prefix yet.** A standalone micromamba is downloaded outside the tier and creates the prefix
  directly at its final `libexec/` path with `micromamba` plus any named packages. A failed create
  removes the prefix, so nothing half-built is ever reported as provisioned.
- **A live prefix with `conda-meta/`.** The prefix's own micromamba updates it in place: named
  packages that are missing are installed, then the named ones (or, with none named, every
  installed package) are updated. `update`, not `install`, because `install` leaves an
  already-satisfied spec at its current version.
- **A live prefix with no `conda-meta/`.** micromamba cannot track such a prefix; it would treat it
  as empty and forget every installed package. It is removed and recreated with the tools that were
  in its `bin/`, plus any named.

Afterwards, whichever path ran: `fusermount3` is linked into `bin/` if only `sbin/` has it,
`clean -a -f -y` reclaims the package cache, the tree is shared with the parent group, the lock
sentinel is ensured, and `verifyToolchain` checks what is installed. A failed update leaves whatever
micromamba left, and re-running the update repairs it; there is no rollback copy.

`micromamba self-update` is not used: it swaps the binary but leaves `conda-meta/` recording the old
version.

### Nothing leaks out of the prefix

`-r <prefix>` alone does not contain micromamba. `runMicromamba` also sets
`CONDA_PKGS_DIRS=<prefix>/pkgs` and `XDG_CACHE_HOME=<prefix>/.cache`, and points `HOME` at a
throwaway directory removed afterwards. Without them micromamba writes to `~/.mamba/pkgs`,
`~/.cache/conda`, and registers the prefix in `~/.conda/environments.txt`, which no variable turns
off. On a shared tier that would also put the cache in whichever user ran the update.

A pre-existing directory at the prefix makes `create` refuse, so the throwaway `HOME` cannot live
inside it. `clean -a -f` is required: `-a` alone leaves the per-package unpacked directories,
because micromamba considers them in use by the very environment they were installed into.

### Why `conda-meta/` stays

Removing it hides the prefix from `conda env list` and IDE scanners, but `install` and `update` then
"succeed" while forgetting every earlier package record, leaving those packages' files untracked.
The cost of keeping it is that the prefix can be offered as an environment.

### Baked-in absolute paths

micromamba's prefix substitution patches some absolute paths into the binaries it installs, using
the `-p` prefix `create` was given — not just library `RPATH`s (which use `$ORIGIN` and survive a
move) but paths a binary execs at runtime. `libfuse3`'s `fusermount3` helper lookup is one: an
absolute `<prefix>/bin/fusermount3`, tried before any `$PATH` search. So the prefix is created at
the path it keeps and never renamed, and `ensureFusermount3InBin` links `bin/fusermount3` to
`../sbin/fusermount3` when the package installed it only in `sbin/`. The path handed to `-p` is
the symlink-resolved one, since that string is what gets baked in.

## Toolchain activation

Every accessor this package exposes (`Path`, `ApptainerPath`, `MicromambaPath`, …) resolves a
binary's absolute path for a caller to invoke directly rather than an environment to `conda
activate` — but the self-provisioned `apptainer` is still run through a wrapper
(`internal/runtime/apptainer`) that sources this prefix's own `etc/conda/activate.d/*.sh` first,
the same environment a real `conda activate` of this prefix would produce. This is required, not
incidental: apptainer mounting a large overlay hung indefinitely when invoked unactivated, and
ran correctly once its own `activate.d` had run first. `CONDA_PREFIX` is set only as a per-command
prefix on each script's own `.` (`CONDA_PREFIX=<prefix> . "$script"`), never exported for the whole
wrapper, so it reaches that script and anything it spawns but not the containerized command
apptainer goes on to run — see `internal/runtime/apptainer/README.md`.

`conda-meta/`, `activate.d`, `deactivate.d`, and `envvars` are all left in place, since any of
them could matter for the same reason `activate.d` does.

## Resolved paths, not logical ones

`Dir()` (and everything built on it — `BinDir`, `Path`, `ApptainerPath`, `MksquashfsPath`,
`UnsquashfsPath`, `SquashfusePath`, `MicromambaPath`, `LockPath`) resolves symlinks before
returning. Every caller embeds the result
as a literal path in a script or command line run inside a container, and
`container.DeduplicateBindPaths` resolves symlinks before constructing the actual `--bind`
argument, using the resolved path on both the host and container side. `$SCRATCH` pointing at
separate real storage is common on HPC, so a caller embedding the *logical* config path would name
something that may not exist inside the container at all.

Callers use the resolved absolute path directly as the command to run
(`micromambaCmd()` in `internal/build/conda.go`), not a `PATH` prepend. Two independent reasons:
Apptainer's `--env` replaces a variable outright rather than merging into it, so a `PATH` prepend
would have to happen inside the script text regardless; and an absolute path cannot be silently
shadowed by some other same-named tool earlier in the container's own `PATH`; a `PATH` prepend
can, with no error at all.

That is the pattern for a *container-bound* call — one that also needs the containing directory
bound in, so the absolute path resolves inside the container too (`internal/build/conda.go`'s conda
install). A *host-side* call — freeze/unfreeze's own packing, `internal/build/squashfs.go`'s
packing (never container-bound at all, any workspace mode — see `internal/build/README.md`, *The
container root*), an ordinary `.sqf` read — has no bind mount to construct, so it only needs a
runnable path, and reaches one a different way: see below.

`ApptainerPath()` is the one remaining container-bound accessor with a caller
(`internal/runtime/apptainer.Normal`), and it is a partial exception to "absolute path, not a
`PATH` prepend": Apptainer needs `unsquashfs`/`mksquashfs` in `PATH` for some of its own operations,
so `runApptainerWithOutput` also adds `BinDir()` to `PATH` for the Apptainer subprocess itself — a
different environment than the one a script running inside the container it launches sees (see
`internal/runtime/apptainer/README.md`, *Apptainer needs squashfs-tools in PATH*).

## Naming a missing tool

`Path(name)` returns a path only when that binary is installed in the provisioned `bin/`, so
`toolpath.Resolve` asking about a name like `debugfs` this package never provisions falls through
cleanly, and so does one that is provisionable but not installed. Neither tells a caller *why*, which
it needs before deciding what to say. `Provides(name)` answers whether `name` is a binary or package
the table can install, and `Installed(name)` whether this tier has it. A package can give more than
one binary: `squashfs-tools` gives `mksquashfs` and `unsquashfs`, and `squashfuse` gives `squashfuse`
and `squashfuse_ll`.

`ErrNotProvisioned` and `NotProvisionedMessage` report a missing tool to different consumers.
`ErrNotProvisioned` is the Go sentinel for a caller that fails before starting a container
(`apptainer.Normal`, `conda.go`'s `micromambaCmd`/`condaExecOpts`); `NotInstalledError(name)`
returns it when no tier is provisioned and otherwise an error naming the package that installs
`name`, since the bare `update --libexec` no longer installs it. `NotProvisionedMessage(name)` is a
plain string for anything else that needs the wording without a Go error to carry it — a host-side
caller that already resolved through `toolpath` and would otherwise need this package as a second
import just for a message. That case is why `toolpath.NotFoundMessage(name)` exists as a thin re-export:
`internal/image/tool`'s `CheckDependencies`, `internal/image/squashfs`'s cached-empty-string case, and
`internal/image/freeze/fuse2fs.go`'s `findSquashfuse`/`findFuse2fs` all go through that, not this
package directly — every one of them already imports `toolpath` for `Resolve` itself, so reaching past
it into `libexec` too would leak this package's own provisioning knowledge into callers that only ever
needed "why wasn't it found," which `toolpath` already owns. `NotProvisionedMessage`'s text is
deliberately single-quoted rather than backtick-quoted like `ErrNotProvisioned`'s own: it is meant to
sit inside a double-quoted shell `echo`, where a backtick would attempt command substitution instead
of printing literally.

## This package is standalone — no dependency on `internal/image`

`Path(name) (string, bool)` names where a binary called `name` would sit in the provisioned
`bin/`, and whether a toolchain is provisioned at all — regardless of whether `name` is actually
one of the tools this package installs. That is deliberately as far as this package's own opinion
goes: it does not know `PATH` or the FHS fallback directories exist, and it does not decide "which
binary wins."

`internal/toolpath` is the package that decides that, for the whole codebase: `toolpath.Resolve`
imports this package directly, checks `Path(name)` first (required, wins even over a same-named
binary already on `PATH`), then falls back to `PATH` and FHS. Every host-side consumer —
`internal/image/squashfs`, `internal/image/freeze`, `internal/conda` — calls `toolpath.Resolve`,
never this package's `Path` directly. See `internal/toolpath`'s own doc comment for why the
combinator lives there and not here: its consumers aren't only this package's own domain (`internal/
conda` isn't under `internal/image` at all), and folding "which binary wins" into this package would
give it an opinion about hosts and PATH search that provisioning has no business holding.

This package importing nothing under `internal/image` is not incidental — it is what lets
`internal/toolpath` import this package directly with no cycle, since `internal/image/squashfs`
(which needs `toolpath.Resolve`) is itself reachable from `internal/image`, and `internal/image`
must never depend on anything that depends back on one of its own children. The lock
below is the one thing this package used to share with `internal/image` (both wanted a plain
non-blocking lock); it now uses `internal/utils.AcquireFileLock` instead of `internal/image.
AcquireLock` — the actual locking mechanics were the only genuinely shared part, so that's the one
piece that moved to a place both packages could reach without an import between them. `internal/
image.Lock` is now a type alias for `utils.FileLock` for exactly this reason: the two are the
same type, so a caller holding a mixed slice of overlay and toolchain locks (`internal/runtime/
exec.Prepare`) needs no conversion between what this package returns and what `internal/image`
returns.

## Locking

Two locks, both non-blocking file locks through `internal/utils.AcquireFileLock`, so a collision fails at
once with a clear message and nobody waits.

- **`<prefix>/.lock`.** A reader (`AcquireUse`, called from `internal/runtime/exec.Prepare` alongside
  the existing overlay/base locks) holds a shared lock for the duration of one apptainer subprocess
  call. `Update` takes an exclusive lock before touching a live prefix and refuses if a reader holds it. This
  mirrors how `.sqf`/`.sif` overlays are protected (CLAUDE.md, *File Locking*). A prefix with no
  `.lock` file gets an empty one before locking: creating it without checking for a race is fine,
  since the lock acquired right after is what serializes access.
- **`<tier>/.libexec.lock`.** Held by `Update` for its whole run, so two updates of one tier cannot
  interleave, including the first one, when there is no prefix yet to hold a lock. When a prefix
  without `conda-meta/` is recreated, the in-prefix lock is released before the removal, because the
  file goes with the directory and unlinking a file this process still has open leaves a silly-rename
  on NFS.

## Verification

`verifyToolchain` checks a prefix's installed packages by output content, not exit status.
`squashfuse` always exits 254 for an argument-free invocation — version and help flags included —
even though it prints its own name and version banner regardless; only a real archive +
mountpoint argument gets exit 0. `apptainer`/`mksquashfs`/`unsquashfs` do exit 0 for their own
version flag, but content-checking all of them the same way needs no per-tool special case.
`micromamba` prints only a bare version, so it is not content-checked.

`versionFlag` is the one place that flag differs per tool: `mksquashfs`/`unsquashfs`
(squashfs-tools) predate GNU-style long options and only recognize `-version`; everything else
takes `--version`.

Two floors, in the package table and checked only for packages that are installed:

- apptainer `>= 1.4` — below that, it cannot mount a zstd-compressed SquashFS. This package's own
  threshold is independent of `internal/runtime/apptainer.CheckZstdSupport` (same number) because
  verifying a prefix must never touch that package's global, live-apptainer state.
- squashfs-tools (checked via `mksquashfs`, since `unsquashfs` always ships the same version)
  `>= 4.4` — below that, `mksquashfs` cannot produce zstd-compressed archives, and `unsquashfs`
  has no `-offset`, which `internal/image/squashfs` relies on to read a SquashFS partition inside
  a SIF in place.

A prefix that fails verification after an update is reported as such, and re-running the update
repairs it; `Dir()` still treats it as provisioned, since its `micromamba` is what runs the repair.
