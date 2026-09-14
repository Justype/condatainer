# internal/libexec

Resolves and provisions CondaTainer's self-provisioned toolchain — `mksquashfs`, `squashfuse`,
and an ordinary (non-fakeroot) `apptainer` — installed via `micromamba` into one of the four
data-directory tiers. It is never a recipe's concern and never contributes to any artifact's
identity: the toolchain is host-local infrastructure, the same status as the Apptainer binary
itself.

## Placement

`libexec/` is a sibling of `images/` and `helper-scripts/` in each of the four data-directory
tiers, using the exact same nearest-read / furthest-write order (`internal/config`'s
`GetLibexecSearchPaths`/`GetWritableLibexecDir`). "One copy serves the whole group" applies at
least as much to a toolchain this size as it does to a single image.

Named `libexec`, not `tools`: this codebase already uses "tool" for other things (`internal/image/
tool`'s archive-read error shaping, `internal/toolpath`'s resolver, `Build.Tools`'s build-provenance
record). `libexec` is the FHS term for exactly this — executables invoked by other programs, never
typed by a user.

## Bootstrap sequence

The bootstrap binary must run from **outside** the prefix it creates. `micromamba create -r X -p
X` refuses outright if `X` already exists — even genuinely empty, not just non-empty — with
"Overwriting root prefix is not permitted." `writableTarget` resolves a tier's directory (which
`config.GetWritableLibexecDir` creates as a side effect, mirroring how an image's write directory
is resolved) and immediately removes it again so `create` sees a path that has never existed.

The create command is exactly:

```
<bootstrap-or-live-micromamba> -r <prefix> --no-rc create -y -p <prefix> -c conda-forge \
    micromamba squashfs-tools squashfuse apptainer
```

- `-r`/`--no-rc` are **global** flags and must precede the subcommand.
- `-r` pinned to the same path as `-p` stops libmamba's default `pkgs_dirs` list from falling
  back to a per-user `$HOME` location regardless of `-p` — a real leak, confirmed by watching it
  happen. On a shared, group-writable tier this is not cosmetic: an unpinned root silently
  defaulting to `$HOME` would break "one copy serves the whole group" for every other user
  hitting the same tier.
- `--no-rc` stops the invoking user's own `~/.condarc` (channels, `channel_priority`, …) from
  influencing the solve.
- The package is `squashfs-tools`, not `mksquashfs`.
- `micromamba` is included in its own install spec so the bootstrap binary is immediately
  superseded by a tracked, updatable package inside the prefix — the same self-management the
  official installer does.

After `create`, `cleanAt` runs `clean -a -f -y`. `-f`/`--force-pkgs-dirs` is required: `-a`/`--all`
alone leaves the per-package unpacked cache directories behind, because micromamba considers them
"in use" by the very environment they were installed into. Measured for the full toolchain: 2964
files / 325 MB with `-a` alone, 938 files / 293 MB with `-a -f`.

`provisionAt` then removes `conda-meta/` entirely and writes this generation's own `.lock`
sentinel (see Locking below). Removing `conda-meta/` is detection-hiding only, not a mutation
safeguard: tested both ways, `install`/`update -p <prefix>` "succeed" regardless — with
`conda-meta/` present it's a normal update, without it micromamba just treats the directory as
fresh and silently relinks over every existing file with warnings instead of errors. The only
effect of removing it is that `conda env list` and IDE Python-interpreter scanners (which walk
the filesystem for `conda-meta/history`) never pick up `libexec` as something to offer.

`provisionAt` also runs `utils.ShareTreeWithParentGroup(prefix)` before the lock sentinel. `create`
writes every file under `prefix` itself, entirely outside this codebase's own `MkdirAllShared`/
`CreateFileWritable` — confirmed against a real provisioned tree, where only `prefix`'s own top
level had ended up group-writable and everything beneath it (`bin/`, `lib/`, `etc/`, …) was still at
micromamba's own umask-derived mode. Left unfixed, `Update`'s later `os.RemoveAll(target)` would
only succeed for whoever happens to own each individual file — breaking "one copy serves the whole
group" for every group member except the one who provisioned this generation.

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

Only `conda-meta/` is stripped from a freshly provisioned prefix (see Bootstrap sequence above);
`activate.d`, `deactivate.d`, and `envvars` are all left in place, since any of them could matter
for the same reason `activate.d` does.

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
(`internal/runtime/apptainer.ResolveBin`), and it is a partial exception to "absolute path, not a
`PATH` prepend": Apptainer needs `unsquashfs`/`mksquashfs` in `PATH` for some of its own operations,
so `runApptainerWithOutput` also adds `BinDir()` to `PATH` for the Apptainer subprocess itself — a
different environment than the one a script running inside the container it launches sees (see
`internal/runtime/apptainer/README.md`, *Apptainer needs squashfs-tools in PATH*).

## Naming a missing tool

`Path`/`BinDir`/`Dir` build a path for any name regardless of whether this package actually installs
it — deliberately, since that is what lets `toolpath.Resolve` ask about a name like `debugfs` this
package never provisions and fall through cleanly. That means none of them can answer "is `name`
actually one of mine," which a caller needs before it decides what to tell the user. `Provides(name)`
answers it, against the exact *binary* names `provision.go`'s package spec installs — not the package
names themselves, since one package can give more than one binary: `squashfs-tools` gives two
(`mksquashfs`, `unsquashfs`), and so does `squashfuse` (`squashfuse`, `squashfuse_ll`).

`ErrNotProvisioned` and `NotProvisionedMessage` both report the same thing — nothing is provisioned,
or `Provides(name)` is false for what was asked — but reach different consumers. `ErrNotProvisioned`
is a Go sentinel for a caller that fails before starting a container (`apptainer.ResolveBin`,
`conda.go`'s `micromambaCmd`/`condaExecOpts`); `NotProvisionedMessage(name)` is a plain string for
anything else that needs the wording without a Go error to carry it — a host-side caller that already
resolved through `toolpath` and would otherwise need this package as a second import just for a
message. That case is why `toolpath.NotFoundMessage(name)` exists as a thin re-export:
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
must never depend on anything that depends back on one of its own children. The generation lock
below is the one thing this package used to share with `internal/image` (both wanted a plain
non-blocking flock); it now uses `internal/utils.AcquireFlock` instead of `internal/image.
AcquireLock` — the actual flock mechanics were the only genuinely shared part, so that's the one
piece that moved to a place both packages could reach without an import between them. `internal/
image.Lock` is now a type alias for `utils.FlockHandle` for exactly this reason: the two are the
same type, so a caller holding a mixed slice of overlay and toolchain locks (`internal/runtime/
exec.Prepare`) needs no conversion between what this package returns and what `internal/image`
returns.

## Locking

Every generation gets its own `.lock` sentinel, created fresh by `provisionAt` — never shared
across generations, since each is a distinct inode from its own `create` call. This mirrors
exactly how `.sqf`/`.sif` overlays are protected (CLAUDE.md, *File Locking*): a reader
(`AcquireUse`, called from `internal/runtime/exec.Prepare` alongside the existing overlay/base
locks) holds `LOCK_SH` for the duration of one apptainer subprocess call; `Update` must take
`LOCK_EX` on the live generation before touching anything, and refuses outright if it can't —
both sides use `internal/utils.AcquireFlock`'s non-blocking flock, so a reader that collides with
an in-progress update fails immediately with a clear message rather than blocking, and `Update`
never waits either.

Holding the lock guarantees nothing *condatainer* is using the outgoing generation — but `Update`
still doesn't remove it in place. It renames it to a fixed `target+".old"`, activates staging, then
makes one best-effort attempt to remove `old`. A rename succeeds even while a file inside is open or
executing; removal does not, so activation is never gated on cleanup.

`liveLock` is released right after that rename, not deferred to the end of `Update`: its own `.lock`
file moves with `old`, so holding it open through the removal attempt would mean unlinking a file
this process still has open — a guaranteed silly-rename on NFS. Releasing here is safe, since
nothing looks up a lock at `target` once staging has taken its place.

If removing `old` fails, `Update` logs a warning and leaves it — no retry loop, no crash-recovery
sweep. The next `Update` call clears it (or tries again) before renaming a fresh `target` aside; a
crash mid-`Update` recovers the same way, by running it again.

A generation with no `.lock` file yet (one written by code that predates this mechanism, or a
sentinel that failed to write) is self-healed by creating an empty one before locking, in both
`Update` and `AcquireUse` — otherwise every reader would report "being updated" forever for a
generation nothing is actually updating. Creating it without checking for a race first is fine:
the flock acquired right after is what actually serializes access, not the file's creation.

## Verification

`verifyToolchain` checks a staged (not yet live) generation by output content, not exit status.
`squashfuse` always exits 254 for an argument-free invocation — version and help flags included —
even though it prints its own name and version banner regardless; only a real archive +
mountpoint argument gets exit 0. `apptainer`/`mksquashfs`/`unsquashfs` do exit 0 for their own
version flag, but content-checking all four the same way needs no per-tool special case.

`versionFlag` is the one place that flag differs per tool: `mksquashfs`/`unsquashfs`
(squashfs-tools) predate GNU-style long options and only recognize `-version`; everything else
provisioned here takes `--version`.

Two floors, both checked by `verifyFloor` against `meetsFloor`:

- apptainer `>= 1.4` — below that, it cannot mount a zstd-compressed SquashFS. This package's own
  threshold is independent of `internal/runtime/apptainer.CheckZstdSupport` (same number) because
  verifying a *staged* binary must never touch that package's global, live-apptainer state.
- squashfs-tools (checked via `mksquashfs`, since `unsquashfs` always ships the same version)
  `>= 4.4` — below that, `mksquashfs` cannot produce zstd-compressed archives, and `unsquashfs`
  has no `-offset`, which `internal/image/squashfs` relies on to read a SquashFS partition inside
  a SIF in place.
