# internal/toolpath

Finds a runnable path for an external host tool by name — the one place in
the codebase that decision is made — and builds a ready-to-run `*exec.Cmd`
for it. Nothing else should resolve a bare tool name on its own.

## Why this is its own package

Three things used to be tangled together, in different packages, for no
reason other than history:

- **Where the self-provisioned toolchain lives and what it installs**
  (`internal/libexec`) — provisioning, locking, "do I have X."
- **Finding a runnable path for a name at all** — `PATH`, then the FHS
  fallback, with the self-provisioned copy checked first when it exists.
- **Image-domain error shaping** (`internal/image/tool`'s `Error`, its
  `analyze()` hint text, the archive sentinels) — genuinely specific to
  interpreting an `ext3`/`squashfs` tool's failure, and used only by
  `ext3`/`freeze`/`squashfs`.

Only the middle one is truly generic — `internal/conda` needs a resolved
`unsquashfs`/`debugfs` path and has no use for `image/tool`'s archive
sentinels or its `Error` shape, and `internal/image` (the parent package)
already imports `internal/image/squashfs`, so `image/squashfs` importing
anything that resolves back through `internal/image` would cycle. Splitting
lookup out into its own package, importable by everything without dragging
in a domain it doesn't need, is what makes it importable at all from
`internal/conda` and from `internal/image` and its children alike.

## Why the preference lives here, not in `internal/libexec`

`internal/libexec` only ever answers "do I have this, and where" (`Path`,
and the named `ApptainerPath`/`MksquashfsPath`/etc. built on it). It has no
opinion about `PATH`, the FHS fallback, or which source should win — that
would give a provisioning package an opinion about the host it has no
business holding. This package is what actually decides "which binary do we
run": `libexec.Path(name)` first — required, wins even over a same-named
binary already on `PATH`, and only when the tool is installed there — then
every `PATH` entry, then the tools bundled with the host apptainer
(`<LIBEXECDIR>/apptainer/bin`, read from `apptainer buildcfg` once per
process), then the FHS fallback directories. A host `mksquashfs` or
`unsquashfs` below the squashfs-tools floor (4.4) is skipped with the reason
kept for the error, so a good copy later in the order can still win.

Running those host binaries (`apptainer buildcfg`, a tool's version flag) is
remembered in a per-user file, `toolpath.json` in the personal cache directory,
so a new process does not repeat it. Each entry records the size and
modification time of the binary it describes and is ignored when they differ,
which is why an unloaded module, an upgrade, or another node's different
binary just misses. Inside a container there is no on-disk cache at all, since the same
path there can name a different binary than on the host. Config is not used for this: a config file can be a group
layer read on hosts with different modules. The cache is never needed for
correctness, and deleting it is harmless.
Nothing this package provisions (`mksquashfs`, `unsquashfs`, `squashfuse`,
`apptainer`) is special-cased; a name libexec never provisions (`debugfs`,
`e2fsck`, `tune2fs`, `mke2fs`, `resize2fs`, `fuse2fs`, all e2fsprogs, host-
only per `internal/image/README.md`) simply isn't found in `Path` and falls
straight through.

## `Resolve` vs. the container-bound accessors

`Resolve`/`Command` are for a **host-side** invocation — one with no bind
mount to construct: `internal/image/freeze`'s own packing (no Apptainer, no
container at all), an ordinary `.sqf` manifest read, `internal/conda`
reading `conda-meta` out of an image directly.

A **container-bound** call — packing inside a build, a conda install running
inside the exec container — needs the self-provisioned tool's *containing
directory* too, so `container.DeduplicateBindPaths` can bind it in and the
absolute path resolves inside the container the same way it does on the
host. That case bypasses this package entirely and calls
`libexec.ApptainerPath`/`MksquashfsPath`/`MicromambaPath` directly
(`internal/build/conda.go`, `squashfs.go`) — see `internal/libexec/README.md`,
*Resolved paths, not logical ones*.

## `Command`

`Command(ctx, name, args...)` is `Resolve` plus `exec.CommandContext` plus a
debug log line — the three steps that were duplicated, verbatim, inside
`internal/image/tool.RunCommand` and would be duplicated again by any other
domain's own run-and-wrap-the-error wrapper. It returns a bare `*exec.Cmd`
with no opinion about stdout/stderr wiring or how a failure should be
reported: the caller wires those and wraps the result in whatever error
shape its own domain needs (`image/tool.Error`, apptainer's
`ApptainerError`, ...). Only the mechanical "find it, then build the command"
step is shared — the error shape never is, because different callers
genuinely need different information back (an image caller wants a
filesystem-failure hint, an apptainer caller wants `ExitCode()`), and forcing
one shape on all of them would mean either losing information or every
caller doing its own type assertion out of a generic wrapper.

## What deliberately does not go through this package

Apptainer's own binary choice (`internal/runtime/apptainer`'s `Normal`, `Fakeroot` and `ForBuild`) is a
policy decision between exactly two named candidates — system/module for
fakeroot, `libexec`'s own otherwise — gated by fakeroot and a zstd version
floor, not a `PATH` walk; `Resolve` has nothing to offer it. Scheduler
(`sbatch`/`qsub`/`bsub`/`condor_submit`) and proxy (`ssh`/`loginctl`) tools
stay on their own resolution (`config.Global.Scheduler.Bin`, module detection)
for a real reason, not just narrower scope: `libexec` will never provision
them, so routing them through here would only ever hit the `PATH`/FHS
fallback — no behavioral difference, just an import for nothing.
