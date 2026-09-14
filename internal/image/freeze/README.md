# freeze

Converts a writable ext3 overlay into an immutable SquashFS artifact
(`Freeze`), and back (`Unfreeze`). Neither needs Apptainer or a base image to
run — both read and write images directly on the host.

## Architecture

```
freeze.go     Freeze: walk, translate deletions, pack, identify, embed metadata
pack.go       Pack/AppendMeta: the mksquashfs runs
identity.go   TreeIdentity: hash the packed payload
unfreeze.go   Unfreeze: rebuild a writable image from an artifact
walk.go       Walk: list an overlay's upper/ via debugfs, no mount
dump.go       dumpUpper: copy upper/ out via debugfs rdump, no mount
mount.go      MountedRun: the one FUSE-mount primitive every mount site uses
sentinel.go   RunSentinel: the process MountedRun re-execs into to survive its own death
fuse2fs.go    Resolve fuse2fs/squashfuse on PATH
tools.go      Build-tool provenance recorded in the manifest
whiteout.go   Whiteout translation
```

## No Apptainer, no base image

`mksquashfs`, `mke2fs`, `debugfs` and `unsquashfs` all run directly on the
host, the same as every other package under `internal/image` already does —
there was never a reason `mksquashfs` alone should be the one tool assumed
absent from the host. Reading an `.img`/`.sqf` without mounting (`Walk`,
`dumpUpper`, `BaseDirLister`) needs nothing further.

The one place a mount is unavoidable is reading a live payload during `Pack`'s
default route, `TreeIdentity`, and `Unfreeze`'s `buildImage` — a mounted
SquashFS is 2.4x cheaper to read than a mounted `.img`'s `fuse2fs`, and
`mke2fs -d` needs a live mount specifically because it can discover whiteout
device nodes via `stat()` that an unprivileged extraction (`unsquashfs -d`)
cannot create at all (see `Unfreeze`'s doc comment). That mount is `mount.go`'s
`MountedRun`, and it needs no container either. Exported — `internal/build`'s
own SquashFS packer (`squashfs.go`'s `packFromScratchImage`) reuses it
directly to read a build's scratch `.img`, the same apptainer-free way, so
mksquashfs's package-time reads never need a container any more than this
package's own do.

## `Build.Tools` records only what actually ran

A frozen artifact's manifest never carries an `Apptainer` entry: freeze never
runs it, on either pack route, so there is nothing to record — not even an
empty `Tool{}` standing in for "didn't run" (`ValidateManifest` rejects a
snapshot manifest that sets one). Micromamba is absent for the same reason
`tools.go`'s own doc comment gives: a freeze packs an environment that
already exists rather than solving one.

`Mksquashfs` and, on the mount route only, `Fuse2fs` come from `Pack` itself
(`pack.go`), captured at the point each binary is resolved rather than
re-resolved by the caller. `Pack` returns them, and `Freeze` (`freeze.go`)
merges them with `buildTools`'s own `Condatainer` field to build the
manifest — which is why `Pack` runs *before* the manifest is built, not
after: the manifest cannot record a tool `Pack` hasn't resolved yet.

## `MountedRun`: the mount, without Apptainer

A bare, unprivileged FUSE mount can be refused outright (`Operation not
permitted`) depending on the host: `fusermount3`'s usual escalation path is a
**setuid-root** binary, and that path is unavailable wherever the setuid bit
doesn't apply — a filesystem mounted `nosuid` (which is exactly how every
container, Apptainer's own included, mounts its root, as a deliberate hardening
measure), a `fusermount3` that was never installed setuid, or an admin policy
disabling it. Apptainer's own `--fusemount` isn't exempt from this: on every
install checked, Apptainer carries no setuid component at all, so its FUSE
mounts (including mounting a bare `.sqf` as a container root in the first
place) succeed by a *different* mechanism — one available to any unprivileged
process, not something special to Apptainer.

That mechanism is `unshare --mount --user --map-root-user`: creating a private
mount+user namespace and mapping namespace-uid 0 to the real caller. The kernel
grants full capability — including `CAP_SYS_ADMIN`, enough for `mount()` — to
whoever is uid 0 *within a namespace nobody else can see*, with no setuid
binary involved at all. Because it doesn't touch the setuid path, it isn't
blocked by `nosuid`, and it works identically whether it's invoked from a bare
login shell or from inside an already-running (even Apptainer) container —
each invocation creates its own fresh namespace rather than inheriting
privilege from whatever it's running inside of. This is the same reason nested
`apptainer exec` works at arbitrary depth: namespaces all the way down, never
inherited.

**The mount is ended by killing the foregrounded FUSE process, never
`umount`/`fusermount`.** Those are the same setuid-root binaries described
above, and a setuid binary's *owning* uid — real root — has no mapping inside
this namespace, so some kernels refuse to even exec it rather than silently
ignoring the bit. Running the FUSE tool with `-f` (foreground) and killing it
tears down its own session directly; the whole namespace, mount included,
disappears the moment nothing is left running in it, so a killed or crashed
run leaks nothing to clean up — as long as something is still alive to do the
killing. That's what the sentinel below is for.

## The sentinel: surviving condatainer itself dying mid-mount

`MountedRun` doesn't call `unshare` directly. It re-execs itself into the
hidden `_mount_sentinel` command (`sentinel.go`'s `RunSentinel`), which then
runs the `unshare`/bash/FUSE tree as its own child, joined into its own
process group.

The reason is `Pdeathsig` (`prctl(PR_SET_PDEATHSIG)`): a process can ask the
kernel to send it a signal the instant its parent dies, for any reason —
SIGKILL, an OOM-kill, a crash — with no polling needed. That would be the
obvious fix for condatainer dying with a mount still open: set it on the
`unshare` process, catch it, kill the group. It doesn't work. Entering the
user namespace clears `Pdeathsig` on whatever process does it — the same
kernel rule that clears it across a setuid exec, since becoming
root-in-namespace is exactly that kind of privilege-elevating credential
change. Confirmed directly against this project's target kernels, not just
read off the man page: a process with `Pdeathsig` set that then runs `unshare
--user` never receives the signal when its parent dies.

So the one process that can reliably be told "your supervisor just died" has
to stay outside the namespace. That's the sentinel: it never escalates
privilege itself, so its own `Pdeathsig` (`SIGTERM`, set by `MountedRun`)
keeps working normally. `unshare`/bash/FUSE are joined into the sentinel's
process group rather than starting one of their own, so when the sentinel's
trap fires, one `kill(-pgid, SIGKILL)` takes the whole tree down with it —
the same call `MountedRun`'s own `cmd.Cancel` already made for the
still-alive-and-choosing-to-cancel case, just triggered from the inside
instead of the outside.

This narrows the leak window, it doesn't close it: if the sentinel itself is
killed directly (not condatainer — the sentinel, by its own pid), it can't
run its trap either, and `unshare`/bash/FUSE orphan exactly as they did
before this existed. What changes is *what* has to die to cause that — a
specific, tiny, short-lived process, rather than condatainer itself for the
entire span of a mount.

## Unfreeze can't get correct ownership from the mount — only after it

`unshare --map-root-user` grants exactly one working identity inside its
namespace: namespace-uid 0. Any other value asserted via `squashfuse -o
uid=<real-uid>` isn't a mapped value in that single-entry namespace, and the
kernel clamps it to the overflow uid (65534) rather than the one asked for.
Widening the mapping to include a second, real identity would need
`newuidmap`/`/etc/subuid` delegation — not guaranteed to exist for an arbitrary
user on an arbitrary cluster, and depending on it would reintroduce the exact
"works here, not there" fragility dropping Apptainer was meant to remove.

So `buildImage` doesn't try: it mounts the artifact with no `-o uid=`
override at all, `mke2fs -d` bakes in whatever raw ownership the `-all-root`
archive already has (`0:0`), and `ext3.ChownRecursively` — the same tool
`overlay chown` already uses — fixes it afterward, on the plain host, where no
namespace mapping applies to get in the way. This is confirmed correct against
a real `-all-root` artifact run through the actual mount pipeline, not just a
synthetic one.

## `BaseDirLister` reads the base directly

Opaque-directory translation needs to know what a directory held in
the base to turn "this directory was replaced wholesale" into one whiteout per
entry that base had. The base is a plain `.sqf`, so this is answered with
`unsquashfs -l -d "" <base> <dirs...>` — one call, however many directories are
asked about, grouping every returned line by its immediate parent to recover
"the direct children of X" for each. A directory the base doesn't have
produces nothing, exactly like one that exists and is empty: neither hides
anything, so the two cases don't need to be told apart.
