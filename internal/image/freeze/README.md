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
mount.go      mountedRun: the one FUSE-mount primitive every mount site uses
fuse2fs.go    Resolve fuse2fs/squashfuse on PATH
tools.go      Build-tool provenance recorded in the manifest
whiteout.go   §2.4a whiteout translation
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
`mountedRun`, and it needs no container either.

## `mountedRun`: the mount, without Apptainer

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
run leaks nothing to clean up.

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

Opaque-directory translation (§2.4a) needs to know what a directory held in
the base to turn "this directory was replaced wholesale" into one whiteout per
entry that base had. The base is a plain `.sqf`, so this is answered with
`unsquashfs -l -d "" <base> <dirs...>` — one call, however many directories are
asked about, grouping every returned line by its immediate parent to recover
"the direct children of X" for each. A directory the base doesn't have
produces nothing, exactly like one that exists and is empty: neither hides
anything, so the two cases don't need to be told apart.
