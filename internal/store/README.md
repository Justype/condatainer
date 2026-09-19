# Overlay store

The store is immutable overflow below each images root. Flat overlays keep the
only bare-name slot; store entries are addressed by artifact name plus a
verified complete identity.

`Begin` first adopts an exact flat or stored artifact when one is already
available. Otherwise it picks the destination: the **flat** bare name when no
readable root holds it, and a collision-safe **store** filename only when some
root already holds that name at a different identity. A flat artifact answers
to `exec -o`, to `list`, and to every other project, so taking a free name
keeps one identity from being rebuilt once per project; the store is overflow
for the conflict case. `store/` is not created unless a conflict sends
something there.

Occupancy is judged across every readable root rather than the destination
alone. Reads are nearest-first while writes go furthest-first, so a flat copy
in a nearer root would shadow a new flat install and leave two identities
answering to one name. The decision is re-checked under the target's producer
lock, and a writer that loses the bare name between the scan and the lock
overflows instead of replacing it.

`BeginOptions.StoreOnly` skips the bare name even when it is free, for a caller
that has already promised it to something else — restore does this for a build
dependency sharing a name with a selection. It decides where a new artifact is
*written*, not what answers to the name: an exact copy already at the flat name
is still adopted.

Having chosen a target, `Begin` acquires the existing producer lock for that
specific path. There is no store-wide or store-directory lock. The producer
writes directly to the prepared sibling returned by the transaction, on the
target filesystem.

`Commit` syncs the prepared file, regenerates its identity and equivalence keys,
checks its name and expected keys, rechecks the target, and publishes with the
same atomic sibling rename used by `BuildObject`. Publication is create-only in
both layouts: an exact target is adopted and a conflicting target is never
replaced. `Abort` removes only that producer's prepared output and releases only
that target's producer lock.

`Detach` hands a reserved target to a producer that outlives this process — a
scheduler job — closing the transaction without releasing the lock. The lock is
rewritten with that job's `producer.Info`, which is what holds the pathname
across the queue wait: the job adopts its own lock by job ID, publishes through
its own transaction, and releases it there. A job that never runs leaves a lock
that goes stale with it, which the next producer clears. `Reserved` reports
whether there is anything to hand over — false means an exact copy was adopted
and nothing has to be produced.

Producer locks coordinate creation because a missing target has no inode to
lock. Once an SQF exists, `exec` and `run` hold a shared inode lock while
reading it; operations that replace or remove an existing file require an
exclusive inode lock.

## Placement and promotion

`InstallFile` is the entry point for an artifact that already exists as a file —
`store add`, and `project restore` for a downloaded prebuilt. The source filename
is read as nothing: the destination name comes from the keys in the file, because
inside the store the filename *is* the address and a name taken from whatever the
sender happened to call it would address a different artifact. It runs the ordinary
transaction, so an exact copy anywhere is adopted rather than copied again — which
callers naming a destination narrow with `SearchDirs`, since "it exists somewhere
else" does not answer a request to have it *here*. It
copies and never moves or links: the source is usually on another filesystem, a
hardlink would make `GC` report freeing bytes a second link still holds and would
share mode bits and locks with the user's own file, and a symlink is refused by
`Commit` outright.

A build reaches the store through `BuildObject.publishToStore`, which reads the
keys back off the packed output and hands the prepared file to a transaction —
one rename, since the build's prepared sibling and the store target are in the
same images root. It always passes `StoreOnly`. That is not a policy choice: the
build holds the producer lock on the flat target for its whole duration
(`Target.Lock` is `Path + ".lock"`), which is the lock `reserveFlat` would take,
so considering flat placement would deadlock the build against its own lock.

`Promote` decides which identity a bare name resolves to. It only renames, and
only inside the candidate's own root: the incumbent goes to that root's `store/`
under its own regenerated identity, and the candidate takes the bare name, both
under the bare name's producer lock so no build can claim it in between. The
demotion completes before the promotion begins, so an interruption leaves both
artifacts in `store/` with the name free — a state re-running repairs, and one no
`.part` sweep can eat.

Confinement to one root is what makes it safe for everyone else. Nothing leaves
the root, so no identity disappears from any reader's set whatever their layer
configuration, and a project lock — which pins name and keys, never a path —
still resolves both. A candidate in a root that a nearer root shadows is refused
with `ErrShadowed` rather than renamed into place, because reads are
nearest-first and the promotion would change nothing; a candidate *nearer* than
the current holder simply wins, and the farther copy is reported as shadowed and
left alone.

The artifact cache is never updated by either. `Lookup` requires the stored
fingerprint to equal a fresh `lstat` over device, inode, size, mtime and ctime,
so a renamed path misses and is recomputed. It could not work anyway: the cache
is per-user, so a promotion in a shared root cannot reach any other reader's copy.
The `Forget` calls are housekeeping for the caller's own file.

## Removal and collection

`Remove` deletes one store entry. Only the store layout — a flat artifact
answers to a bare name and belongs to `remove`, and taking a name out of service
through a command that reads as identity housekeeping would be a surprise.

`GC` reports what the selected stores could give back, and with `Apply` set,
takes it. Report and apply are one traversal rather than two phases, which is
what closes the gap between deciding and deleting: the exclusive lock taken to
judge an entry is still held when it is unlinked. A later apply run repeats the
whole judgement, because the report it was shown may be hours old.

An entry is collectable when its exclusive lock is free, it is older than the
grace, and it validated during the scan. Every uncertainty retains: an entry
that cannot be opened, locked, stat'd or validated stays, because a file that
cannot be shown to be garbage is not garbage.

Age is `max(atime, mtime, ctime)`, and the report names which one decided. atime
is the signal that tracks use and is what HPC scratch purging is built on; taking
the newest is what makes the degradations safe. A spuriously refreshed atime — a
backup, an indexer, a tree walk — makes an entry look newer and retains it, and a
`noatime` mount freezes atime at creation so the other two answer, which is no
worse than reading ctime alone. Neither deletes. Reporting the basis is what
distinguishes a `noatime` filesystem from genuinely cold artifacts.

The grace is `store_gc_grace` in days, default 30, overridden by `--grace`. It is
configurable because a store in a group root serves everyone using that install
and its owner knows its turnover — but it resolves per invocation, not per store,
which is why `--apply` requires `--dir` or `--layer`. Abandoned `.part` files are
collected on a separate fixed 24-hour window: that is crash recovery, and the only
question it answers is whether a producer could still be writing.

Nothing here consults reachability. Locks and project symlinks are not roots, a
manifest edge is not liveness — an artifact bakes its inputs in — and no build
dependency is ever installed, so every entry in the store is something someone
asked for.

The artifact cache is a disposable acceleration layer. Every candidate is
accepted only after a matching fingerprint and prior complete-key regeneration;
explicit validation bypasses it.
