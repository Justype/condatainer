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
lock. Once an SQF exists, `exec` and `run` hold a shared inode `flock` while
reading it; operations that replace or remove an existing file require an
exclusive inode `flock`.

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
