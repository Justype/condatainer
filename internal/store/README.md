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

Producer locks coordinate creation because a missing target has no inode to
lock. Once an SQF exists, `exec` and `run` hold a shared inode `flock` while
reading it; operations that replace or remove an existing file require an
exclusive inode `flock`.

The artifact cache is a disposable acceleration layer. Every candidate is
accepted only after a matching fingerprint and prior complete-key regeneration;
explicit validation bypasses it.
