# Overlay store

The store is immutable overflow below each images root. Flat overlays keep the
only bare-name slot; store entries are addressed by artifact name plus a
verified complete identity.

`Begin` first adopts an exact flat or stored artifact when one is already
available. Otherwise it selects a collision-safe store filename and acquires
the existing producer lock for that specific target. There is no store-wide or
store-directory lock. The producer writes directly to the prepared sibling
returned by the transaction, on the target filesystem.

`Commit` syncs the prepared file, regenerates its identity and equivalence keys,
checks its name and expected keys, rechecks the target, and publishes with the
same atomic sibling rename used by `BuildObject`. Store publication is
create-only: an exact target is adopted and a conflicting target is never
replaced. `Abort` removes only that producer's prepared output and releases only
that target's producer lock.

Producer locks coordinate creation because a missing target has no inode to
lock. Once an SQF exists, `exec` and `run` hold a shared inode `flock` while
reading it; operations that replace or remove an existing file require an
exclusive inode `flock`.

The artifact cache is a disposable acceleration layer. Every candidate is
accepted only after a matching fingerprint and prior complete-key regeneration;
explicit validation bypasses it.
