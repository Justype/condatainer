# Project lock

`cnt-lock/` is the tracked record of which exact artifact identity satisfies
each dependency a project declares. It holds the selection map plus the vendored
manifests and rebuild sources — everything Git should carry, and nothing
machine-local: no payload, no absolute path, no hostname, no restore result. A
checkout is a complete rebuild specification on a machine that has never run
CondaTainer.

```text
project/
  analysis.sh
  cnt-lock/
    lock.json
    artifacts/
      star--2.7.11b@a31f902c12ab/
        manifest.json
        recipe
```

## What is stored, and what is not

Only the mapping and non-derivable acquisition locations are written. Requests
come from rescanning the scripts; identity, type, platform, sources and
dependency edges come from the vendored manifests; closure membership comes from
following those edges. Nothing derivable is serialized, so a stale copy cannot
disagree with the truth.

`Marshal` publishes fixed field order, sorted keys, two-space indentation and one
final newline, so equal content produces equal bytes and a lock appears in a diff
only when it changed. `Unmarshal` rejects unknown fields rather than dropping
them: a lock written by a newer build may mean something this one would discard
on the next write.

Every path in a lock is hostile input — it names a directory a later step will
read — so `artifacts/<entry>` is required to be relative, clean, slash-separated
and non-escaping before anything touches the filesystem.

## Origins

`origins` maps an artifact path to the exact places it can be fetched from:
a repository coordinate and a **platform manifest digest**, never a mutable tag.
It is absent until something records a location, and it is keyed by artifact
rather than nested under a selection because a closure-only dependency can have
origins too.

Order is retry priority, so `AddOrigin` deduplicates and appends rather than
reordering. An origin is a location, never artifact metadata — nothing in it is
compared against a payload, which carries its own keys and is verified on
arrival.

**Origins are written by machines, never typed.** One is recorded only when
something has confirmed the artifact is there: a pull endpoint the recipe
collection declares advertises this exact identity, or `project push` has just
put it there. There is no flag for typing a coordinate by hand — a digest nobody
verified is a lock entry that fails on someone else's machine.

An origin cannot be recovered from the artifact either. It is
`repository@sha256:<digest over the pushed content>`, so embedding one would need
the digest before the bytes it covers exist; and a record written at pull time
would live in the per-user cache while images roots are shared, so it would exist
only for whoever ran the pull.

Nothing writes an origin yet, so restore rebuilds. The map is carried now because
the lock format is tracked in Git and adding it later would be a schema change.

## Scanning

A declaration counts wherever it is written. `catalog.ScanAnnotations` is the
shared tokenizer, so the scanner and the runtime read a script the same way and
moving a `#DEP:` cannot make the lock and `run` disagree about it. The cost is
that a heredoc writing another script contributes that script's declarations
too.

Discovery reads `.sh` and `.bash` files and extensionless files with a shell
shebang. It never follows a symlinked directory and never reads a symlinked
script: either can point outside the checkout, and a lock describes the checkout.
`.git/` and `cnt-lock/` are always skipped; anything else must be named in
`ScanOptions.ExcludeDirs`, because a silent skip is how a declaration goes
missing.

**An analysis script must name an exact version.** `#DEP: star/2.7.11b>=2.7.0`
is a build-recipe feature and is rejected here as a finding. A recipe declares a
range so a build can reuse a satisfying version already installed rather than
producing another — and for reference data the tool version often does not
matter, since `samtools faidx` writes the same `.fai` whichever recent samtools
ran. An analysis reports results from what it mounts, where "any version in this
range" is not a claim worth attaching to a result.

It is also what keeps one module from having two keys. The key is
`NameVersion() + Op + Min`, so `star/2.7.11b` and `star/2.7.11b>=2.7.0` would be
separate selections that could point at different artifacts, with nothing in the
lock to say which a given script meant.

A declaration is one of four kinds, and `Kind.Pinnable()` is the split that
matters:

| kind | pinnable | why |
|---|---|---|
| **name** | yes | resolved by identity across the images roots |
| **path** — a project-relative `.sqf` | yes | restore owns that path and writes the artifact there |
| **writable** — an `.img` | no | nothing to pin: it has no identity |
| **external** — a `.sqf` absolute or above the root | no | nowhere to put an answer: restore does not own that path |

The two unpinnable kinds fail for the same underlying reason — a selection
records where an artifact will be, and neither has such a place. An external
`.sqf` stays perfectly usable at run time; it just cannot be locked, because
locking it would record a promise restore could not keep without writing to
someone else's file.

## One anchor: the project root

Inside a project a relative path is relative to the **project root** — never to
the declaring script, never to the process working directory. The scanner keys a
`path:` request on the cleaned declaration text, restore materializes it at
`<root>/<key>`, and `run` resolves it the same way. So `../overlays/tool.sqf` in
a subdirectory script points outside the project and is external, which the
finding says outright, because that path looks script-relative and is not.

`project.WorkDir` is the other half. `container.ResolveOverlayPaths` resolves a
relative overlay against the *process* working directory, and a scheduler
chooses that: the submission directory under SLURM and LSF, `$HOME` under PBS —
and a script's own `#SBATCH --chdir` / `#PBS -d` / `#BSUB -cwd` is parsed and
re-emitted, so an author can relocate their own job. A project therefore states
the root when nothing was declared and refuses a declared directory that is not
the root.

Refusing rather than overriding: absolute overlay paths would survive any
directory, but the script's *own* relative paths — inputs, outputs — resolve
against the working directory. A job running elsewhere splits the two anchors.
A declared `--chdir` was written on purpose, so only its author can resolve the
conflict.

An unpinnable declaration must say so with `## unpinned` — see **The unpinned
marker** below.

## Where a restored artifact lands

What asked for an artifact decides where it goes, and the three are not
interchangeable. Depth in the dependency graph decides nothing: a directly
selected artifact is user-owned however deep it sits, and one reached only
through edges is scaffolding however shallow.

A **named** selection is addressed by identity. It is satisfied by a copy in any
readable images root and a new one lands wherever the store's destination rule
puts it — flat when the name is free, `store/` on conflict.

A **`path:`** selection is a project output. It must be materialized at exactly
the declared path, and a copy in an images root is *not* a substitute: adopting
one would report a successful restore while the declared path stays empty and
the script still fails at mount time. Such an artifact never enters the store.

A **build dependency** — a closure node no selection names, meaning a `#DEP:`
overlay and nothing to do with `#INPUT:` — is not installed at all.
It is produced in a restore-scoped directory under the stable writable tmp root
and removed when the restore ends, on failure and cancellation alike. Nothing
asked for it by name, and nothing needs it afterwards: an artifact bakes its
inputs in, so a dependency image is needed to rebuild it and never to use it.
One that is *already* installed is adopted in place and never copied, so only a
genuine miss is transient. `--keep-build-deps` installs newly produced ones
through the named rule instead.

The tree is one directory per restore, holding one directory per artifact:

```text
<writable tmp root>/cnt-restore-deps-<random>/
  samtools--1.21@c28134957e0b/artifact.sqf
  zlib--1.3@ff019a2bd4c1/artifact.sqf
```

The per-restore directory is created on first use and removed whole, so an
artifact two dependents need is produced once and mounted twice. Each staging
directory is named for the artifact — the encoded `name/version` and a short
identity, the same form `cnt-lock/artifacts/` uses — because this tree is what
someone reads when a build fails, and one entry per identity cannot collide
within a restore. The writable tmp root is the stable one, deliberately not the
`CNT_TMPDIR` fast root: several conda environments is more than node-local
scratch holds, and a job sweep must not remove it mid-restore. Nothing else sees
any of it — it is not an images root, so no scan finds it and `list` cannot show
it.

Two paths selecting one artifact are two files to produce, so they plan as two
steps. An artifact selected *both* by name and at a path has two destinations
and no way to choose; that is refused as a lock to fix.

The declared path is validated where the lock is parsed, not where it is used:
it must be relative, clean, slash-separated, inside the project, outside
`cnt-lock/`, and a `.sqf`. It is a write target, so an escaping one would put a
file outside the checkout entirely.

## What may satisfy a build dependency

A selection answers for its own keys: someone asked for it by name. A build dependency
only has to leave its dependent's equivalence key unchanged, and the role the
edge records is the scheme's statement about that.

| role of the edge | in the dependent's equivalence preimage | may be satisfied by |
|---|---|---|
| `data` | the dependency's equivalence key | anything carrying that equivalence |
| `app` | the dependency's `name/version` | that version, any build of it |
| `history` | nothing | any version of that name |

The exact recorded identity is preferred in every case, because it yields the
exact identity for the dependent too — the result that satisfies every other
project that locked it. Where two dependents disagree about a shared dependency,
the stricter role wins. Under `--match identity` nothing but the exact identity
is accepted, and not by policy: a script identity hashes dependency identities,
so a substitution produces a dependent the mode rejects anyway.

The recipe's `#DEP:` range is never consulted. It governed which version the
artifact was *first* built against, and the outcome of that decision is the
manifest edge; restore reproduces recorded outcomes rather than re-running
authoring decisions.

## Verification

`Verify` is strictly checkout-local: no installed overlay, no store, no catalog,
no build configuration, no network, and no payload. That is the CI contract — a
clean checkout on a machine that has never run CondaTainer either is a complete,
internally consistent rebuild specification or is not, and this says which.

Every entry is regenerated rather than trusted. The manifest's recorded keys are
recomputed from the vendored sources, and the directory name is checked against
the result: the name is addressing, never identity, so renaming a directory
cannot make a different artifact answer to it. Each dependency edge is followed
by name *and* complete identity to a child that must agree on both.

Reachability is the closure, not the selection set. A dependency directory is
reached through its parent's manifest edges, so anything computed from
selections alone would miss exactly the entries a rebuild needs. Problems are
collected rather than raised at the first failure, because someone fixing a lock
wants the whole list.

Reads are bounded and refuse symlinks. A lock is untrusted input: it names files
a later step will open, so nothing in it can escape the checkout or exhaust
memory.

## Publication

`Publish` marshals (which validates), writes a complete document to a dotted
sibling, renames it over the old lock, and only then prunes. A failure before
the rename leaves the previous lock authoritative and at worst leaves harmless
unreferenced directories; pruning never runs against a lock that was not
published.

`StageArtifact` builds an artifact directory as a dotted sibling on the same
filesystem and renames it into place, so a reader never sees a half-written
artifact. An existing entry is left alone — one `(name, identity)` has one set
of records, and rewriting them could only substitute different bytes claiming
the same key. A concurrent writer that wins the rename simply wins; its bytes
are ours by construction.

Dotted names are load-bearing: a write in progress is invisible to the artifact
listing, so an interrupted write is not a spurious validation problem.
`SweepStaging` clears what a crash left behind, and `Prune` leaves unreadable
entries alone — an entry that cannot be read cannot be shown to be unreachable,
and deleting on a read error would turn a corrupt file into data loss.

## The unpinned marker

An unpinnable dependency must say so:

```bash
#DEP: env.img  ## unpinned — scratch environment, rebuilt per machine
```

The marker is required and the reason is optional. Requiring prose would be
enforcement theater — nothing can tell a real reason from a placeholder — while
the step that matters, someone deliberately writing it on that line, has already
happened and `git blame` attributes it. A note that is not the marker is an
ordinary comment and claims nothing. The claim is about the artifact, so one
script marking a dependency settles it for the project.

Without the marker the scanner emits a **finding**, so `project validate` fails
rather than passing while the project mounts something unpinned; `project lock`
reports it as a warning and still publishes. Because the marker merges across
scripts, the finding is decided only once every script has been read.
