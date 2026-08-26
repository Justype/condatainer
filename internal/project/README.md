# Project lock

`cnt-lock/` is the tracked record of which exact artifact identity satisfies
each dependency a project declares. It holds the pins plus the vendored
manifests and rebuild sources — everything Git should carry, and nothing
machine-local: no payload, no absolute path, no hostname, no restore result. A
checkout is a complete rebuild specification on a machine that has never run
CondaTainer.

```text
project/
  analysis.sh
  cnt-lock/
    lock.json
    provenance/
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
read — so `provenance/<entry>` is required to be relative, clean, slash-separated
and non-escaping before anything touches the filesystem.

One entry directory is one artifact's record set: `manifest.json` plus exactly
the sources that manifest names. It is the same directory an image carries at
`/.cnt/provenance/`, down to the name `capsule.EntryName` produces, and it is
read by the same code — `capsule.ReadRecord`, given this package's bounded,
symlink-refusing reader. Two readers would be free to drift into a lock that
verifies against an image that does not.

## Remotes

`remotes` maps an artifact path to the exact places it can be fetched from:
a repository coordinate and a **platform manifest digest**, never a mutable tag.
It is absent until something records a location, and it is keyed by artifact
rather than nested under a pin because a closure-only dependency can have
remotes too.

Order is retry priority, so `AddRemote` deduplicates and appends rather than
reordering. A remote is a location, never artifact metadata — nothing in it is
compared against a payload, which carries its own keys and is verified on
arrival.

`remotes` lives at the top of `lock.json`, keyed *by* an entry path, and never
inside the entry directory. That is the hashing boundary, not a layout
preference: an entry's bytes are the key preimage, so anything written in there
changes the identity it is addressed by. A registry coordinate is mutable
information about where a copy happens to sit — it must be able to change
without the artifact becoming a different artifact.

**Remotes are written by machines, never typed.** One is recorded only when
something has confirmed the artifact is there: a pull endpoint the recipe
collection declares advertises this exact identity, or `project push` has just
put it there. There is no flag for typing a coordinate by hand — a digest nobody
verified is a lock entry that fails on someone else's machine.

A remote cannot be recovered from the artifact either. It is
`repository@sha256:<digest over the pushed content>`, so embedding one would need
the digest before the bytes it covers exist; and a record written at pull time
would live in the per-user cache while images roots are shared, so it would exist
only for whoever ran the pull.

Exactly two writers exist, and a project normally uses both. `project lock` and
`project pin` record, best effort, any endpoint of the recipe collection that already
publishes the pinned artifact's **complete** identity — free, because the bytes
are already there. `project push` uploads the rest. Order is retry priority and
comes out right without a rule: a pin records its upstream remote before
any push runs, so the free location is tried first.

`Prune` is not what reconciles the map — `Publish` does, before it marshals, so
the bytes that get written never name an artifact the closure no longer reaches.
Pruning is filesystem-only and runs after the rename.

## Where a project publishes

`Lock.OCI` is the project's own destination: one repository coordinate and the
`audience` claim that decides what may be published there. It is tracked rather
than machine-local for the same reason a pin is — every collaborator has to
push to the same package, or the recorded fetch locations describe a set of
places *some* of the artifacts are.

Its coordinate grammar is duplicated here rather than imported from
`internal/registry`. That is what keeps this package free of the transport, which
is what makes `Verify` checkout-local by construction: a lock has to validate on
a machine with no network, no credentials and no registry, and a package that
cannot reach one cannot accidentally start.

`Lock.Source` is the project's code repository, derived once from the git origin
when the destination is set and recorded, never re-derived at push time — two
collaborators with different remote spellings would otherwise publish two
different source annotations for one project.

## Publishing the payload

`publish` decides and performs the upload. Planning splits in two so `--dry-run`
can be honest about cost:

- `Build` works from the checkout alone — the publish set, the tags, and the
  refusals. It covers every pin, and `--closure` adds the build
  dependencies reached through the manifest edges.
- `Refine` asks the network the two questions a checkout cannot answer: whether a
  collection already serves an artifact at this identity, and whether the
  destination already holds it. Both only ever *remove* work, so an unreachable
  registry leaves the plan as computed rather than failing a push nobody could
  complete offline. `--all` is a flag on this half, not on `Build`: it suppresses
  the upstream check rather than widening the set.

Push never builds. An artifact missing at its locked identity is a refusal naming
the restore that would produce it, because publishing something this checkout did
not already describe would put bytes in a registry that no lock vouches for.

Recording is **one lock transaction per artifact**, not one per push: load,
`AddRemote`, publish, next. An interrupted push therefore leaves every artifact
that did land recorded, and re-running skips them. Batching would trade that for
a single write nobody needs.

## Scanning

A declaration counts wherever it is written. `catalog.ScanAnnotations` is the
shared tokenizer, so the scanner and the runtime read a script the same way and
moving a `#DEP:` cannot make the lock and `run` disagree about it. The cost is
that a heredoc writing another script contributes that script's declarations
too.

Discovery reads `.sh` and `.bash` files and nothing else. `run` executes every
project script with `/bin/bash`, so no other shell's script could run here, and
a file with no extension is not read at all: sniffing a shebang to catch one was
a substring test that fired on any interpreter path containing `sh`, which on an
HPC filesystem means `/home/shared/…` and every user called `josh`.

`overlays/` is skipped: it is where a project keeps its own overlays and the
recipes that built them, and a recipe's `#DEP:` are its overlay's build
dependencies — already recorded in that overlay's provenance — not declarations
the project mounts. Declaring one from a project script is unaffected, since
`#DEP: overlays/tool.sqf` lives in the script rather than in `overlays/`.

A directory is the rule because only a stated convention answers the same every
time. The two alternatives both make a file's meaning depend on something
outside it: a sibling `.sqf` appears when `create -f` runs, so the scan would
change on a build, and deriving recipe-ness from which overlays other scripts
declare means deleting one script silently reclassifies another.

It never follows a symlinked directory and never reads a symlinked
script: either can point outside the checkout, and a lock describes the checkout.
`cnt-lock/`, `overlays/` and every dot-directory are always skipped; anything
else must be named in `ScanOptions.ExcludeDirs`.

The dot rule is the one inherited exclusion, and it is worth the silent skip it
costs. A dot-directory is tool state rather than project source — `.venv`,
`.tox`, `.snakemake`, a local conda env — and each ships shell scripts carrying
declarations that are not this project's. Since `project lock` now pins what it
finds, reading them turns a stray `#DEP:` into a failed lock rather than a note.

**An analysis script must name an exact version.** `#DEP: star/2.7.11b>=2.7.0`
is a build-recipe feature and is rejected here as a finding. A recipe declares a
range so a build can reuse a satisfying version already installed rather than
producing another — and for reference data the tool version often does not
matter, since `samtools faidx` writes the same `.fai` whichever recent samtools
ran. An analysis reports results from what it mounts, where "any version in this
range" is not a claim worth attaching to a result.

It is also what keeps one module from having two keys. The key is
`NameVersion() + Op + Min`, so `star/2.7.11b` and `star/2.7.11b>=2.7.0` would be
separate pins that could point at different artifacts, with nothing in the
lock to say which a given script meant.

A declaration is one of four kinds, and `Kind.Pinnable()` is the split that
matters:

| kind | pinnable | why |
|---|---|---|
| **name** | yes | resolved by identity across the images roots |
| **path** — a project-relative `.sqf` | yes | restore owns that path and writes the artifact there |
| **writable** — an `.img` | no | nothing to pin: it has no identity |
| **external** — a `.sqf` absolute or above the root | no | nowhere to put an answer: restore does not own that path |

The two unpinnable kinds fail for the same underlying reason — a pin
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
pinned artifact is user-owned however deep it sits, and one reached only
through edges is scaffolding however shallow.

A **named** pin is addressed by identity. It is satisfied by a copy in any
readable images root and a new one lands wherever the store's destination rule
puts it — flat when the name is free, `store/` on conflict.

A **`path:`** pin is a project output. It must be materialized at exactly
the declared path, and a copy in an images root is *not* a substitute: adopting
one would report a successful restore while the declared path stays empty and
the script still fails at mount time. Such an artifact never enters the store.

A **build dependency** — a closure node no pin names, meaning a `#DEP:`
overlay and nothing to do with `#INPUT:` — is not installed at all.
It is produced in a restore-scoped directory under the stable writable tmp root
and removed when the restore ends, on failure and cancellation alike. Nothing
asked for it by name, and nothing needs it afterwards: an artifact bakes its
inputs in, so a dependency image is needed to rebuild it and never to use it.
One that is *already* installed is adopted in place and never copied, so only a
genuine miss is transient. `--keep-build-deps` installs newly produced ones
through the named rule instead — with one difference: where a pin in the
same plan carries the same name at a different identity, the pin keeps the
bare name and the dependency goes to `store/`. The flat name answers `exec -o`,
`list`, and every other checkout, so it belongs to what someone asked for by name
rather than to what a rebuild happened to need. `yieldSharedNames` decides it from
the plan, because leaving it to the store means the first installer wins and
dependencies are built first — exactly backwards.

The tree is one directory per restore, holding one directory per artifact:

```text
<writable tmp root>/cnt-restore-deps-<random>/
  samtools--1.21@c28134957e0b/artifact.sqf
  zlib--1.3@ff019a2bd4c1/artifact.sqf
```

The per-restore directory is created on first use and removed whole, so an
artifact two dependents need is produced once and mounted twice. Each staging
directory is named for the artifact — the encoded `name/version` and a short
identity, the same form `cnt-lock/provenance/` uses — because this tree is what
someone reads when a build fails, and one entry per identity cannot collide
within a restore. The writable tmp root is the stable one, deliberately not the
`CNT_TMPDIR` fast root: several conda environments is more than node-local
scratch holds, and a job sweep must not remove it mid-restore. Nothing else sees
any of it — it is not an images root, so no scan finds it and `list` cannot show
it.

A destination step whose path already holds something that does not answer the
lock records it as `Step.Replaces`, read from the file itself during planning.
`placeAt` renames over whatever regular `.sqf` it finds, so without that a plan
would preview an overwrite as an ordinary build — indistinguishable from an
empty path. Replacing is right (a restore exists to make the checkout match the
lock, and refusing would leave a drifted project repairable only by hand), but
silently replacing is not.

Two paths pinning one artifact are two files to produce, so they plan as two
steps. An artifact pinned *both* by name and at a path has two destinations
and no way to choose; that is refused as a lock to fix.

The declared path is validated where the lock is parsed, not where it is used:
it must be relative, clean, slash-separated, inside the project, outside
`cnt-lock/`, and a `.sqf`. It is a write target, so an escaping one would put a
file outside the checkout entirely.

## What may satisfy a build dependency

A pin answers for its own keys: someone asked for it by name. A build dependency
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

## Scheduler submission

A rebuild whose vendored recipe carries `#SBATCH`/`#PBS`/`#BSUB` directives is
handed to the scheduler instead of run in place, and the job re-enters restore as
`project restore --project <root> --only <artifact>`. Nothing about the build is
serialized into the job script: the lock is the specification and the node reads
the same checkout, so the job needs only to be told which artifact it is for.

`partition` decides the split before the first step runs, because it cannot be
decided step by step — a build dependency comes before its dependent in the order,
and whether it runs here or inside that dependent's job is only knowable once the
dependent's fate is known. Three rules produce it:

- **A build dependency is never its own job.** It lives in a restore-scoped
  directory, and a directory cannot cross a job boundary; the job that needs it
  produces it inline. `--keep-build-deps` installs it instead, which makes it an
  ordinary step with its own job and its own directives — the escape hatch when a
  dependency needs more of a machine than its dependent does.
- **A step needs the scheduler when its own recipe declares directives, or when a
  build dependency it produces inline declares them.** Otherwise a heavy
  dependency would drag its light dependent onto the login node. Directives carry
  up a whole inline chain, not one level.
- **A step whose dependency was submitted is submitted too**, waiting on it with
  an `afterok` edge. Its input does not exist yet, so it cannot run here whatever
  it declares.

Each declared project path is one output and therefore one job, even when two
paths pin the same artifact.

### Reservation and resume

The submitting process reserves where the job will publish and leaves an ordinary
producer lock there carrying the job ID — through `store.Transaction.Detach` for a
store-addressed artifact, and directly on the file for a project path. There is no
second mechanism and no project-local job-state file.

That one lock answers every question resume has to ask, because
`producer.AcquireLocal` already knows how to read it:

| lock at the target | meaning | what restore does |
|---|---|---|
| absent | never submitted | submit |
| held, job alive | submitted by an earlier restore | report the job, submit nothing |
| held, job gone | the job died without publishing | clear it and the orphaned output, resubmit |
| held, job is *this* one | the worker reached its own reservation | adopt and publish |

Re-running the restore is the resume operation. The command exits with the
jobs-submitted code, having made nothing available yet; a later run adopts what
the finished jobs published.

The job's working directory is the project root or the submission is refused
(`project.WorkDir`) — a project's relative paths resolve against the root, and a
job that runs anywhere else splits CondaTainer's overlay resolution from the
script's own paths.

## Acting in a project

`run`, `check` and `exec` act *in* whatever project the caller is standing in.
There is no flag for it either way: `--project DIR` belongs to the commands that
act *on* a project, and standing somewhere else is the opt-out.

`exec -o` and `e -o` share one hook, `cmd.projectOverlays`. An argument classifies
through `lock.ParseDeclaration` — the same grammar a `#DEP:` uses — so a name
typed on the command line and the same text in a script cannot mean different
things. That means a version constraint is refused here too, and a project path
answers to the lock rather than being mounted on sight: inside a project
`overlays/tool.sqf` is a restore output, and `LookupAt` verifies the file there
against the locked keys. A writable `.img` and an external `.sqf` stay unpinnable
and are mounted as written.

The hook sits *above* `container.ResolveOverlayPaths` rather than inside it, even
though that is the one place a name becomes a path. Five of its callers resolve a
build's own `#DEP:` or a base image, and none of them may pick up the lock of
whatever directory the user happened to be standing in.

Like `run`, it resolves and never acquires: an absent artifact is an error naming
`project restore`, never a fetch, a build, or a fallback to whatever currently
answers to the name. That fallback is the failure a lock exists to prevent.

## Verification

`Verify` is strictly checkout-local: no installed overlay, no store, no catalog,
no build configuration, no network, and no payload. A fresh clone that has
restored nothing either is a complete, internally consistent rebuild
specification or is not, and this says which. Whether an artifact is actually
there to run is `Resolve`'s question, not this one.

Every entry is regenerated rather than trusted. The manifest's recorded keys are
recomputed from the vendored sources, and the directory name is checked against
the result: the name is addressing, never identity, so renaming a directory
cannot make a different artifact answer to it. Each dependency edge is followed
by name *and* complete identity to a child that must agree on both.

Reachability is the closure, not the pin set. A dependency directory is
reached through its parent's manifest edges, so anything computed from
pins alone would miss exactly the entries a rebuild needs. Problems are
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
