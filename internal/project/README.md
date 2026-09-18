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
    .gitignore
    .helper-history/
      rstudio-server-root-1737654321000000000.json
```

`cnt-lock/.helper-history/` is the one deliberate exception to "everything
under `cnt-lock/` is tracked": it is gitignored, not committed. `cnt-lock/.gitignore`
(`/.*/`, written by `lock.EnsureIgnore` whenever the lock is published or a
history entry is recorded, and left alone once it exists) ignores every dot
directory there, so machine-local state needs no edit to the project's own
`.gitignore`. It holds one
JSON file per distinct overlay combination a helper has been run with in
this project — shared across every user via the same group-writable
`cnt-lock/`, so a teammate running `rstudio-server` from this project sees
what already worked instead of rediscovering it — but the data is
machine-local recency (a file's own mtime says which combination is
newest), not part of what a checkout needs to rebuild. A `git checkout`
restamping every file's mtime would silently break that, which is exactly
what tracking it would risk. See `internal/helper/README.md`,
"Shared, cross-user helper-setup history" for the format and the write path.

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

## Manual pins

`Pins` is otherwise a projection of the scan: `Reconcile` deletes every key the
current scan does not produce, which is what makes the map safe to rebuild from
nothing. A **manual** pin opts out of that one rule and changes nothing else.

It exists because a project depends on artifacts no `#DEP:` names:

- a helper declares `#REQUIRED_OVERLAYS: r{POSIT_R} rstudio-server
  build-essential`, which the project scanner never reads — it reads `#DEP:` in
  project scripts, and a helper's overlay list is a different declaration
  entirely. `Standing.ResolveNames` (`standing.go`) is how those names still
  reach the pins once launched with a project as their working directory:
  `internal/helper.CheckRequiredOverlays` calls `project.StandingAt` first,
  falling back to its ordinary on-disk check and auto-install only when there
  is no project to resolve against;
- a **frozen environment** is named by nobody. It is what downstream analysis
  runs *in* rather than something a script consumes, so nothing selects one per
  script and no `#DEP:` ever mentions it.

Everything else about such a pin is ordinary. It vendors its provenance, records
remotes in `Lock.Remotes`, restores and publishes through the same steps — a
frozen environment included, since its identity is hashed from the packed payload
and lives in its manifest like every other artifact's.

`project pin` sets the flag by asking whether anything declares the request.
Derived once, at pin time, and then **recorded**: "not currently declared" is
exactly what the sweep tests, so deriving it on every run would turn the sweep
into a no-op. A scan that cannot run answers "declared", which leaves the pin
sweepable rather than silently permanent — the cautious answer is the one a later
`project lock` can still correct.

Re-pinning never clears the flag. What changed is which artifact answers, not why
the pin is in the lock.

`project unpin` is how one goes, and the only way: the sweep cannot reach it by
construction. It deletes the key and republishes, which is the whole
implementation — `Publish` already prunes the provenance directories nothing
reaches and the remotes addressing them. A pin the scan produces is refused there
rather than removed, because the next `Reconcile` re-pins the same declaration,
possibly at a different identity.

No scan-shaped view can show a manual pin, since no scan produces one. `project
list` is the lock-shaped view: it reads `Pins` and the records beside them and
marks which are manual, so the pin key `unpin` takes is discoverable. A pin whose
artifact does not verify is listed as unreadable rather than dropped from the
listing — the key is what addresses it, and hiding it would hide the thing to
fix.

### Used but not pinned, and who uses a manual pin

`UnpinnedHelperOverlays` and `ManualPinUsage` (`lock_cmd.go`) both read
`usageIndex`, a reverse index built from the project's shared helper-setup
history (`internal/helper/README.md`): which pin key each recorded overlay
combination would address, mapped to the helper names that used it.
`UnpinnedHelperOverlays` is that index's keys with no matching pin — what
answers "why does a helper's `#REQUIRED_OVERLAYS:` refuse instead of
resolving ambiently" concretely: the strict resolution stays exactly as
strict, and this only adds discoverability, surfaced in `project
status`/`project lock`'s report under its own heading, never merged into the
`#DEP:`-derived unpinned list and never written as a pin. `ManualPinUsage` is
the same index read the other way — for each manual pin, which helpers are
recorded using it, the manual-pin equivalent of a `#DEP:` pin's own
`Request.Scripts`. An empty result is reported as a fact ("no recorded helper
usage"), never a suggestion to unpin: nothing about a helper that has not run
since this shipped, or one that runs rarely, is distinguishable from "no
longer needed" by usage alone, and this package's principle — usage is never
evidence of importance — applies to its absence too.

## The project's root

`lock.BaseKey` (`"base:"`) is a reserved `Pins` key holding the artifact
`restore` treats as the container root, instead of falling back to whatever
`config default_distro` says on the machine running it. It is an ordinary
`PinEntry` — same shape, same `Manual` meaning as every other pin — addressed
by a fixed key rather than by name, because there is no manifest type left
that would say "this one is a root": `catalog.TypeOS` covers both an ordinary
`os` overlay and whatever plays root. The key's trailing colon is what a
rendered `#DEP:` request can never produce, the same guarantee `PathPrefix`
gives a path pin.

Unlike every other pin, it is written unconditionally: `project lock` always
derives it, whether or not anything in the project actually reaches outside
`/cnt_env`. There is no closure to walk deciding whether one is needed, and no
cost recording one when nothing does — see `plan/base-provenance.md` for why
that changed once Apptainer and Micromamba stopped being something a root
provided. Derivation has exactly one source, `config.ResolvedDefaultDistro()`;
`project select-distro <distro>` overrides it with `Manual` set, and
`--auto` clears the override back to derived. `Reconcile`'s sweep leaves the
key alone unconditionally — it is never something the `#DEP:` scan produces,
so treating an absent scan hit as staleness would delete it every run.

Restore gives it two guarantees an ordinary pin does not need:

- **It runs first.** `reorderBaseFirst` (`restore/plan.go`) moves its step to
  the front of the plan. Nothing in the manifest graph points at it — an `os`
  artifact may not carry `#DEP:` — so topological order alone would leave it
  wherever it falls; every conda or script rebuild needs its resolved path
  before it can build.
- **`--only` cannot drop it.** A submitted job re-enters `Compute` restricted
  to one artifact and its dependency closure, and the base is in nobody's
  closure. `restrict` keeps it anyway, so a job rebuilding one artifact on a
  compute node resolves its root from the lock the same way an unrestricted
  restore does, rather than falling back to that node's own configuration.

Once its step has run, `Run` reads its resolved path out of the same
`available` map every other dependency edge uses, and hands it to
`build.LockedSpec.Base` for whichever conda or script rebuild follows —
`internal/build`'s `resolveBase` early-returns when `Spec.Base` is already
set, so a locked rebuild never reaches `config.GetBaseImage()`. A `.def`
rebuild ignores it: it bootstraps its own root and never reads `Spec.Base`.

`cmd.ensureRootBaseImage` gives `exec`, `e` and `run` the same guarantee for an
ordinary session, not just a restore: standing in a project, if nothing the
caller requested is itself root-eligible, `cmd.projectBaseImage` calls
`Standing.Base`, which resolves the reserved base pin the same way
`Standing.ResolveComplete` resolves any other pin — strictly, so an unresolved
root is a refusal naming `project restore`, never a silent fall back to this
machine's `default_distro`. All three commands reach this through the one
function, so none of them can disagree about which root a project means.

`cmd.projectDefaultDistro` extends the same substitution to every bare-name
shortcut: `avail`, `list`, `remove`, `info`, `overlay`, shell completion, and
`create`'s own bare-name expansion all called `config.ResolvedDefaultDistro()`
directly, which meant `build-essential` could mean a different artifact
depending on whose machine typed it — this plan's failure arriving through the
name table rather than through the mount. Standing in a project,
`cmd.projectSelectedDistro` calls `Standing.SelectedDistro`, which reads the
distro straight out of the base pin's vendored manifest name —
`select-distro`/`DeriveBase` only ever compose `<distro>/base`, so splitting on
the first `/` is exact — without resolving a local copy, so it answers the
same in a fresh clone that has restored nothing. Unlike `Base` this never
refuses: a lookup that fails for any reason falls back to the configured
`default_distro`, which is the right default for a completion or display path.

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

Push gates on `project validate`, not on `lock.Verify` alone. Verify reads
`cnt-lock/` and answers whether the lock is internally consistent; only a scan of
the scripts sees a declaration nothing can pin or a declaration with no pin, and
neither ever becomes a pin for Verify to object to. A project with a writable
`.img` in a `#DEP:` therefore had a perfectly valid lock and published clean,
while failing validate — a published project that no other checkout can restore,
which is the one thing publishing is for. `cmd.projectProblems` is the single
definition both commands read, so the two cannot drift apart again.

A dry run reports those problems and still builds the plan: its job is to say
what would happen, and refusing to answer is less useful than answering with the
gap named.

`publish` decides and performs the upload. Planning splits in two so `--dry-run`
can be honest about cost:

- `Build` works from the checkout alone — the publish set, the tags, and the
  refusals. It covers every pin, and `--closure` adds the build
  dependencies reached through the manifest edges. Two pins sharing a manifest
  name lose the plain tag here rather than in `Refine`, because the collision
  check runs immediately after: a shared plain tag left in place would refuse
  the whole push offline, and two frozen environments — both named `env` — make
  that ordinary rather than exotic.
- `Refine` asks the network the two questions a checkout cannot answer: whether a
  collection already serves an artifact at this identity, and whether the
  destination already holds it. Both only ever *remove* work, so an unreachable
  registry leaves the plan as computed rather than failing a push nobody could
  complete offline. `--all` is a flag on this half, not on `Build`: it suppresses
  the upstream check rather than widening the set.

Push never builds. An artifact missing at its locked identity is a refusal naming
the restore that would produce it, because publishing something this checkout did
not already describe would put bytes in a registry that no lock vouches for.

**A frozen environment inverts the usual fallback.** Every other artifact can be
rebuilt from what the lock vendors, so a registry is an optimisation; a snapshot
vendors nothing and was captured rather than built, so a registry copy is the
only thing that can produce it. `classify` therefore refuses it at plan time as
`ActionUnavailable` rather than planning a build that must fail during
acquisition — the plan is where a restore is allowed to be honest about what it
cannot do. Push is what closes the gap, which is why the refusal names it.

`SkipPrebuilt` is deliberately not consulted for one: `--no-prebuilt` chooses
building over fetching, and where building is impossible there is no preference
left to express. Honouring it would refuse a restore that a recorded remote can
satisfy, and refuse it with a message naming the `project push` that had already
been done.

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

**A build recipe is skipped, and what identifies one is `$CNT_PREFIX`.** A
recipe's `#DEP:` are its artifact's build dependencies — already recorded in that
artifact's provenance — not declarations the project mounts, so reading them
would make the project pin what it never mounts. Declaring the built overlay is
unaffected: `#DEP: overlays/tool.sqf` lives in an analysis script.

`$CNT_PREFIX` is where a build writes its payload, so only a recipe expands it
and every recipe must. That makes it a property of the file, which is what the
rule needs: a directory convention classified `overlays/run-qc.sh` as a recipe
for sitting in the wrong folder, and left a recipe elsewhere in the tree read as
an analysis script. Comments are stripped first, so prose naming the variable is
not a declaration of intent, and both `$CNT_PREFIX` and `${CNT_PREFIX}` count
because both are ordinary shell.

The two rejected alternatives both make a file's meaning depend on something
outside it: a sibling `.sqf` appears when `create -f` runs, so the scan would
change on a build, and deriving recipe-ness from which overlays other scripts
declare means deleting one script silently reclassifies another. A skipped
recipe is not counted among the project's scripts at all.

`ScanScript` applies the same rule, because a run and a lock must not disagree
about what a script declares.

It never follows a symlinked directory and never reads a symlinked
script: either can point outside the checkout, and a lock describes the checkout.
`cnt-lock/` and every dot-directory are always skipped; anything else must be
named in `ScanOptions.ExcludeDirs`.

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
the process working directory. The scanner keys a `path:` request on the
cleaned declaration text, restore materializes it at `<root>/<key>`, and `run`
resolves it the same way.

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

An unpinnable declaration is a finding — see **Unpinnable declarations** below.

### A declaring script's own directory is a fallback address, never the anchor

The rule above governs one thing: what a declared path *is* — pinnable,
external, or in bounds. What it might *mean when a real file has to be found
for it* is a separate question, and `Request.PathCandidates()`
(`lock/scan.go`) is where that lives. For a `path:` request, candidate 0 is
always the literal, root-relative key — a script written in the recommended
layout (an analysis script at the root, `overlays/` beside it) never sees
anything else. Candidates 1+ are the same declared suffix taken relative to
each script that declared it instead, read off `Request.Scripts`. A script
paired one-for-one with its own overlay — `steps1/run.sh` declaring
`#DEP: xxx.sqf` for a file that lives at `steps1/xxx.sqf` — resolves without
its author having to spell out `steps1/` themselves.

The two places this list gets consumed check two different things, and
neither may check the other's:

- **`PinAll` checks the filesystem**, trying each candidate until one answers
  to a real `.sqf` it can hash — sound only here, because a file has to exist
  already to be pinned at all.
- **`Reconcile` and `Resolve` check `l.Pins` membership only**, via
  `lock.MatchPin`, never the filesystem — a fresh checkout has to resolve
  correctly before anything has been restored, and neither may invent an
  answer from files that are not there yet.

For a `name:` request, `MatchPin` has one more fallback past its
`PathCandidates`. With no pin under its own key, it matches an existing pin
whose name agrees and whose version the request's `catalog.Dep` admits — a bare
name admits any version, a partial one any version it is a dot-component prefix
of — taking the newest when several qualify. This keeps a project consistent
when two scripts name the same dependency at different precision: whichever
one a prior `project lock` pinned is what the other means too, not an
independent re-resolve that could land on a different version. It is
filesystem-free. A name nothing in `l.Pins` matches is answered live — see
**Live resolution of an unpinned name** in *Acting in a project*.

Whichever candidate a real file answered to becomes the **stored** key —
`path:steps1/xxx.sqf`, never the literal text as declared if that is not the
one that matched. A `path:` pin is always root-relative once recorded; the
ambiguity exists only in how a bare declaration gets interpreted once, at the
moment something real answers for it.

One more piece follows from this: a `../`-leading declaration is not
classified purely from its own text any more, for a *scanned* (script-owned)
declaration. `../overlays/tool.sqf` written in `steps1/run.sh` escapes when
read as root-relative, but relative to `steps1/` it means
`overlays/tool.sqf` — safely inside the project. `finalize`'s
`reclassifyEscaped` (`lock/scan.go`) tries each declaring script's directory
before giving up and reporting the declaration unpinnable, and only promotes
it to `KindPath` when the resolved key is not already claimed by a different
declaration — two different declaration texts must never end up sharing one
identity by accident. `exec -o` and manual `project pin` are untouched by any
of this: with no declaring script there is exactly one candidate, root-relative,
forever — the same grammar `ParseDeclaration` has always given a name typed on
the command line.

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

`run`, `check`, `exec` and a helper's `#REQUIRED_OVERLAYS:` all act *in*
whatever project the caller is standing in. There is no flag required either
way, by default: standing somewhere else is the ordinary opt-out.

`Standing` (`standing.go`) is what each of them stands on: `StandingAt(cwd)`
finds the project rooted at or *above* `cwd` — walking up through ancestors —
or answers `(nil, nil)` when no ancestor has one, so a caller's own ordinary
resolution takes over unchanged. Every project-aware entry point in the tree
is a method on it — `ResolveComplete`/`ResolveNames` for overlays,
`Base`/`SelectedDistro` for the root — so `cmd`'s exec/run hooks and
`internal/helper.CheckRequiredOverlays` share one implementation of "find the
root, load the lock, resolve, refuse with the same message" instead of three
packages repeating it.

`StandingAt` walking upward means one directory deeper than the root is still
standing in the project, not an ambient fallback — `exec scripts/align.sh`
behaves the same whether it's typed from the root or from `scripts/`. Only
`Standing` does this: `lock.RootFor`'s cwd-fallback branch, used by the
`project` subcommand family when `--project` is omitted, keeps calling the
strict, non-walking `lock.RootAt` — a write command reaching an ancestor
project from a subdirectory that looks empty is a hazard an ambient read
never is, and `project lock`'s own "no lock here yet, so make one" behavior
depends on "no lock here" meaning literally the named directory.

`exec`, `e`, `run` and `check` each carry two more flags on top of the ambient
default, registered together by `RegisterProjectFlags`
(`cmd/project_flags.go`): `--project DIR` relocates the whole invocation —
`os.Chdir(DIR)` once, before anything else runs, so `StandingAt` and every
other relative argument (a script path, an `-o` local file, a `--bind` path)
resolve against `DIR` exactly as if the caller had `cd`'d there themselves,
never as a second anchor kept alongside the real one. `--no-project` instead
disables lookup outright, behaving as if `StandingAt` had found nothing,
regardless of what standing there would otherwise resolve — the deliberate
opt-out that walking upward would otherwise remove. The two are refused
together. `internal/helper.RunOptions.NoProject` gives helper launches (CLI
and dashboard alike, since both share `PlanRun`) the same opt-out; helper
already had the relocation flag's equivalent, since a launch's `cwd` is
always given explicitly rather than inherited from an ambient process
directory.

`exec -o` and `e -o` share one hook, `cmd.projectOverlays`, and `run`'s script
scan uses `cmd.projectRunContext` — both call `Standing.ResolveComplete` after
building their own `[]lock.Request`. A request classifies through
`lock.ParseDeclaration` — the same grammar a `#DEP:` uses — so a name typed on
the command line and the same text in a script cannot mean different things.
That means a version constraint is refused here too, and a project path
answers to the lock rather than being mounted on sight: inside a project
`overlays/tool.sqf` is a restore output, and `LookupAt` verifies the file there
against the locked keys. A writable `.img` and an external `.sqf` stay unpinnable
and are mounted as written. Both also set `ResolveOptions.LiveResolve` and
`Distro: standing.SelectedDistro()` — see **Live resolution of an unpinned name**
below for what that changes. `check` (`cmd.projectCheck`) sets the same two
fields, calling `project.Resolve` directly, so it answers "can this script run"
with exactly what `run` would do.

`Standing` sits *above* `container.ResolveOverlayPaths` rather than inside it,
even though that is the one place a name becomes a path. Five of its callers
resolve a build's own `#DEP:` or a base image, and none of them may pick up the
lock of whatever directory the user happened to be standing in.

Like `run`, every method resolves and never acquires: an absent artifact is an
error naming `project restore`, never a fetch or a build, in every case and
regardless of `LiveResolve`. Whether an *unpinned* name refuses is a separate
question, answered per caller — see below. `Standing.ResolveComplete`'s
refusal is the one message every caller shares when it does refuse, so a `-o`
typo, an unpinned helper overlay and an unresolved root all fail the same
recognizable way.

### Live resolution of an unpinned name

A `path:` declaration is a claim to be pinned, so one with no matching pin is
always that refusal. A `name:` request — bare, partial or exact — is held to an
exact identity when it is *pinned*, not when it runs: `ResolveOptions.LiveResolve`
lets `Resolve` answer one no pin matches with `catalog.SolveInstalled`, exactly
as an unpinned name resolves outside a project: an installed version first, the
catalog's own candidates only once nothing is installed, and never a build.
Nothing installed answering the name is the same "declared but not pinned"
error, and a catalog-only hit with nothing on disk is still unresolved.

`exec -o`, `run`'s script scan and `check` set `LiveResolve: true`, matching
the fact that `-o` and `#DEP:` already classify through the identical grammar
and are meant to behave alike whether or not a project happens to be standing.
A helper's `#REQUIRED_OVERLAYS:` (`Standing.ResolveNames`) and the project's
own root (`Standing.Base`) leave it unset: both are fixed requirements a
project locks rather than lets float, and their own callers already document
why — see **Manual pins** and **The project's root**. A resolved mount that
came from `LiveResolve` carries `Mount.Live` rather than a recorded
`Identity`, since nothing was pinned to compare against; every caller that
sets `LiveResolve` prints a note naming what it resolved to, so the outcome is
visible without needing to be reconstructed.

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

## Unpinnable declarations

A declaration nothing can pin is always a **finding**. `project validate` fails
on it and `project lock` reports it as a warning and still publishes what it
could. Nothing written in a script silences one.

There is no opt-out because an unpinnable dependency is exactly the case the
lock exists to prevent: a project that mounts something no one else can obtain
reproduces its results nowhere, and a lock that records the rest looks complete
while the environment the analysis actually ran in is missing from it. A marker
that quieted the warning would make the hole a formality rather than a problem,
and the artifact it hides is usually the environment itself.

Each kind names what closes it, since the two are unpinnable for different
reasons and only one remedy applies:

- **A writable `.img`** has no identity to hash — its content changes under any
  reader. `overlay freeze` packs it into a `.sqf` with a `snapshot-env-v1`
  identity, and a project-relative `.sqf` is `KindPath`, which pins. This is the
  case freeze exists for.
- **An external `.sqf`** is already immutable; the problem is that restore does
  not own that path and must never write there. Copying it under the project
  makes it a path restore can own.

Because a declaration merges across the scripts that make it, one finding is
reported per dependency rather than per line, which is why `finalize` decides
this after every script has been read rather than as each line is parsed.
