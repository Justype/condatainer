# catalog

Resolves a module name to a recipe, and a set of names to a build order.

It is a **getter and a dependency solver**. It never builds, never chooses, and
never talks to conda — every decision its answers feed belongs to the caller. 
Shared between two tools, so anything opinionated about *what to do* with a 
result stays out.

A collection is a directory or a URL holding `recipes/`, `index/`, `helpers/` and
a `source.json` descriptor. Configured sources are ordered and **first match
wins**, like `PATH`.

## The boundary

Conda is deliberately absent. A name no source provides is not an error — it is
the normal path for a conda package, and the node comes back with a nil entry and
its constraint intact for the caller's fallback. Nothing here pre-solves conda:
micromamba solves against the whole environment, and a second opinion computed
beforehand could only agree or be wrong.

Likewise absent: the scheduler. Recipes read normalized `$NCPUS`/`$MEM`, which the
scheduler package produces. The edge runs from each tool to both packages, never
between them.

## `SolveName`: a raw name to a concrete module

`Lookup` and `Resolve` require an exact index key or a fully filled template
target; `SolveName` is the same boundary (conda absent, not-found an outcome)
extended to a name that is bare, carries a partial version, or is missing its
distro prefix — what a person actually types, or what `#REQUIRED_OVERLAYS:`
names.

Two attempts, in order: `raw` exactly as given, then, only once that comes up
empty and `raw` has at most one slash, the same attempt under `distro/raw`.
Trying the literal name first is not a guess; prepending a distro is, so it
only happens after the unprefixed form has already failed.

Each attempt is itself three steps, tried in order and stopping at the first
hit:

1. **A complete module, no splitting at all.** `Lookup(raw)` directly — an
   exact index key, or a template with its placeholder already filled. A hit
   on a *bare* template (the key exists, nothing filled its axis) does not
   count here; that names a family, not a module, and falls through to the
   next step. This alone resolves anything already fully specified, of any
   type — `star/2.7.11b`, `ubuntu24/xfce4`, `grch38/genome/ucsc` — with no
   name/version guessing at all.
2. **App pool**, only tried when `raw` has at most one slash: split at the
   one possible `/` (`name`, or `name/version`). `have` is checked first,
   and if it already answers, the catalog is never consulted at all —
   installed beats a catalog lookup entirely, not just beats a newer
   version, the same "check what's already installed first, no network
   call" rule the load side follows. Only once nothing is installed does the
   catalog decide: a flat sibling one segment below `name`, or a `#PH:` axis
   on a template named exactly `name`. Always safe to pick newest here,
   because an app's candidates, flat or templated, are the same tool at
   different points in time by construction, never a different tool sharing
   a name prefix.
3. **OS pool**, only tried when `raw` has at most two slashes: the first
   segment is `distro`, the rest is `name` or `name/version`. Same
   have-first order as the app pool; only once nothing is installed does the
   entry at `distro/name` have to exist and declare the `#PH:` axis a
   version query is checked against — never a scan of the distro's other
   children. The catalog side of this is a direct, exact-key lookup, not a
   scan, so it costs the same whether the catalog holds a dozen entries or a
   hundred thousand.

An installed overlay resolves through either pool even when the catalog has
nothing backing it — an empty source, or a recipe that has since moved or
been dropped. `Resolved.Entry` is nil in that case; there is simply nothing
to attach.

**Data gets no third step of its own — and that omission is what "no autofill
for data" means in practice.** `grch38/genome` (bare) fails step 1 (not a
literal key), fails step 2 (`grch38` matches no `TypeApp` candidate), and
fails step 3 (`grch38` is not a distro, so no entry sits at `grch38/genome`)
— it simply runs out of steps, with no scan of `grch38/genome`'s own flat
children (`ucsc`, `ensembl`) ever happening. That scan would be unsound
regardless: those are alternative sources, not versions of one another, and
comparing them via `CompareVersions` would be comparing two arbitrary
strings. `grch38/genome/ucsc` typed in full still resolves normally, via step
1 — only the *choice* among data siblings is refused, never a name that
already names one exactly.

**A bare distro name has no versions of its own**, for the same reason:
`"ubuntu24"` alone is one segment, so step 2 finds no app named that and step
3 has no `/` to split into distro+name — it is not "matched too many
candidates and picked wrong," there is no candidate step it can even reach.
`condatainer create ubuntu24` means nothing on its own and reports not found
rather than guessing. A misparsed multi-segment bare name recovers the same
way any other bare name does: `"ubuntu24/rstudio-server"` fails step 1 (bare
template) and step 2 (not app-shaped, more than one slash), then step 3 finds
the one entry at `ubuntu24/rstudio-server` and its `#PH:` axis directly — no
separate "is this actually a distro" check needed, because step 3 never scans
for what a distro's children *are*, it only ever looks up the one address it
was given.

`have` and `distro` are always parameters, never derived from ambient state.
That is what lets `SolveName` mean the same thing regardless of which package
calls it — `cmd` (a project's selected root, or config `default_distro`),
`internal/helper`, `internal/project/lock` — without any of them needing to
import back into `catalog` to supply it, and without the name resolving to a
different artifact depending on whose machine typed it.

## The header boundary

`ScanAnnotations` is the one place a `#KEY: value ## note` line is tokenized.
Recipes, user scripts and project scanning all read through it, so a key means
the same thing wherever it appears and adding one happens once.

A line qualifies when, after leading blanks and tabs, it begins with `#`, a key
of upper-case letters, digits or underscores, then `:`. Position carries no
meaning — there is no header block and nothing bounds the scan. Requiring an
upper-case key is what keeps ordinary prose out: `# note: rerun weekly` is a
comment, `#DEP: star/2.7.11b` is a declaration.

The consequence is deliberate and worth stating: a recipe that writes a job
script in a heredoc also declares whatever that script declares, because nothing
distinguishes the two without parsing shell.

Scheduler directives are not annotations. `#SBATCH --time=01:00:00` carries
colons in its value, so cutting on the first one would split it; they are
matched by prefix and handed to the scheduler packages verbatim.

`StripComments` is what the keys hash. It removes every whole-line comment from
the recipe — the header along with the rest — so both keys see what the recipe
*does*, and editing a comment moves neither.

Whole-line only, and blank lines stay. A trailing comment cannot be removed
without parsing shell, for the same reason the boundary does not: `echo "a # b"`
and `${v#pre}` both carry a `#` that is not one. Neither function touches what is
stored or what runs — `/.cnt/recipe` is the original bytes, and so is the script
the build executes.

## Validation and lint

`Recipe.Validate` reports headers a recipe of that type may not declare, and
`Recipe.Lint` reports declarations that are legal but probably wrong. Both are
returned, never printed — nothing here writes to a terminal — and both run when a
recipe is fetched for a build rather than while indexing, so one bad recipe stops
its own build instead of taking a whole collection out of every listing.

| rule | kind |
|---|---|
| `#DEP:` on anything but data | error |
| `#ARCH:` on an OS, or a value other than `native`/`noarch` | error |
| a dependency mentioned in the name but not as whole components | lint |

Only data has build dependencies: an app is prebuilt and self-contained, and an
OS is self-contained by definition. A recipe that genuinely needs a compiler is
an `os` artifact providing that toolchain, not an app depending on one.

`HasComponents` is the matching rule the lint is built on: a dependency's
slash-separated components must occur as a contiguous run of the artifact name's.
Components compare as exact strings and versions are never compared semantically,
so `star/2.7.11b` matches `grch38/star/2.7.11b/gencode49` and not
`grch38/star2.7.11b/gencode49`. The second is the near-miss the lint reports — the
author meant the version to be load-bearing, and the `#TARGET:` quietly stopped it
counting.

## What may surprise you

Everything below is a decision, not an accident. A function read on its own will
not reveal any of it.

**Sorting happens once, and only where the separator is visible.** `#PH:` value
lists are ordered by their separator — comma sorts newest-first, pipe keeps the
author's order — because `values[0]` is the default offered to users. A directory
source reads the source form, where the separator is still there, and sorts. An
HTTP source reads the published index, which stores the finished array and
**records no separator**, so it takes the values verbatim and must: a comma list
and a pipe list are indistinguishable once written, so re-sorting would silently
reorder every author-ordered list. The generator and this package therefore have
to agree on one rule, and are pinned to a shared table of cases rather than to
each other's source.

**An unreachable source keeps its place.** Dropping it would let first-wins
silently promote the next source's recipes, and the build would look normal. It
stays in the list carrying its error, lookups skip it, and the caller reports it
once. A catalog where every source failed is empty rather than failing, for the
same reason not-provided is an outcome: empty is a state callers already handle.

**A stale index beats a failed fetch.** An expired cache entry with no route out
is served with a flag set. That is a compute node, not an error.

**The selected source owns its artifact endpoints.** `source.json` may declare
one OCI `push`, ordered `pull` endpoints, and an `audience`. A build retains the
exact source that won first-match lookup and tries only its endpoints. A stale
cached descriptor may still be used, but a pulled artifact must match the
equivalence key derived from that same cached recipe before it is installed.
An invalid descriptor disables its endpoints without hiding its recipes.

**`>=` is for reuse, not for widening.** `samtools/1.23.1>=1.10` admits
`[1.10, 1.23.1]` — the preferred version is the implicit upper bound. The lower
bound exists so an artifact already on disk can satisfy the dep; it never reaches
a fresh solve, because handing micromamba a range would let two machines resolve
one recipe to different versions.

**Installed beats newer.** With no preferred version, the newest *installed*
version in range wins over the newest *available* one, so a bare `#DEP:` does not
rebuild the moment upstream moves.

**Placeholder key order is not stored.** It is derivable from the `#TARGET:`
token order, so storing it would be a second copy free to drift. Value order *is*
stored — it picks the default.

**Expanding a template leaves `Text` alone.** `Expand` sets `Rendered` — the
substituted copy a build runs, reached through `Script()` — and keeps `Text` as
the template it was fetched as, tokens and all. The template is what an artifact
embeds and what a rebuild starts from, so every variant shares one recipe digest
and is told apart by the placeholder values recorded beside it. `PH` carries
those: after expansion it holds exactly one value per name.

**A template match is not a wildcard.** Names are matched against the declared
`#PH:` values; only an explicit `*` is permissive. A half-filled name is not a
resolvable form — it does not determine its remaining values, so there is no
partial fill.

## Build dependencies

`ValidateDeps` holds both rules a *build's* `#DEP:` answers to, and its only
callers are `Recipe.Validate` and `build.FromExternalSource` — never `run` or
`check`. It lives here rather than on `Recipe` because an external build
(`create -p -f <script>`) answers to the same rules and has no recipe path to be
parsed from.

**Only a `data` recipe may declare one.** An `app` is self-contained — a conda
env, or a prebuilt package carrying its own libraries — and an `os` is
self-contained by definition. Producing an index needs the producing tool,
which is why data is the type with deps. The dep may name an app, data or os.

**A build's `#DEP:` is a `name/version`, never a path; a running script's may be
either.** The asymmetry is the point: a running script mounts what it names and
records nothing, while a build's declaration becomes an edge in an artifact that
must mean the same thing on another machine — and a path is neither resolvable
there nor a key anything can regenerate. So `#DEP: ./overlays/x.sqf` stays valid
in an analysis script and is rejected in a recipe.

`IsPathDep` is the single answer to "is this dep a path", because this package
owns the `Normalize`/`ParseDep` grammar that applies only to names. Its extension
set must stay `utils.IsOverlay`'s, `.ext3` included, plus `.sif` — a root-only,
literal-path reference, never an overlay: it never stacks and carries no build
identity, so `ValidateDeps` refuses it in a recipe the same as any other path,
and it is only ever useful in a running script's own `#DEP:`.

## `#TARGET:` names the artifact

Without a `#PH:` beside it, `#TARGET:` is not a template — it is how an external
build gets a name at all. Absent it, the name is the `-p` basename, a single
component, so `key.Role`'s component match never fires and every `#DEP:` is
silently downgraded to build history. `FromExternalSource` therefore refuses an
external script that declares a dep without one.

The name and the file path are separate namespaces: `#TARGET:` fixes
`/cnt/<name>` and the role classification, `-p` fixes where the `.sqf` lands.
That is why a path-addressed artifact's filename carries no naming claim —
`project.LookupAt` matches it by manifest name alone — while a flat or store scan
still requires the filename to encode the name, because there the filename *is*
the address.

## Invariants worth not breaking

- One spelling of a name. `Normalize` decides what a name *is*; index keys use
  that form, so a second implementation elsewhere would make one string resolve
  two ways. Same for `CompareVersions` and constraint satisfaction.
- The two backends must return identical entries for the same collection. A test
  builds one on disk, serves it over HTTP, and diffs the results — that test is
  the contract, not the prose.
- The caller owns the cache directory and TTL; the package owns everything inside
  it.
