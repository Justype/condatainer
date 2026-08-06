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

**Expanding a template rewrites the body, and leaves the headers.** The recipe
hash is taken over expanded text, so an unexpanded template would give every
variant one key. `#PH:`/`#TARGET:` survive in the output because they are
comments to `bash` and the recipe is embedded in the artifact for provenance.

**A template match is not a wildcard.** Names are matched against the declared
`#PH:` values; only an explicit `*` is permissive. A half-filled name is not a
resolvable form — it does not determine its remaining values, so there is no
partial fill.

## Invariants worth not breaking

- One spelling of a name. `Normalize` decides what a name *is*; index keys use
  that form, so a second implementation elsewhere would make one string resolve
  two ways. Same for `CompareVersions` and constraint satisfaction.
- The two backends must return identical entries for the same collection. A test
  builds one on disk, serves it over HTTP, and diffs the results — that test is
  the contract, not the prose.
- The caller owns the cache directory and TTL; the package owns everything inside
  it.
