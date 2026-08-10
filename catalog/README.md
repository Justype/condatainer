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

## The header boundary

`headerBoundary` decides where a recipe's headers end and its body begins.
Scanning from offset 0, a line is header when its content — after removing the
line terminator and any leading spaces or tabs — is empty or begins with `#`. The
body starts at the first line that fails that test.

Its only job is bounding what the recipe parser reads, which is what keeps a
`#DEP:` in a heredoc — or a `#SBATCH` in a file the recipe writes — inert. It
feeds no key. Nothing parses shell: a body carries `#` inside strings and
`${x#y}` expansions, so a comment-aware rule could move the boundary on a body
that never changed.

Two consequences are worth stating: a leading `#!/bin/bash` is header, and so is a
blank line before the first command.

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
| `#ARCH:` on an OS or a base, or a value other than `native`/`noarch` | error |
| a dependency mentioned in the name but not as whole components | lint |

Only data has build dependencies: an app is prebuilt and self-contained, an OS is
self-contained by definition, and a base *is* the build environment. A recipe that
genuinely needs a compiler is an `os` artifact providing that toolchain, not an
app depending on one. `Resolve` enforces the other half — a dep that resolves to a
base is refused, though a base may still be a root of the walk, which is how a
base gets built.

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

## Invariants worth not breaking

- One spelling of a name. `Normalize` decides what a name *is*; index keys use
  that form, so a second implementation elsewhere would make one string resolve
  two ways. Same for `CompareVersions` and constraint satisfaction.
- The two backends must return identical entries for the same collection. A test
  builds one on disk, serves it over HTTP, and diffs the results — that test is
  the contract, not the prose.
- The caller owns the cache directory and TTL; the package owns everything inside
  it.
