# Reproducible Projects

A project pins the exact overlays its scripts use, so a collaborator, a cluster, or you in two years can run the same analysis with the same tools. The pin is a folder, `cnt-lock/`, that you commit with the code.

This page follows one small project from the first lock to a published, restorable copy. The [Project section of the manual](../manuals/condatainer.md#project) lists every command and flag.

## The example

```
rnaseq/
└── src/
    └── count_genes.sh
```

```bash
#!/bin/bash
#DEP: grch38/gtf-gencode/49
#DEP: openjdk/25.0.2
set -euo pipefail
awk '$3=="gene"' "$ANNOTATION_GTF" | wc -l
```

Each `#DEP:` names an exact `name/version`. A project script never uses a version range.

## Lock it

Install what the scripts need, then lock from the project folder:

```bash
condatainer check src -a       # installs any missing overlay
condatainer project lock
```

```
[CNT] New project: /home/me/rnaseq
[CNT]   base → provenance/ubuntu24--base@2136efcde6c8
[CNT]   grch38/gtf-gencode/49 → provenance/grch38--gtf-gencode--49@8108eb46cf1a
[CNT]   openjdk/25.0.2 → provenance/openjdk--25.0.2@1585b1d6a848
[CNT✓] Every declaration is pinned (3).
```

Every `#DEP:` is pinned to the build you have installed, and so is the project's base image. Lock fails if a declared overlay is not installed.

`cnt-lock/` holds the pins and, for each pinned artifact, its manifest and the recipe that rebuilds it. It holds no data and no absolute paths. Commit it.

```bash
condatainer project list       # what is pinned
condatainer project validate   # is the lock complete and consistent
```

Run `project lock` again whenever a script's `#DEP:` lines change. It adds new pins and drops the ones nothing declares.

Inside the project folder, `condatainer run`, `exec` and `check` mount the pinned builds.

## Identity and equivalence

Every artifact has two keys, and the lock records both.

- **Identity** names one exact build.
- **Equivalence** groups builds that differ only in ways that should not change the result, so any of them can substitute for the others.

Some things change only identity:

- an upstream `#SOURCE:` file was re-released under the same name;
- a Conda solve picked the same package versions with other build strings;
- a dependency was rebuilt as an equivalent build.

Editing a recipe's commands, or choosing another `#PH:` value, changes both. The [manual](../manuals/condatainer.md#identity-and-equivalence) lists every case.

The project decides which one it requires. The default is `equivalence`: an equivalent build is used, and noted where it is. To require the exact builds:

```bash
condatainer project select-match identity
```

This is stored in the lock, so commit it. From then on `restore` refuses a substitute, and `run` and `exec` stop when only a substitute is installed.

## Restore it somewhere else

On another machine, check out the project and ask what a restore would do:

```bash
git clone <your repository> && cd rnaseq
condatainer project restore --dry-run
```

On a machine with nothing installed:

```
[CNT]   build       ubuntu24/base → store or flat
[CNT]   build       grch38/gtf-gencode/49 → store or flat
[CNT]   build       openjdk/25.0.2 → store or flat
[CNT✓] 3 artifact(s) to acquire, 0 already available.
```

On the machine that locked it, every step reads `adopt`, because it is already installed:

```
[CNT]   adopt       grch38/gtf-gencode/49 → store or flat
[CNT✓] 0 artifact(s) to acquire, 3 already available.
```

Each artifact takes the first route that works:

| step | what happens |
|---|---|
| `adopt` | it is installed here, as the pinned build or, under `equivalence`, an equivalent one |
| `fetch` | a registry recorded in the lock has it, and it is downloaded |
| `build` | it is rebuilt from the recipe in `cnt-lock/` |

Run `condatainer project restore` to do it. What to expect:

- A rebuild downloads its `#SOURCE:` files again. If an upstream file has changed since, the result is a different build.
- Under `equivalence`, an equivalent result is kept and the line says so: `(equivalent, not sha256:41ab1c2d3e4f; differs: src:gtf)`. Under `identity`, the restore fails instead.
- A recipe with scheduler directives is submitted as a job, and restore exits with a status that says jobs are pending. Run it again after they finish.
- `--no-prebuilt` builds locally instead of downloading, `--keep-build-deps` keeps the dependencies a rebuild needed, and `condatainer project validate --installed` checks the result.

## Publish it

Restoring by rebuilding is slow for large data. Publish the built artifacts to an OCI registry, and restore downloads them instead.

### Record where

```bash
condatainer project registry set ghcr.io/my-lab/rnaseq/cnt
```

```
[CNT✓] Publishing to ghcr.io/my-lab/rnaseq/cnt
[CNT]   audience public
```

The value is `<registry>/<owner>/<repository>`, with no tag. Every artifact goes into that one repository, named by its tag, so the whole project is a single package. The setting is stored in the lock.

**Audience** sets what a push may publish to this registry. It does not control who can pull; that is set at the registry. At the default `public`, an app whose recipe does not allow redistribution is refused. A `restricted` registry takes anything, so choose it only when the registry really is limited to people who may receive everything:

```bash
condatainer project registry set ghcr.io/my-lab/rnaseq/cnt --audience restricted
```

Nothing checks the claim. It is not GitHub's package visibility setting, which CondaTainer never reads or changes.

### Log in

```bash
printf '%s\n' "$TOKEN" | condatainer registry login ghcr.io --username "$USER" --password-stdin
```

For GHCR, a token with `write:packages` can push. A restore from a private package needs `read:packages`.

### Push

Preview first:

```bash
condatainer project push --dry-run
```

```
[CNT] Publishing to ghcr.io/my-lab/rnaseq/cnt (public)
[CNT]   upload   grch38/gtf-gencode/49  grch38--gtf-gencode--49__8108eb46cf1a grch38--gtf-gencode--49
[CNT]   upload   openjdk/25.0.2  openjdk--25.0.2__1585b1d6a848 openjdk--25.0.2
[CNT]            channels: conda-forge, bioconda
[CNT]   upload   ubuntu24/base  ubuntu24--base__2136efcde6c8 ubuntu24--base
[CNT] 3 to upload
```

Then `condatainer project push`. Each upload is recorded in the lock, so commit the lock again afterwards, and a later `restore` downloads instead of rebuilding. What to expect:

- **It never builds.** An artifact that is not installed at its locked identity is reported with the restore that would produce it. The project must also pass `validate`.
- **Artifacts a public collection already serves are skipped** (`upstream`), since restore tries that location first. `--all` uploads them too, for a project that has to outlive the collection.
- **Each artifact gets two tags:** its name, and its name plus identity. The name tag moves when you re-pin; the identity tag never does, which keeps an older commit's lock restorable.
- **The GitHub package is separate from the repository.** Its visibility is its own setting, and CondaTainer never changes it, so check it on GitHub after the first push. Collaborators need to be able to pull the package. `--source <url>` records the repository the package links back to.

### When a push is refused

At a `public` endpoint, whether an artifact may be published depends on what its recipe declares:

```
[CNT]   refused  cellranger/9.0.1  app artifacts are not published to a public registry…
```

| the recipe says | at a public endpoint |
|---|---|
| `#REDISTRIBUTE: no` | never |
| `#REDISTRIBUTE: yes` | published |
| nothing, and it is `data`, `os` or a Conda build | published |
| nothing, and it is an `app` | refused |

An app is refused because it is someone else's software with unstated terms, and unknown is not permission. A refusal stops the whole push before anything uploads. To proceed:

- If you own the recipe, add `#REDISTRIBUTE: yes` (and `#LICENSE:`) when you are allowed to.
- If the registry really has a limited audience, set `--audience restricted`, which lifts these checks.
- Otherwise leave the artifact out of the published project.

A recipe from someone else's collection cannot be declared from your project. Use a `restricted` endpoint for those.

## Related

- [Project, in the manual](../manuals/condatainer.md#project) — every command and flag
- [Identity and equivalence](../manuals/condatainer.md#identity-and-equivalence) — what changes each key
- [Distributing Artifacts](../deployment/distribution.md) — why a project is one package, and how tags keep old locks restorable
- [Publishing rules](../manuals/condatainer.md#publishing-rules) — the full table
