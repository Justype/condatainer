# Custom Data Build Scripts

Any dataset your work depends on can be an overlay: a search database, a model checkpoint, an annotation release, a pretrained embedding, a genome index. They share the same three problems on HPC:

- The dataset is many thousands of small files (hostile to inode quotas).
- Everyone keeps a private copy.
- Six months later nobody can say which release or which tool version produced it.

A data recipe fixes all three. The dataset becomes one `.sqf` file (**one inode**), installed and referred to by name, with the inputs that produced it recorded in that name.

Data recipes are [app recipes](./custom_app.md) with two differences: they're named by dataset rather than by tool, and they usually **depend on other overlays** instead of downloading everything themselves.

See the [Build Script Manual](../manuals/build_script.md#data) for the complete reference.

## Naming

```
<collection>/<datatype>/<version>
```

`<collection>` is whatever groups the data — a genome assembly, a project, a database name, a model family. The file path under `recipes/` *is* the name:

| Path | Installs as |
|---|---|
| `grch38/genome/gencode` | GRCh38 primary assembly FASTA |
| `grch38/gtf-gencode/47` | GENCODE 47 annotation |
| `grch38/star/2.7.11b/gencode47-101` | STAR index — tool version **and** annotation in the name |

```{important}
If the data was **derived** by a tool rather than downloaded as-is, put the tool version in the name. Derived artifacts — search indexes, quantized weights, preprocessed caches — are usually not portable across versions of the tool that wrote them, and `grch38/star/gencode47` gives you no way to tell which one you have.
```

## Example 1: Downloaded Data

The simplest shape, and the same whatever the payload is — fetch it, prepare it, expose the path. The example below downloads a genome FASTA; a model checkpoint or a database dump differs only in the URL and the preparation step.

```bash
#!/usr/bin/env bash
#DESC:GRCh38 primary assembly genome FASTA, GENCODE naming
#URL:https://www.gencodegenes.org/human/

#DEP:samtools/1.23.1>=1.10
#SOURCE:genome https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz

#ENV:GENOME_FASTA={prefix}/GRCh38.primary_assembly.genome.fa   ## GRCh38 primary assembly FASTA

cd "$CNT_PREFIX"

echo "Decompressing the genome" >&2
pigz -dc -p "${NCPUS:-4}" "$CNT_SRC_genome" > GRCh38.primary_assembly.genome.fa

echo "Indexing with samtools faidx" >&2
samtools faidx GRCh38.primary_assembly.genome.fa
```

Four things to copy:

1. **`cd "$CNT_PREFIX"` first.** For data you usually want to write directly into the target instead of staging in `$CNT_TMP` and copying a multi-gigabyte file twice.
2. **`#DEP:samtools/1.23.1>=1.10`** — whatever tool you need to unpack, validate, or index the download, declare it. CondaTainer installs it and puts it on `$PATH` before running the recipe; don't assume it exists on the node and don't `module load` it. The `>=1.10` makes any release in `[1.10, 1.23.1]` acceptable, so an existing install is reused instead of pulling another copy.
3. **`#SOURCE:` for the download, decompressed at build time.** The file arrives read-only as `$CNT_SRC_genome`, and its SHA-256 goes into the overlay's identity, so a re-released upstream file is a different overlay instead of a silent duplicate. Decompress into the payload with `pigz -dc` rather than keeping the `.gz`: the overlay is compressed anyway, and a `.gz` inside it costs two decompressions on every read.
4. **`#ENV:GENOME_FASTA`** — the whole point. Consumers get the path without knowing the filename or layout, and, as the next example shows, so do other recipes.

## Example 2: Derived Data, as a Template

When a tool has to *process* the data — build a search index, quantize weights, precompute a cache — two things change: the recipe consumes other overlays instead of downloading, and the work needs a scheduler job.

### Headers

```bash
#!/usr/bin/env bash
#DESC:STAR {star_version} index for GRCh38 GENCODE {gencode_version}, read length {read_length}
#URL:https://github.com/alexdobin/STAR
#TARGET:grch38/star/{star_version}/gencode{gencode_version}-{read_length}

#PH:star_version:2.7.11b,2.7.11a,2.7.10b,2.7.10a,2.7.9a,2.7.8a
#PH:gencode_version:49,48,47,46,45,44,43,42,41,40
#PH:read_length:101,151,*

#AUTOUPDATE:star_version:conda bioconda::star min=2.7.8a
#AUTOUPDATE:gencode_version:dep from=grch38/gtf-gencode

#DEP:star/{star_version}
#DEP:grch38/genome/gencode
#DEP:grch38/gtf-gencode/{gencode_version}

#ENV:STAR_INDEX_DIR={prefix}   ## STAR index, GRCh38 GENCODE {gencode_version}, read length {read_length}

#SBATCH --cpus-per-task=16
#SBATCH --mem=42G
#SBATCH --time=2:00:00
#SBATCH --job-name=star-index
```

One file covers hundreds of concrete outputs — a template is worth it whenever the same recipe is re-run across tool versions, data releases, or a build parameter. Three placeholder styles appear here:

- `star_version` — an explicit comma list, sorted newest-first, kept current by `#AUTOUPDATE:` from a package channel. Use this for **the tool that produces the data**.
- `gencode_version` — a list of **numbered upstream releases**, kept in step with the annotation recipe it depends on. A pure-integer range such as `22-49` expands to every value in between.
- `read_length:101,151,*` — the `*` makes it **open-ended**: 101 and 151 are suggested, but any value is accepted. Use this for **build parameters** that are a property of the user's data, not a fixed catalog.

Placeholders flow into `#DEP:` too, so the dependency graph is version-correct per expansion: building `gencode47-101` pulls `grch38/gtf-gencode/47`, not some other release.

### Consuming Your Dependencies

```bash
cd "$CNT_TMP"

echo "Building STAR index for $CNT_NAME" >&2
STAR --runThreadN "$NCPUS" \
    --runMode genomeGenerate \
    --genomeDir "$CNT_PREFIX" \
    --genomeFastaFiles "$GENOME_FASTA" \
    --sjdbGTFfile "$ANNOTATION_GTF" \
    --sjdbOverhang $(( {read_length} - 1 ))
```

Note what is **not** here: no download, no paths. `$GENOME_FASTA` and `$ANNOTATION_GTF` come from the `#ENV:` lines of the two `#DEP:` overlays, already mounted. This is why `#ENV:` matters — it's the interface between recipes, not just a convenience for end users. Your recipe names *what* it needs, never where it lives.

The rest is ordinary care:

- **`$NCPUS`** is set from your `#SBATCH --cpus-per-task`, so the same recipe is correct on a laptop and on a 16-core allocation.
- **`{read_length}`** is substituted textually before bash sees it — `$(( {read_length} - 1 ))` becomes `$(( 101 - 1 ))`.

### Scheduler Directives

Plain downloads can run where you launch `create`; **anything compute- or memory-hungry should not**. The `#SBATCH` block is what makes `condatainer create` submit a job instead of hogging a shared login node — 16 cores and 42 GB here, because that's what this particular build genuinely needs.

Directives are translated across schedulers, so PBS and LSF users get the equivalent job from the same `#SBATCH` lines. If your recipe only downloads, omit the block entirely.

## Naming Environment Variables

`#ENV:` names are the contract with everything downstream, so keep them predictable — a consumer that reads `$GENOME_FASTA` should keep working when you swap one overlay for a comparable one:

| Pattern | Use | Example |
|---|---|---|
| Standard name | Well-known data types | `GENOME_FASTA`, `ANNOTATION_GTF` |
| `*_GZ` suffix | The file is compressed | `ANNOTATION_GTF_GZ` |
| `<TOOL>_*_DIR` | The artifact is a directory | `STAR_INDEX_DIR` |
| `<TOOL>_PREFIX` | The artifact is a file prefix | `BOWTIE2_PREFIX` |

Prefix with the tool name whenever the data is only meaningful to that tool. `{prefix}` expands to the image's install prefix at load time. See [ENV Naming Guidelines](../manuals/build_script.md#env-naming-guidelines) for the full list.

## Build and Verify

One command builds the entire graph — CondaTainer resolves the `#DEP:` tree, fetches both upstream datasets, installs the tool, then submits the build job:

```bash
condatainer create grch38/star/2.7.11b/gencode47-101
```

Or let it prompt you through the placeholders:

```bash
condatainer create grch38/star-gencode
```

Check what the overlay exposes:

```bash
condatainer info grch38/star/2.7.11b/gencode47-101
# Environment
#  - STAR_INDEX_DIR=/cnt/grch38/star/2.7.11b/gencode47-101
#    # STAR index for GRCh38 GENCODE v47 with read length 101
```

Then use it from an analysis script — declare the dependency and read the variable:

```bash
#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=42G
#DEP: star/2.7.11b
#DEP: grch38/star/2.7.11b/gencode47-101

STAR --genomeDir "$STAR_INDEX_DIR" --runThreadN "$NCPUS" ...
```

```bash
condatainer check align.sh -a   # install anything missing
condatainer run align.sh        # submit with overlays mounted
```

## Share It

Data is the strongest case for sharing: it's large, it's identical for everyone, and one build serves the whole group.

Put the recipe in a collection that your group's config lists under `sources`, and build the overlay into the shared images directory — group members then get it from `condatainer avail` with no rebuild. See [Configuration](../manuals/configuration.md) for how sources are merged, and [Sharing Your Recipes](../deployment/share_scripts.md) for hosting a collection.

## Related

- [Build Script Manual](../manuals/build_script.md#data) — full header/variable reference
- [Sharing Your Recipes](../deployment/share_scripts.md) — upstreaming or hosting your own collection
- [Custom App Build Scripts](./custom_app.md) — packaging the tools that build the data
- [Module Overlays](../user_guide/module_overlays.md) — installing and using data overlays
- [Custom Bundle Overlays](./custom_bundle.md) — read-only Conda environments
- [Scheduler Integration](../qa/scheduler.md) — how `#SBATCH` directives are handled
