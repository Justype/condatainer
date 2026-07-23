# Custom Data Build Scripts

Any dataset your work depends on can be an overlay: a search database, a model checkpoint, an annotation release, a pretrained embedding, a genome index. They share the same three problems on HPC:

- The dataset is many thousands of small files (hostile to inode quotas).
- Everyone keeps a private copy.
- Six months later nobody can say which release or which tool version produced it.

A data build script fixes all three. The dataset becomes one `.sqf` file (**one inode**), installed and referred to by name, with the inputs that produced it recorded in that name.

Data scripts are [app build scripts](./custom_app.md) with two differences: they're named by dataset rather than by tool, and they usually **depend on other overlays** instead of downloading everything themselves.

See the [Build Script Manual](../manuals/build_script.md#data) for the complete reference.

## Naming

```
<collection>/<datatype>/<version>
```

`<collection>` is whatever groups the data — a genome assembly, a project, a database name, a model family. The file path under `build-scripts/` *is* the name:

| Path | Installs as |
|---|---|
| `grch38/genome/gencode` | GRCh38 primary assembly FASTA |
| `grch38/gtf-gencode/47` | GENCODE 47 annotation |
| `grch38/star/2.7.11b/gencode47-101` | STAR index — tool version **and** annotation in the name |

```{important}
If the data was **derived** by a tool rather than downloaded as-is, put the tool version in the name. Derived artifacts — search indexes, quantized weights, preprocessed caches — are usually not portable across versions of the tool that wrote them, and `grch38/star/gencode47` gives you no way to tell which one you have.
```

Start from the [`assets/build-template-ref`](https://github.com/Justype/cnt-scripts/blob/main/assets/build-template-ref) template — it carries the same shared boilerplate as the apps template, plus `smart_copy_or_link` for staging large files.

## Example 1: Downloaded Data

The simplest shape, and the same whatever the payload is — fetch it, prepare it, expose the path. The example below downloads a genome FASTA; a model checkpoint or a database dump differs only in the URL and the preparation step.

Source: [`grch38/genome/gencode`](https://github.com/Justype/cnt-scripts/blob/main/build-scripts/grch38/genome/gencode)

```bash
#!/usr/bin/bash
#DEP:samtools/1.23.1>=1.10
#AUTOUPDATE:samtools:bioconda:samtools
#WHATIS:GRCh38 reference genome FASTA
#URL:https://www.gencodegenes.org/human/
#ENV:GENOME_FASTA=$app_root/GRCh38.primary_assembly.genome.fa
#ENVNOTE:GRCh38 reference genome

install_app() {
    cd "$target_dir"

    url='https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz'
    print_stderr "Downloading ${YELLOW}${app_name_version}${NC}"
    curl -fsSL -o "GRCh38.primary_assembly.genome.fa.gz" "$url"

    print_stderr "Extracting ${YELLOW}${app_name_version}${NC}"
    pigz_or_gunzip "GRCh38.primary_assembly.genome.fa.gz"

    print_stderr "Indexing ${YELLOW}${app_name_version}${NC} with samtools faidx"
    samtools faidx "GRCh38.primary_assembly.genome.fa"
}
```

Four things to copy:

1. **`cd "$target_dir"` first.** `install_app()` starts in `$tmp_dir`. For data you usually want to write directly into the target instead of staging in tmp and copying a multi-gigabyte file twice.
2. **`#DEP:samtools/1.23.1>=1.10`** — whatever tool you need to unpack, validate, or index the download, declare it. CondaTainer installs it and puts it on `$PATH` before running the script; don't assume it exists on the node and don't `module load` it. The `>=1.10` makes any release in `[1.10, 1.23.1]` acceptable, so an existing install is reused instead of pulling another copy.
3. **`pigz_or_gunzip`** instead of plain `gunzip` — uses all `$NCPUS` cores when `pigz` is present. Worth it on any large archive.
4. **`#ENV:GENOME_FASTA`** — the whole point. Consumers get the path without knowing the filename or layout, and, as the next example shows, so do other build scripts.

## Example 2: Derived Data, as a Template

When a tool has to *process* the data — build a search index, quantize weights, precompute a cache — two things change: the script consumes other overlays instead of downloading, and the work needs a scheduler job.

Source: [`grch38/star-gencode`](https://github.com/Justype/cnt-scripts/blob/main/build-scripts/grch38/star-gencode)

### Headers

```bash
#!/usr/bin/bash
#PL:star_version:2.7.0b,...,2.7.11a,2.7.11b
#AUTOUPDATE:star_version:bioconda:star>=2.7.0b
#PL:gencode_version:22-49
#PL:read_length:101,151,*
#TARGET:grch38/star/{star_version}/gencode{gencode_version}-{read_length}

#DEP:grch38/genome/gencode
#DEP:grch38/gtf-gencode/{gencode_version}
#DEP:star/{star_version}

#WHATIS:STAR GRCh38 GENCODE{gencode_version} index for read length {read_length}
#URL:https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
#ENV:STAR_INDEX_DIR=$app_root
#ENVNOTE:STAR index for GRCh38 GENCODE v{gencode_version} with read length {read_length}

#SBATCH --cpus-per-task=16
#SBATCH --mem=42G
#SBATCH --time=2:00:00
#SBATCH --job-name=star-index
#SBATCH --output=%x-%j.log
```

One file covers hundreds of concrete outputs — a template is worth it whenever the same recipe is re-run across tool versions, data releases, or a build parameter. Three placeholder styles appear here:

- `star_version` — an explicit comma list, sorted newest-first, kept current by `#AUTOUPDATE:` from a package channel. Use this for **the tool that produces the data**.
- `gencode_version:22-49` — a pure-integer **range**, expanded to every value in between. Use this for **numbered upstream releases**.
- `read_length:101,151,*` — the `*` makes it **open-ended**: 101 and 151 are suggested, but any value is accepted. Use this for **build parameters** that are a property of the user's data, not a fixed catalog.

Placeholders flow into `#DEP:` too, so the dependency graph is version-correct per expansion: building `gencode47-101` pulls `grch38/gtf-gencode/47`, not some other release.

### Consuming Your Dependencies

```bash
install_app() {
    if [ -z "$ANNOTATION_GTF" ]; then
        pigz_or_gunzip_pipe "$ANNOTATION_GTF_GZ" > annotation.gtf
        ANNOTATION_GTF="annotation.gtf"
    fi

    print_stderr "Building STAR index for ${YELLOW}${app_name_version}${NC}"
    STAR --runThreadN $NCPUS \
        --runMode genomeGenerate \
        --genomeDir "$target_dir" \
        --genomeFastaFiles "$GENOME_FASTA" \
        --sjdbGTFfile "$ANNOTATION_GTF" \
        --sjdbOverhang $(({read_length} - 1))
}
```

Note what is **not** here: no download, no paths. `$GENOME_FASTA` and `$ANNOTATION_GTF_GZ` come from the `#ENV:` lines of the two `#DEP:` overlays, already mounted. This is why `#ENV:` matters — it's the interface between build scripts, not just a convenience for end users. Your script names *what* it needs, never where it lives.

The rest is ordinary care:

- **Do the expensive step only if needed.** The `-z "$ANNOTATION_GTF"` guard reuses an uncompressed copy when an upstream overlay already provides one, and only falls back to unpacking the `.gz`.
- **`$NCPUS`** is set from your `#SBATCH --cpus-per-task`, so the same script is correct on a laptop and on a 16-core allocation.
- **`{read_length}`** is substituted textually before bash sees it — `$(({read_length} - 1))` becomes `$((101 - 1))`.

### Scheduler Directives

Plain downloads can run on the login node; **anything compute- or memory-hungry must not**. The `#SBATCH` block is what makes `condatainer create` submit a job instead of hogging a shared login node — 16 cores and 42 GB here, because that's what this particular build genuinely needs.

Directives are translated across schedulers, so PBS and LSF users get the equivalent job from the same `#SBATCH` lines. If your script only downloads, omit the block entirely.

## Naming Environment Variables

`#ENV:` names are the contract with everything downstream, so keep them predictable — a consumer that reads `$GENOME_FASTA` should keep working when you swap one overlay for a comparable one:

| Pattern | Use | Example |
|---|---|---|
| Standard name | Well-known data types | `GENOME_FASTA`, `ANNOTATION_GTF` |
| `*_GZ` suffix | The file is compressed | `ANNOTATION_GTF_GZ` |
| `<TOOL>_*_DIR` | The artifact is a directory | `STAR_INDEX_DIR` |
| `<TOOL>_PREFIX` | The artifact is a file prefix | `BOWTIE2_PREFIX` |

Prefix with the tool name whenever the data is only meaningful to that tool. `$app_root` expands to the overlay's mount path at load time. See [ENV Naming Guidelines](../manuals/build_script.md#env-naming-guidelines) for the full list.

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

Put the script in a `build-scripts/` directory under `$CNT_ROOT` (app root) or `$CNT_EXTRA_ROOT` (group layer) and the resulting `.sqf` in the matching `images/` directory — group members then get it from `condatainer avail` with no rebuild. For a shared script repo or a PR upstream, see [Sharing Your Scripts](../deployment/share_scripts.md).

## Related

- [Build Script Manual](../manuals/build_script.md#data) — full header/variable reference
- [Sharing Your Scripts](../deployment/share_scripts.md) — upstreaming or hosting your own source
- [Custom App Build Scripts](./custom_app.md) — packaging the tools that build the data
- [Module Overlays](../user_guide/module_overlays.md) — installing and using data overlays
- [Custom Bundle Overlays](./custom_bundle.md) — read-only Conda environments
- [Scheduler Integration](../qa/scheduler.md) — how `#SBATCH` directives are handled
