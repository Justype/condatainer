# Build Script Manual

This document gives instructions on how to create your own build scripts for **CondaTainer**.

## Table of Contents

- [Naming Conventions](#naming-conventions)
- [Available Variables](#available-variables)
- [Headers](#headers)
  - [Description and URL](#description-and-url)
  - [Set Dependencies](#set-dependencies)
  - [Scheduler Parameters](#scheduler-parameters)
  - [Type Tag](#type-tag)
  - [Template Tags](#template-tags)
  - [Auto-Update Tag](#auto-update-tag)
  - [Environment Variables](#environment-variables) and [ENV Naming Guidelines](#env-naming-guidelines)
  - [Source Tag](#source-tag)
  - [Input Tag](#input-tag)
- [Apps](#apps)
- [Data](#data)
- [OS](#os)

## Naming Conventions

The file path must follow the naming convention below to be recognized by CondaTainer:

`recipes/<name_conversion>` inside a collection (no `.sh` suffix)

Where `<name_conversion>` is defined as:

### OS

Apptainer definition files for distro-level system tools.

* **Format:** `<distro>/<name>`
* **Structure:**
  * **distro**: The base OS distribution (e.g., `ubuntu24`).
  * **name**: The tool name (e.g., `build-essential`, `r4.4.3`).
* **Example:** `ubuntu24/build-essential`

### Apps

Apps not available as conda packages, or specific versions not in conda.

**Single version per file:**

* **Format:** `<name>/<version>`
* **Example:** `cellranger/9.0.1`

**Template (multiple versions, one file):**

* **Format:** `<name>` (a single file with `#PH:` and `#TARGET:` headers)
* **Example:** `cytoscape` → expands to `cytoscape/3.10.3`, `cytoscape/3.10.4`, etc.
* Use this when the install logic is identical across versions and only the download URL changes.

### Data

Any data, including genome reference indexes.

* **Format:** `<assembly|project>/<datatype>/<version>`
* **Structure:**
  * **assembly/project**: The genome assembly or project name (e.g., `grch38`).
  * **datatype**: The type of data (e.g., `gtf-gencode`).
  * **version**: The release or build version (e.g., `47`).
* **Example:** `grch38/genome/gencode`

**Template (multiple versions, one file):**

* **Format:** `<assembly|project>/<datatype>` (a single file with `#PH:` and `#TARGET:` headers)
* **Example:** `grch38/star-gencode` → expands to `grch38/star/2.7.11b/gencode47-101`, `grch38/star/2.7.11b/gencode50-151`, etc.
* Use this when the install/create logic is identical across versions.

### Example

- `cellranger/9.0.1` (App)
- `cytoscape` (App template)
- `grch38/cellranger/2024-A` (Data)
- `grch38/star-gencode` (Data templdate)

## Available Variables

A recipe runs top to bottom as `bash -euo pipefail`. There is no wrapper function and no helper functions: call `curl`, `tar` and `pigz` directly, and print progress with `echo ... >&2`.

| Variable | Description |
| -------- | ----------- |
| `CNT_NAME` | The complete module name, e.g. `cellranger/9.0.1` |
| `CNT_TYPE` | `app`, `data` or `os` |
| `CNT_PREFIX` | Where the payload goes. Exactly this directory is packed into the overlay |
| `CNT_TMP` | Temporary working directory, also exported as `TMPDIR` (managed by **CondaTainer**) |
| `CNT_SRC_<name>` | A file fetched by a [`#SOURCE:`](#source-tag) header, read-only |
| `NCPUS` | Number of CPUs (from script directives, or `build.ncpus` config setting) |
| `MEM` | Memory per task in MB (from script directives, or `build.mem` config setting) |
| `MEM_GB` | Memory per task in GB (integer) |

During a script build, `micromamba` is on `PATH` and `MAMBA_ROOT_PREFIX` is the build's scratch, with `MAMBA_NO_RC=true` so your `.condarc` is ignored. `$CNT_PREFIX` already exists as an empty directory, so a recipe that installs an environment there first creates `$CNT_PREFIX/conda-meta/history` (see [Custom Bundle Overlays](../advanced_usage/custom_bundle.md)); micromamba refuses a directory that is not an environment. A script that solves an environment with it is identified by its text, not by the packages it resolves, so pin versions.

> `$CNT_TMP` paths:
>
> - `app` module, `.def` and the base image: fast local scratch — `$CNT_TMPDIR` → scheduler tmp → `$TMPDIR` → `/tmp`
> - `data` module: the stable condatainer data-dir tmp
> - A `.def` is always on the fast root: it needs `--fakeroot`, which NFS, Lustre, GPFS and PanFS do not support
> - External build (`-f`): `app` and `.def` use fast scratch; `data` builds next to the target dir
> - `$CNT_TMPDIR` selects the fast root only — it never moves a build off the stable root or off the target dir

## Headers

Headers are special comments at the beginning of build scripts that provide metadata and instructions for CondaTainer.

**Example Header**: `grch38/bowtie2/ucsc_no_alt`

```bash
#!/usr/bin/env bash
#DEP:bowtie2/2.5.5>=2.3
#AUTOUPDATE:bowtie2:bioconda:bowtie2
#DEP:grch38/genome/ucsc_no_alt

#DESC:bowtie2 index for GRCh38 UCSC no alt reference genome
#URL:https://bowtie-bio.sourceforge.net/bowtie2/index.shtml

#ENV:BOWTIE2_PREFIX={prefix}/GRCh38_no_alt_analysis_set   ## GRCh38 reference genome

#SBATCH --cpus-per-task=12
#SBATCH --mem=12G
#SBATCH --time=2:00:00
#SBATCH --job-name=bowtie2-build
#SBATCH --output=%x-%j.log

bowtie2-build --threads "$NCPUS" "$GENOME_FASTA" "$CNT_PREFIX/GRCh38_no_alt_analysis_set"
```

### Description and URL

`#DESC:` is a one-line description of what the overlay provides; `#URL:` points to its homepage or documentation.

**CondaTainer** displays `#DESC:` in:

- `condatainer avail --description` — shows the description beside each build script.
- `condatainer info <overlay>` — shows the description of a built overlay. This is read from the `#DESC:` line of the build script embedded inside the overlay (a sidecar `.env`, if present, overrides it).
- `condatainer create` — shown during interactive template resolution.

The description supplies the generated module's description (Lmod calls this `whatis`), while the
URL supplies its help link.

### Set Dependencies

`#DEP:` lines specify dependencies that must be installed before building the current overlay.

When **CondaTainer** processes the build script, it will ensure that all specified dependencies are available and load them in the same order as listed.

**Only a `#TYPE:data` recipe may declare `#DEP:`.** A build that finds one on an
app or an OS layer stops with an error. A `#DEP:` is what must be
*mounted while the recipe runs*, and only a dataset genuinely needs that —
producing an index requires the tool that produces it. An app is prebuilt and
self-contained (a Conda environment, or a package carrying its own libraries), and an
OS layer is self-contained by definition. A recipe that genuinely needs a
compiler should be an OS layer providing that toolchain, not an app depending
on one.

A `#DEP:` may name an app, a dataset, or an OS layer.

**In a build script a `#DEP:` is a `name/version`, never an overlay path.** `#DEP:overlays/tool.sqf`
and `#DEP:env.img` are rejected. A build's declaration becomes an edge recorded in the finished
artifact as a name plus a complete identity, and a path is neither — nothing could re-resolve it on
another machine, and no key could be regenerated from it.

This is the one `#DEP:` rule that does **not** apply to a script you *run*. A running script mounts
what it names and records nothing, so `#DEP:./overlays/tool.sqf` and `#DEP:env.img` are
perfectly valid there. Both rules above are enforced only when something is being built — including
an external build (`condatainer create -p <path> -f <script>.sh`), which answers to them exactly as a
catalog recipe does.

**Basic (exact version):**

```bash
#DEP:samtools/1.21
```

Requires exactly `samtools/1.21`. If not installed, builds it.

**With version constraint:**

```bash
#DEP:samtools/1.22.1>=1.10
```

- **Preferred version**: `1.22.1` — built if no compatible version is installed.
- **Minimum version** (`>=1.10`): any installed version in the range `[1.10, 1.22.1]` is accepted.
- **Upper bound**: the preferred version acts as an implicit upper bound. A version higher than `1.22.1` (e.g., `2.0`) will not be used and the preferred version will be built instead.

`>` is also supported for a strict lower bound:

```bash
#DEP:samtools/1.22.1>1.10
```

**Version format**: partial versions are accepted — `1`, `1.10`, `1.22.1`. Missing components are treated as `0` (so `1.10` matches `1.10.0`, `1.10.3`, etc.).

When multiple installed versions satisfy the constraint, the latest one is used.

### Scheduler Parameters

Scheduler directive lines (`#SBATCH`, `#PBS`, or `#BSUB`) allow you to specify job parameters for the build process. (HTCondor uses native `.sub` submit files, so its directives is not supported.)

If a supported scheduler is available, **CondaTainer** will submit the build job with the specified parameters.

**Slurm Example:**

```bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=42G
#SBATCH --time=2:00:00
#SBATCH --job-name=star-index
#SBATCH --output=%x-%j.log
```

- `--cpus-per-task`, `--mem`, and `--time` should be set according to the expected resource requirements.
- `--nodes`, `--ntasks`: **must not be set** (always single-task).
- `--output`: will always be overwritten to point to the `logs` directory.

**PBS Example:**

```bash
#PBS -l select=1:ncpus=16:mem=42gb
#PBS -l walltime=2:00:00
#PBS -N star-index
```

**LSF Example:**

```bash
#BSUB -n 16
#BSUB -R "span[hosts=1] rusage[mem=2688MB]"
#BSUB -W 2:00
#BSUB -J star-index
```

- `span[hosts=1]` is **required** — it keeps all 16 CPUs on a single node. Without it, `-n 16` means 16 MPI tasks and will be randomly distributed.
- `rusage[mem=]` is **per slot** (here, per CPU): `2688MB` × 16 = 42 GB total.

### Type Tag

`#TYPE:` sets the recipe's `catalog.Type`, which decides the install prefix, the
SquashFS block size, and where the build works.

- `#TYPE:app` (default): fast local scratch, honouring `$CNT_TMPDIR`.
- `#TYPE:data`: the stable data-dir tmp, or — for an external build — alongside the target prefix.

Only `app` and `data` are accepted, and only on a script recipe: a `.def` is
always `os`, so declaring `#TYPE:` on one is an error.

```bash
#TYPE:app
#TYPE:data
```

An unrecognized value is **not** an error — it falls through to the default,
which is `data` at two or more module-name components and `app` below that. So a
misspelled type silently builds with the wrong install prefix and block size;
the recipe collection's validator rejects the value for exactly that reason.

### Architecture Tag

`#ARCH:` declares that a payload runs anywhere, rather than only on the
architecture it was built on. It applies to script recipes of either type; a
Conda app and an OS layer may not declare it.

```bash
#ARCH:noarch     ## jar files; the JRE comes from the runtime
```

Two values: `native` — the default, so omit the header for anything
architecture-specific — and `noarch`.

The default is strict on purpose. Every mainstream HPC architecture is 64-bit
little-endian, so most binary indexes *would* in fact port, which is the trap:
getting it wrong does not crash, it silently returns wrong answers. Only the
author knows whether a build produced machine code, a tool-specific binary dump,
or bytes that mean the same everywhere.

- **Portable**: a Java tool shipped as jars, a pure-Python one, a FASTA, a GTF, a VCF.
- **Not portable**: anything compiled, a STAR index, a Kraken2 hash table.

The value lands in the image's recorded platform and is checked when the image is
mounted. Adding `#ARCH:noarch` to a recipe later does not invalidate images
already built from it — they keep claiming `native`, which is honest, and the new
claim applies from the next build onward.

### Licence and Redistribution

Two headers say what may be done with the artifact a recipe produces. Both are
recorded in the image's manifest and travel with every copy of it.

```bash
#LICENSE: MIT
#REDISTRIBUTE: yes
```

`#REDISTRIBUTE:` takes `yes` or `no` and is the **answer** to whether the built
artifact may be published to a registry anyone can pull from. Anything else is a
validation error; omitting it leaves the question unanswered, and the recipe's
`#TYPE:` supplies the default:

| the recipe declared | may publish publicly |
|---|---|
| `#REDISTRIBUTE: no` | **never**, whatever the type, and no flag overrides it |
| `#REDISTRIBUTE: yes` | yes, whatever the type |
| nothing, and it is `os` | yes — a container root, and packages from a public distribution |
| nothing, and it is `data` | yes — the type asserts public reference data and the indexes built from it |
| nothing, and it is an `app` | **no** — someone else's software with unstated terms |

So an `app` recipe that installs a vendor tarball is refused at a public endpoint
until someone reads the licence and writes the answer down. Declaring it here
rather than at push time is the point: it is authored once, reviewed in a commit,
and applies to every artifact ever built from the recipe, including on someone
else's machine. Nothing verifies it.

A vendor download that requires a login or an accepted EULA is the usual reason
to write `#REDISTRIBUTE: no`.

`#LICENSE:` is an [SPDX expression](https://spdx.org/licenses/) — `MIT`,
`GPL-3.0-only`, `Apache-2.0 OR MIT`, `LicenseRef-10x-Genomics-EULA` — published
verbatim as the image's `org.opencontainers.image.licenses` annotation. It is
**never parsed** and gates nothing: deriving redistribution permission from a
licence expression is a judgement a tool gets wrong in the permissive direction,
and that direction cannot be taken back. It is documentation for whoever has to
make that call.

Neither header affects the artifact's identity — like every whole-line comment,
they are stripped before the recipe is hashed — so adding one to a recipe does
not invalidate images already built from it. Their values are read from each
image's own manifest at push time.

See [Distributing Artifacts](../deployment/distribution.md) for how this plays
out across a collection's endpoint and a project's own.

### Template Tags

`#PH:` and `#TARGET:` turn a single script into a parameterized **template**. Running `condatainer create <template-path>` prompts the user for each placeholder value and then creates the concrete overlay.

- `#PH:<name>:<values>` — declares a placeholder. The separator sets the ordering:
  - `,` — comma list, **sorted descending** (use for versions): `3.10,3.9,3.8`
  - `|` — pipe list, **author's order preserved** (use for labels): `kasm | turbo`

  In either form:
  - Integer ranges (`a-b`) expand to every integer between `a` and `b`, inclusive (`22-25` → `22,23,24,25`). Only pure digits are a range, so `2024-A` stays literal.
  - A `*` entry makes the list open-ended (any value accepted); it is always listed last.
  - Duplicates are dropped and whitespace around entries is ignored.
- `#TARGET:<pattern>` — the **module path** of the resulting overlay, once the placeholders are filled.
  Use `{name}` tokens matching the `#PH:` names.
  This is the name users refer to from then on — what they pass to `condatainer create`/`build`, what `#DEP:` lines point at, and the modulefile path [ModGen](https://github.com/Justype/condatainer/blob/main/assets/modgen/manual.md) generates. It is independent of the script's own file name.

Every `#PH:` name must appear as a `{name}` token in `#TARGET:`, and every token must have a matching `#PH:`. Otherwise **CondaTainer** warns and skips the expansion, because unused placeholders would silently collapse every value onto the same target.

#### `#TARGET:` in an external build

`#TARGET:` also works **without** `#PH:`, and there it is not a template — it simply names the artifact. This matters for an external build (`condatainer create -p <path> -f <script>.sh`), where the two are otherwise decided by different things:

| | comes from | governs |
|---|---|---|
| artifact name | `#TARGET:`, else the `-p` basename | the payload's `/cnt/<name>` prefix, and which `#DEP:` lines count toward equivalence |
| file location | `-p` | where the `.sqf` is written |

They are deliberately unrelated: an external overlay is mounted by its path, so its filename makes no claim about its name. `create -p ./overlays/idx -f build.sh` with `#TARGET:star/2.7.11b/index` writes `./overlays/idx.sqf` whose payload lives at `/cnt/star/2.7.11b/index`.

Declaring the name matters because an app or OS `#DEP:` counts toward the artifact's equivalence only when its components appear in the artifact's own name. `star/2.7.11b/index` contains `star/2.7.11b`, so that dependency's version is binding; a basename like `idx` contains nothing, so every dependency would silently be treated as build history instead. That is why **an external script that declares `#DEP:` must also declare `#TARGET:`** — otherwise the classification would depend on the path someone happened to type, and `create` refuses it.

A `{placeholder}` is rejected here: an external build has no `#PH:` declarations and no requested name to select values from, so the pattern could never be filled in.

**Example** — STAR index template `grch38/star-gencode`:

```bash
#!/usr/bin/env bash
#PH:star_version:2.7.0b,...,2.7.11a,2.7.11b
#PH:gencode_version:22-49
#PH:read_length:101,151,*
#TARGET:grch38/star/{star_version}/gencode{gencode_version}-{read_length}

#DEP:grch38/genome/gencode
#DEP:grch38/gtf-gencode/{gencode_version}
#DEP:star/{star_version}

#DESC:STAR GRCh38 GENCODE{gencode_version} index for read length {read_length}
#URL:https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

#ENV:STAR_INDEX_DIR={prefix}   ## STAR index for GRCh38 GENCODE v{gencode_version} with read length {read_length}

#SBATCH --cpus-per-task=16
#SBATCH --mem=42G
#SBATCH --time=2:00:00
#SBATCH --job-name=star-index

STAR --runThreadN "$NCPUS" --runMode genomeGenerate --genomeDir "$CNT_PREFIX" \
    --genomeFastaFiles "$GENOME_FASTA" --sjdbGTFfile "$ANNOTATION_GTF" \
    --sjdbOverhang $(( {read_length} - 1 ))
```

When the user runs `condatainer create grch38/star-gencode`, **CondaTainer** shows the `#DESC:` description, the `#TARGET:` pattern (with `{placeholder}` tokens highlighted), then prompts for each placeholder in declaration order:

```
[CNT] Placeholder template: grch38/star-gencode
[CNT] STAR GRCh38 GENCODE{gencode_version} index for read length {read_length}
  Target: grch38/star/{star_version}/gencode{gencode_version}-{read_length}
  star_version [2.7.0b-2.7.11b] (default: 2.7.11b):
  gencode_version [22-49] (default: 49): 47
  read_length [suggested: 151, 101, or any value] (default: 151): 
  → Creating grch38/star/2.7.11b/gencode47-151
```

If the user already has a compatible dependency installed (e.g. `star/2.7.10` is installed), the default for `star_version` will be `2.7.10` instead of the latest available `2.7.11b`.

See `grch38/star-gencode` for a real template example.

User can directly use the target to fill the placeholders:

```bash
condatainer create grch38/star/2.7.11b/gencode47-101
# Will directly set:
#   star_version=2.7.11b
#   gencode_version=47
#   read_length=101
```

### Auto-Update Tag

`#AUTOUPDATE:` opts a script into automatic version maintenance. A CI workflow runs twice every month, fetches the latest versions from the specified source, and rewrites the version list in place.

```
#AUTOUPDATE:{key}:{source}:{identifier}[>={min}][<{max}|<={max}]
```

`{key}` must match an existing `#PH:`, `#DEP:`, or (in helpers) `#VALUE:` header in the same file. The target type is detected automatically:

| Header matched | Behavior |
|---|---|
| `#PH:key:` | Rewrites the full version list (all versions ≥ min) |
| `#DEP:key/` | Rewrites the pinned version to latest only; preserves `>=constraint` |

**Supported sources:**

| Source | Format | Example |
|---|---|---|
| `github` | `github:{org}/{repo}` | `github:cytoscape/cytoscape>=3.9.0` |
| `bioconda` | `bioconda:{package}` | `bioconda:star>=2.7.0b` |
| `conda-forge` | `conda-forge:{package}` | `conda-forge:r-base>=4.0.0` |
| `docker` | `docker:{image}:{tag_regex}` | `docker:posit/r-base:^(\d+\.\d+\.\d+)-noble(?:-[^-]+)?$>=4.0.0` |

The `docker` source requires a full capture-group regex as the tag pattern — the first capture group is extracted as the version string.

**Examples:**

```bash
# PH template — full version list, all 3.9+
#PH:cytoscape_version:3.9.0,3.9.1,3.10.0,3.10.3,3.10.4
#AUTOUPDATE:cytoscape_version:github:cytoscape/cytoscape>=3.9.0
#TARGET:cytoscape/{cytoscape_version}

# DEP pin — latest samtools; min constraint preserved from the #DEP: line
#DEP:samtools/1.23.1>=1.10
#AUTOUPDATE:samtools:bioconda:samtools

# DEP pin with upper bound — stay on openjdk 17.x, never upgrade to 18+
#DEP:openjdk/17.0.12>=17
#AUTOUPDATE:openjdk:bioconda:openjdk>=17<18
```

Place `#AUTOUPDATE:` immediately after the `#PH:` or `#DEP:` line it manages.

### Environment Variables

- `#ENV:` lines define environment variables to be set when the image is loaded.
- An inline `## ` note after the value describes the variable, and is shown by
  `condatainer info` and in the modulefile help text.

```
#ENV:ORA_REF_PATH={prefix}/oradata   ## reference search path
```

These declarations are captured into the image's metadata at build time, so the
image is self-contained and no sidecar file travels with it. A writable `.img`
is the exception: it carries no embedded metadata, and reads a `<overlay>.env` sidecar
instead, one `KEY=value ## note` per line.

`{prefix}` is a special placeholder that is replaced **at load time** with the image's install prefix (e.g. `/cnt/orad/2.7.0`) whenever the image is loaded.

**Example:**

```bash
#ENV:CELLRANGER_REF_DIR={prefix}   ## cellranger reference dir
#ENV:GENOME_FASTA={prefix}/fasta/genome.fa   ## genome fasta
#ENV:ANNOTATION_GTF_GZ={prefix}/genes/genes.gtf.gz   ## 10X modified gtf
```

#### ENV Naming Guidelines

For common data: genome fasta, gtf, etc., use standard variable names like `GENOME_FASTA`, `ANNOTATION_GTF_GZ`.

- If the file is compressed, add `_GZ` suffix. e.g. `ANNOTATION_GTF` for uncompressed gtf, `ANNOTATION_GTF_GZ` for gzipped gtf.

For tool-specific references, use the tool name as a prefix.

- If the index is a directory, use `_DIR` suffix.
- If the index is a file prefix, use `_PREFIX` suffix.
- If the index is a specific file, use appropriate suffix based on file type.

**Examples:**

- `CELLRANGER_REF_DIR` for Cellranger references.
- `STAR_INDEX_DIR` for STAR indices.
- `BOWTIE2_PREFIX` for Bowtie2 indices.
- `BWA_MEM2_FASTA` for BWA-MEM2 genome fasta with `bwa-mem2` indices.

### Source Tag

`#SOURCE:<name> <url>` declares a file the build downloads. It can be repeated, once per file.

- **CondaTainer** downloads each source before the script runs and exposes it as `$CNT_SRC_<name>`. The name may hold letters, digits and `_`.
- A source takes one URL. `{placeholders}` from `#PH:` are substituted into it.
- The file is read-only. To change it, copy it into `$CNT_TMP` first.
- The SHA-256 of each file is recorded in the overlay and is part of its identity, so an upstream file that is re-released under the same name gives a different identity. The URL is not recorded.
- A download made inside the script, with `curl` or `wget`, is not recorded.
- A `.def` cannot declare `#SOURCE:`.

**Example:**

```bash
#!/usr/bin/env bash
#DESC:GENCODE {gencode_version} comprehensive gene annotation GTF for GRCh38
#TARGET:grch38/gtf-gencode/{gencode_version}

#PH:gencode_version:49,48,47

#SOURCE:gtf https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{gencode_version}/gencode.v{gencode_version}.primary_assembly.annotation.gtf.gz

#ENV:ANNOTATION_GTF={prefix}/gencode.v{gencode_version}.primary_assembly.annotation.gtf

cd "$CNT_PREFIX"
pigz -dc "$CNT_SRC_gtf" > "gencode.v{gencode_version}.primary_assembly.annotation.gtf"
```

#### Links only you have

`#SOURCE:<name> ask:<prompt>` asks for the link instead of writing one down. Use it for a download that expires or is tied to your account.

- **CondaTainer** shows the prompt before building, fetches the link you paste, and records only the file's hash. The link is never recorded or printed.
- A submitted job carries the answer with it, so `always_submit_data` needs nothing extra. `--yes` cannot supply a link, and an empty answer stops the build.
- If the link does not mention the overlay's version, a warning is shown. It is only a warning.
- You can use `\n` to add new lines in the prompt message.

**Example:** `cellranger/9.0.1`

```bash
#SOURCE:crx ask:10x download links expire after one day and are per-user.\nOpen https://www.10xgenomics.com/support/software/cell-ranger/downloads/previous-versions and paste the cellranger-9.0.1.tar.gz link

tar -xf "$CNT_SRC_crx" -C "$CNT_PREFIX" --strip-components=1 --no-same-owner
```

### Input Tag

`#INPUT:<Prompt>` asks for a value the script needs during the build, such as a licence acceptance or a setting.

- **CondaTainer** prompts before building and gives the answers to the script on stdin, one line each, in order: `IFS= read -r ANSWER`.
- The answer is never recorded and must not change what is built. A value that selects which overlay to build is a `#PH:` placeholder instead.
- You can use `\n` to add new lines in the prompt message.
- `#INPUT:` prompts are asked first and `#SOURCE: … ask:` prompts after them, each in the order written.
- A `.def` cannot declare `#INPUT:`.

## Apps

- Do not try to manually download apps that are already available via conda-forge or bioconda.
- Also, I don't recommend compiling apps from source unless absolutely necessary.
  - HPC systems often lack required build tools or dependencies unless you load specific modules.
  - To maximize compatibility (**CondaTainer**), it's better to rely on pre-compiled packages.

### Tips

Download with [`#SOURCE:`](#source-tag) so the file is recorded, and extract it into `$CNT_PREFIX` with `tar` and `--use-compress-program="pigz -d -p ${NCPUS:-4}"`. `pigz` uses every core the job has.

If the app requires specific environment variables to function properly, make sure to add them using `#ENV:` tags with inline `## ` notes. e.g. `orad/2.7.0`, walked through in [Custom App Recipes](../advanced_usage/custom_app.md).

### Examples

- `cellranger/9.0.1`
- `orad/2.7.0`

## Data

- Data often require downloading large files from external sources.
- Indices may need to be built using specific versions of software.
  - If indices are version dependent, ensure the app version is included in the name. e.g. `grch38/star/2.7.11b/gencode47-101`
  - If indices require building, ensure you have the scheduler parameters (e.g. `#SBATCH`) set appropriately to allocate sufficient resources.
- Always add environment variables using `#ENV:` with inline `## ` notes to help users locate and understand the reference data.

### Tips

Download with [`#SOURCE:`](#source-tag) and decompress into `$CNT_PREFIX` with `pigz -dc "$CNT_SRC_<name>" > <file>`. The overlay is compressed already, so keep the payload uncompressed instead of leaving a `.gz` inside it.

### Examples

- `grch38/genome/ucsc_no_alt`
- `grch38/transcript-gencode` (template)
- `grch38/star-gencode`
- `grcm39/salmon-gencode`

## OS

OS scripts are Apptainer definition files (`.def`) for distro-level system tools, built into a `.sqf` image. The `.def` suffix is not part of the overlay name — `ubuntu24/build-essential.def` is referred to as `ubuntu24/build-essential`.

- Use these for tools that are not available as conda packages and need a full distro environment.
- Prefer an existing upstream container image (via `Bootstrap: docker`) over building from source.

### Choosing a container root

Any `os` overlay can serve as the container root — root selection happens per
invocation, at run time, not by declaring a special type. `config base` names
the one CondaTainer builds by default and falls back to when nothing else in
the requested overlays is root-eligible; a project's own recipe collection may
still name it `<distro>/base` for readability, but that name carries no
special meaning to the type system.

Before packing any `.def` build, CondaTainer checks that its finished sandbox
can serve a build that might later run inside it, and refuses one that cannot:

| Tool | From | |
| --- | --- | --- |
| `/bin/bash` | bash | required, **at that exact path** — every script CondaTainer runs inside a chosen root is launched as `/bin/bash`, a recipe included |
| `apptainer` | — | optional — needed only to mount an image inside the container; a sandbox without it warns |

Packing (`mksquashfs`) and conda installs (`micromamba`) run through
CondaTainer's own self-provisioned toolchain, not whatever the root happens to
carry, so neither is checked here. Tools CondaTainer runs on the host instead —
`unsquashfs`, `debugfs`, `e2fsck`, `resize2fs`, `mke2fs` — are host
prerequisites, not a root's concern either.

### OS Templates

`.def` files support [Template Tags](#template-tags) too. Placeholders are substituted **throughout the whole file** at build time — including the `Bootstrap`/`From` header — so one file can cover every upstream image tag.

`ubuntu24/posit-r.def` builds every R version from a single definition:

```
#PH:version:4.4.3,4.5.0,4.5.1
#AUTOUPDATE:version:docker:posit/r-base:^(\d+\.\d+\.\d+)-noble(?:-[^-]+)?$>=3.1.3

#TARGET:ubuntu24/r{version}
#DESC: Ubuntu noble R {version}

Bootstrap: docker
From: posit/r-base:{version}-noble
```

- `#TARGET:ubuntu24/r{version}` expands to `ubuntu24/r4.4.3`, `ubuntu24/r4.5.0`, …
- `{version}` in `From:` is substituted at build time, so each expansion pulls its own upstream tag.
- `#AUTOUPDATE:` can track a Docker tag to keep the `#PH:` list current. See [Auto-Update Tag](#auto-update-tag).

### Examples

- `ubuntu24/code-server.def`
- `ubuntu24/posit-r.def` (template)
- `ubuntu24/xfce4.def`
