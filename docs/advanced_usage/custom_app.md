# Custom App Recipes

Most tools need no recipe at all. If the app is on a Conda channel, `condatainer create <name>/<version>` builds a [module overlay](../user_guide/concepts.md#overlay-types) straight from Conda:

```bash
condatainer create samtools/1.22.1   # resolved from bioconda, no recipe needed
condatainer create openjdk/17.0.18   # resolved from conda-forge
```

You need a recipe when the app is **not packaged for Conda** (vendor tarball, installer), or when you need a version Conda doesn't carry.

This page walks through two:

- **[orad/2.7.0](#example-1-orad--a-single-version-recipe)** — a plain single-version recipe.
- **[cytoscape](#example-2-cytoscape--a-version-template)** — a template covering every release from one file, used by the [`cytoscape` helper](./custom_helper.md#example-2-cytoscape-an-app-not-on-conda).

See the [Build Script Manual](../manuals/build_script.md) for the complete header and variable reference.

## Where Recipes Live

Recipes come from **sources**: an ordered list of recipe collections in your config, where the first match wins. A source is a directory (or URL) holding a `recipes/` folder:

```
/scratch/me/recipes/
└── recipes/
    ├── orad/2.7.0
    └── cytoscape
```

Add yours ahead of the others:

```bash
condatainer config prepend sources mine=/scratch/me/recipes
condatainer avail cytoscape   # confirm CondaTainer sees the recipe
```

The **path under `recipes/` is the module name**, with no `.sh` suffix:

- `recipes/cellranger/9.0.1` => `cellranger/9.0.1`
- `recipes/cytoscape` => `cytoscape`

See [Configuration](../manuals/configuration.md) for how sources are merged across config files.

## What a Recipe Is

A recipe is a shell script with a header block of `#KEY:` lines. The body runs top to bottom as `bash -euo pipefail`. There is no wrapper function and no injected helper: call `curl`, `tar` and `pigz` directly, and print progress with `echo ... >&2`.

| Variable | Meaning |
|---|---|
| `$CNT_PREFIX` | where the payload goes; exactly this directory is packed into the overlay |
| `$CNT_NAME` | the complete module name, e.g. `orad/2.7.0` |
| `$CNT_TMP` | scratch directory, also exported as `$TMPDIR` |
| `$CNT_SRC_<name>` | a file fetched by a [`#SOURCE:`](../manuals/build_script.md#source-tag) header, read-only |
| `$NCPUS`, `$MEM`, `$MEM_GB` | resources, taken from your scheduler directives |

## Which Shape to Write

| | Single-version | Template |
|---|---|---|
| **File path** | `<name>/<version>` | `<name>` |
| **Headers** | none required | `#PH:` + `#TARGET:` |
| **Write it when** | install differs between versions,<br/> or you only need one | the logic is identical |
| **Example** | [`orad/2.7.0`](#example-1-orad--a-single-version-recipe) | [`cytoscape`](#example-2-cytoscape--a-version-template) |

Start single-version. Promote to a template once you've copied the same recipe to a second version and changed nothing but the URL.

## Example 1: orad — A Single-Version Recipe

Illumina's ORA decompressor is a vendor tarball behind a download page (never going to be on Conda). The file lives at `recipes/orad/2.7.0`, so it installs as `orad/2.7.0`.

```bash
#!/usr/bin/env bash
#DESC:Illumina ORA Decompressor
#URL:https://support.illumina.com/sequencing/sequencing_software/DRAGENORA/software-downloads.html

#SOURCE:orad https://s3.amazonaws.com/webdata.illumina.com/downloads/software/dragen-decompression/orad.2.7.0.linux.tar.gz

#ENV:ORA_REF_PATH={prefix}/oradata   ## Illumina ORA decompressor reference search path

echo "Extracting $CNT_NAME" >&2
tar -xf "$CNT_SRC_orad" -C "$CNT_PREFIX" --strip-components=1 \
    --use-compress-program="pigz -d -p ${NCPUS:-4}"

# orad ships its binary at the archive root, so move it to bin/
mkdir -p "$CNT_PREFIX/bin"
mv "$CNT_PREFIX/orad" "$CNT_PREFIX/bin"
```

`#DESC:` is what users see in `condatainer avail` and `condatainer info` — always write it.

Four things to take from it:

1. **Everything installs under `$CNT_PREFIX`.** CondaTainer packs exactly that directory into the `.sqf` overlay — anything written elsewhere is lost.
2. **`$CNT_PREFIX/bin` goes on `$PATH`** when the overlay is loaded.
3. **`#SOURCE:` downloads the file** before the recipe runs and records its SHA-256 in the overlay, so an upstream file that is re-released under the same name gives a different identity. It is read-only; copy it into `$CNT_TMP` if you need to change it.
4. The recipe never says where the download came from or how it got there. It reads `$CNT_SRC_orad`.

````{note}
The extraction command depends on the tarball structure. Most tarballs wrap everything in a single top-level `<name>/` (or `<name>-<version>/`) directory, so extracting straight into `$CNT_PREFIX` with `--strip-components=1` is enough.

Layouts vary, though. Check before you write the line:

```bash
tar -tzf <file>.tar.gz | head
```

Then fix up whatever the archive got wrong — `orad` ships its binary at the root rather than in `bin/`, hence the `mv` above.
````

### Setting Variables the App Needs

**Most apps need nothing here.** Once the binary is on `$PATH`, they run. Skip this section unless the app fails without a variable set.

`orad` is one that does: it looks up its decompression reference through `$ORA_REF_PATH` and errors out if no reference is found. Rather than making every user export it by hand, the recipe declares it:

```bash
#ENV:ORA_REF_PATH={prefix}/oradata   ## Illumina ORA decompressor reference search path
```

- `{prefix}` is replaced at load time with the image's install prefix (`/cnt/orad/2.7.0`) whenever the image is loaded.
- The `## ` note after the value becomes the description shown by `condatainer info`.

### Build and Verify

```bash
condatainer create orad/2.7.0
```

Load the overlay and check the binary and the variable:

```bash
condatainer exec -o orad/2.7.0 \
    bash -c 'command -v orad && echo "$ORA_REF_PATH"'
```

`info` shows the same thing without entering a container:

```bash
condatainer info orad/2.7.0
# Environment
#  - ORA_REF_PATH=/cnt/orad/2.7.0/oradata
#    # Illumina ORA decompressor reference search path
```

## Example 2: Cytoscape — A Version Template

Cytoscape publishes a pre-built Linux tarball for every release at a predictable URL. Writing one file per version would mean a dozen near-identical recipes, so this is a **template**: a single file at `recipes/cytoscape` that expands into `cytoscape/3.10.4`, `cytoscape/3.9.1`, and so on.

```bash
#!/usr/bin/env bash
#DESC:Cytoscape {cytoscape_version} — network biology visualization platform (needs Java)
#URL:https://github.com/cytoscape/cytoscape/releases
#TARGET:cytoscape/{cytoscape_version}

#PH:cytoscape_version:3.10.4,3.10.3,3.10.2,3.10.1,3.10.0,3.9.1,3.9.0

#SOURCE:cytoscape https://github.com/cytoscape/cytoscape/releases/download/{cytoscape_version}/cytoscape-unix-{cytoscape_version}.tar.gz

tar -xf "$CNT_SRC_cytoscape" -C "$CNT_PREFIX" --strip-components=1 \
    --use-compress-program="pigz -d -p ${NCPUS:-4}"

mkdir -p "$CNT_PREFIX/bin"
ln -s ../cytoscape.sh "$CNT_PREFIX/bin/cytoscape"
```

| Header | Role |
|---|---|
| `#PH:` | Declares the `cytoscape_version` placeholder and its allowed values |
| `#TARGET:` | Module name pattern — expands to `cytoscape/3.10.4`, `cytoscape/3.9.1`, … |
| `#SOURCE:` | The download; its URL takes the placeholder too |
| `#DESC:` | Shown in `condatainer avail` and `condatainer info` |

Every `#PH:` name must appear as a `{name}` token in `#TARGET:` and vice versa — otherwise every value would collapse onto the same target.

`{cytoscape_version}` tokens are substituted **before** the recipe runs — in the headers *and* in the body. It is not a bash variable, so make sure there is no `${cytoscape_version}` in the script.

### Installing From a Template

Users can name the concrete target directly, or let CondaTainer prompt:

```bash
condatainer create cytoscape          # prompts for each placeholder
condatainer create cytoscape/3.10.4   # direct — fills cytoscape_version=3.10.4
```

`condatainer avail` shows the template collapsed, with its placeholders listed:

```
cytoscape  [7 variants]
  Cytoscape {cytoscape_version} — network biology visualization platform (needs Java)
  → cytoscape/{cytoscape_version}
  - cytoscape_version:  3.9.0-3.10.4  (7 values)
```

```{tip}
At the placeholder prompt, hit <kbd>Tab</kbd> to see the available values and autocomplete.
```

A template can declare several placeholders, and they aren't limited to versions. See [Template Tags](../manuals/build_script.md#template-tags) for the full syntax, and [Custom Data Build Scripts](./custom_data.md#example-2-derived-data-as-a-template) for a template with three placeholders.

### Build and Verify

```bash
condatainer create cytoscape/3.10.4
```

Cytoscape needs a JVM 17, which can come from Conda — no recipe required:

```bash
condatainer create openjdk/17.0.18
```

Load both overlays and check the binary resolves:

```bash
condatainer exec -o openjdk/17.0.18 -o cytoscape/3.10.4 \
    bash -c 'command -v cytoscape && java -version'
```

```{note}
Overlay order is the stack order — later `-o` flags sit on top. Keep the app last so its `bin/` wins on conflicts.
```

## Links Only You Have

Some vendor downloads have no stable URL — Cell Ranger's link is signed, tied to your account, and expires after a day. Declare it with the `ask:` form of `#SOURCE:` and CondaTainer prompts for it before the build starts:

```bash
#!/usr/bin/env bash
#DESC:10x Genomics Cell Ranger, single-cell gene expression analysis
#URL:https://www.10xgenomics.com/support/software/cell-ranger/downloads/previous-versions

#SOURCE:crx ask:10x download links expire after one day and are per-user.\nOpen https://www.10xgenomics.com/support/software/cell-ranger/downloads/previous-versions and paste the cellranger-9.0.1.tar.gz link

tar -xf "$CNT_SRC_crx" -C "$CNT_PREFIX" --strip-components=1 --no-same-owner \
    --use-compress-program="pigz -d -p ${NCPUS:-4}"
```

CondaTainer prints the prompt (`\n` splits it across lines), fetches the link you paste, and gives the recipe the file as `$CNT_SRC_crx`. The link itself is never recorded or printed, and a warning appears if it doesn't mention the version. See [Source Tag](../manuals/build_script.md#source-tag) for the rules.

A value that isn't a download — a licence acceptance, say — uses [`#INPUT:`](../manuals/build_script.md#input-tag) instead, which the recipe reads on stdin.

Expected output:

```
$ condatainer create cellranger/9.0.1
[CNT◇] 10x download links expire after one day and are per-user.
[CNT◇] Open https://www.10xgenomics.com/support/software/cell-ranger/downloads/previous-versions and paste the cellranger-9.0.1.tar.gz link
Enter here: 
```

## Share It

Put the recipe in a collection other people's config lists as a source. See [Sharing Your Recipes](../deployment/share_scripts.md).

## Related

- [Build Script Manual](../manuals/build_script.md) — full header/variable reference
- [Sharing Your Recipes](../deployment/share_scripts.md) — hosting your own collection
- [Custom Data Build Scripts](./custom_data.md) — packaging datasets, databases, and derived artifacts
- [Custom Helper Scripts](./custom_helper.md) — launching the app as a browser service
- [Custom OS Overlays](./custom_os.md) — when you need system packages (`apt`) instead
- [Custom Bundle Overlays](./custom_bundle.md) — read-only Conda environments
