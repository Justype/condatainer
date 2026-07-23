# Custom App Build Scripts

Most tools need no build script at all. If the app is on a Conda channel, `condatainer create <name>/<version>` builds a [module overlay](../user_guide/concepts.md#-overlay-types) straight from Conda:

```bash
condatainer create samtools/1.22.1   # resolved from bioconda, no script needed
condatainer create openjdk/17.0.18   # resolved from conda-forge
```

You need a build script when the app is **not packaged for Conda** (vendor tarball, installer), or when you need a version Conda doesn't carry.

This page walks through two real ones:

- **[orad/2.7.0](#example-1-orad--a-single-version-script)** — a plain single-version script.
- **[cytoscape](#example-2-cytoscape--a-version-template)** — a template covering every release from one file, used by the [`cytoscape` helper](./custom_helper.md#example-2-cytoscape-an-app-not-on-conda).

See the [Build Script Manual](../manuals/build_script.md) for the complete header and variable reference.

## Where Build Scripts Live

Build scripts are looked up in your `build-scripts/` search paths:

```bash
condatainer config paths      # show build-scripts/ helper-scripts/ search dirs
condatainer avail cytoscape   # confirm CondaTainer sees the script
```

Drop your script into the first writable `build-scripts/` directory (no `.sh` suffix). The **relative file path is the module name**: 

- `build-scripts/cellranger/9.0.1` => `cellranger/9.0.1`
- `build-scripts/cytoscape` => `cytoscape`

## Which Shape to Write

| | Single-version | Template |
|---|---|---|
| **File path** | `<name>/<version>` | `<name>` |
| **Headers** | none required | `#PL:` + `#TARGET:` |
| **Write it when** | install differs between versions,<br/> or you only need one | the logic is identical |
| **Example** | [`orad/2.7.0`](#example-1-orad--a-single-version-script) | [`cytoscape`](#example-2-cytoscape--a-version-template) |

Start single-version. Promote to a template once you've copied the same script to a second version and changed nothing but the URL.

Either way, copy [`assets/build-template-apps`](https://github.com/Justype/cnt-scripts/blob/main/assets/build-template-apps), which carries the [shared boilerplate](#the-shared-boilerplate) plus skeletons for the two common formats (tarball and single binary).

## Example 1: orad — A Single-Version Script

Source: [`build-scripts/orad/2.7.0`](https://github.com/Justype/cnt-scripts/blob/main/build-scripts/orad/2.7.0)

Illumina's ORA decompressor is a vendor tarball behind a download page (never going to be on Conda). The file lives at `build-scripts/orad/2.7.0`, so it installs as `orad/2.7.0`.

### Headers

```bash
#!/usr/bin/bash
#WHATIS:Illumina ORA Decompressor
#URL:https://support.illumina.com/sequencing/sequencing_software/DRAGENORA/software-downloads.html

#ENV:ORA_REF_PATH=$app_root/oradata
#ENVNOTE:Illumina ORA decompressor reference search path
```

`#WHATIS:` is what users see in `condatainer avail` and `condatainer info` — always write them.

The `#ENV:` pair is app-specific, covered [below](#setting-variables-the-app-needs).

### `install_app()`

```bash
install_app() {
    # By default, $target_dir and $tmp_dir are created and we are now in $tmp_dir
    os=$(uname -s)
    case $os in
        Linux)  os="linux";;
        Darwin) os="mac";;
    esac

    print_stderr "Downloading ${YELLOW}${app_name_version}${NC} binaries"
    url="https://s3.amazonaws.com/webdata.illumina.com/downloads/software/dragen-decompression/orad.${version}.${os}.tar.gz"
    curl -fsSL -o "${app_name}_${version}.tar.gz" "$url"

    print_stderr "Extracting ${YELLOW}${app_name_version}${NC}"
    tar_xf_pigz "${app_name}_${version}.tar.gz" -C "$target_dir" --strip-components=1

    # orad is under module folder, move it to bin
    mkdir -p "$target_dir/bin" # modules/apps/name/version/bin
    mv "$target_dir/orad" "$target_dir/bin"
}
```

`install_app()` is the only part you write. Four things to take from it:

1. **Everything installs under `$target_dir`.** CondaTainer packs exactly that directory into the `.sqf` overlay — anything written elsewhere is lost.
2. **`$target_dir/bin` goes on `$PATH`** when the overlay is loaded.
3. `tar_xf_pigz` is a helper function: use `pigz` when available.

````{note}
The extraction command depends on the tarball structure. Most tarballs wrap everything in a single top-level `<name>/` (or `<name>-<version>/`) directory, so extracting straight into `$target_dir` with `--strip-components=1` is enough.

Layouts vary, though. Check before you write the line:

```bash
tar -tzf <file>.tar.gz | head
```

Then fix up whatever the archive got wrong — `orad` ships its binary at the root rather than in `bin/`, hence the `mv` above.
````

### Setting Variables the App Needs

**Most apps need nothing here.** Once the binary is on `$PATH`, they run. Skip this section unless the app fails without a variable set.

`orad` is one that does: it looks up its decompression reference through `$ORA_REF_PATH` and errors out if no reference found. Rather than making every user export it by hand, the script declares it:

```bash
#ENV:ORA_REF_PATH=$app_root/oradata
#ENVNOTE:Illumina ORA decompressor reference search path
```

- `$app_root` is replaced at install time with the overlay's path (`/cnt/orad/2.7.0`). The variable is then set whenever the overlay is loaded.
- `#ENVNOTE:` must directly follow its `#ENV:` line; it becomes the description shown by `condatainer info`.

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

Source: [`build-scripts/cytoscape`](https://github.com/Justype/cnt-scripts/blob/main/build-scripts/cytoscape)

Cytoscape publishes a pre-built Linux tarball for every release at a predictable URL. Writing one file per version would mean a dozen near-identical scripts, so this is a **template**: a single file at `build-scripts/cytoscape` that expands into `cytoscape/3.10.4`, `cytoscape/3.9.1`, and so on.

### Headers

```bash
#!/usr/bin/bash
#PL:cytoscape_version:3.9.0,3.9.1,3.10.0,3.10.1,3.10.2,3.10.3,3.10.4
#AUTOUPDATE:cytoscape_version:github:cytoscape/cytoscape>=3.9.0

#TARGET:cytoscape/{cytoscape_version}
#WHATIS:Cytoscape {cytoscape_version} — network biology visualization platform (needs Java)
#URL:https://github.com/cytoscape/cytoscape/releases
```

| Header | Role |
|---|---|
| `#PL:` | Declares the `cytoscape_version` placeholder and its allowed values |
| `#TARGET:` | Module name pattern — expands to `cytoscape/3.10.4`, `cytoscape/3.9.1`, … |
| `#AUTOUPDATE:` | Lets CI refresh the `#PL:` list from GitHub releases (`>=3.9.0` is the floor) |
| `#WHATIS:` | Shown in `condatainer avail` and `condatainer info` |

Every `#PL:` name must appear as a `{name}` token in `#TARGET:` and vice versa — otherwise every value would collapse onto the same target, so CondaTainer warns and skips the expansion.

`{cytoscape_version}` tokens are substituted **before** the script runs — in the headers *and* in the body. (It is not a bash variable, make sure no `${cytoscape_version}` in the script)

### `install_app()` in template

```bash
install_app() {
    local url="https://github.com/cytoscape/cytoscape/releases/download/{cytoscape_version}/cytoscape-unix-{cytoscape_version}.tar.gz"

    print_stderr "Downloading ${YELLOW}${app_name_version}${NC}"
    curl -fsSL -o "cytoscape.tar.gz" "$url"

    print_stderr "Extracting ${YELLOW}${app_name_version}${NC}"
    tar_xf_pigz "cytoscape.tar.gz" -C "$target_dir" --strip-components=1

    mkdir -p "$target_dir/bin"
    ln -s "$target_dir/cytoscape.sh" "$target_dir/bin/cytoscape"
}
```

Structurally identical to orad — download, extract to `$target_dir`. The only template-specific part is `{cytoscape_version}` in the URL, substituted textually before executed.

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

Cytoscape needs a JVM 17, which can come from Conda — no script required:

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

## Prompting for Input

Some vendor downloads have no stable URL — Cell Ranger's link is signed and expires after a day. When there's nothing reliable to hardcode, an `#INTERACTIVE:` header asks the user for the value before the build starts (see [`cellranger/9.0.1`](https://github.com/Justype/cnt-scripts/blob/main/build-scripts/cellranger/9.0.1)):

```bash
#INTERACTIVE:⚠️ 10X links only valid for one day. Please go to the link below and get tar.gz link.\nhttps://www.10xgenomics.com/support/software/cell-ranger/downloads/previous-versions
```

CondaTainer prints the prompt (`\n` splits it across lines) and feeds the user's answer to the script on **stdin**, so read it with `read -r`. Declare one `#INTERACTIVE:` per value, read in order:

```bash
install_app() {
    url='https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.1.tar.gz'
    # A live 10X link carries an Expires= token; a bare one has gone stale.
    if [[ "$url" != *"Expires="* ]]; then
        print_stderr "Enter the download link: "
        read -r url
    fi
    # Validate the pasted URL before trusting it.
    if [[ "$url" != *"cellranger-$version.tar.gz"* ]]; then
        print_stderr "❌ The download link is invalid or not for version $version."
        return 1
    fi
    curl -fsSL -o "${app_name}_${version}.tar.gz" "$url"
    tar_xf_pigz "${app_name}_${version}.tar.gz" -C "$target_dir" --strip-components=1
}
```

Expected output:

```
$ condatainer i cellranger/9.0.1
[CNT◇] Installing to /path/condatainer/images (app-root)
[CNT◇] ⚠️ 10X links only valid for one day. Please go to the link below and get tar.gz link.
[CNT◇] https://www.10xgenomics.com/support/software/cell-ranger/downloads/previous-versions
Enter here: 
```

## The Shared Boilerplate

Everything below `install_app()` in both scripts (`print_stderr`, the `install()` driver, the `pigz_*` helpers, `set -e`, the final `install` call) is a block copied verbatim between scripts — it validates `$target_dir`/`$tmp_dir`, refuses to overwrite an existing target, then `cd`s into `$tmp_dir` and calls your `install_app()`.

Don't write it from scratch:

- [`assets/build-template-apps`](https://github.com/Justype/cnt-scripts/blob/main/assets/build-template-apps) — apps
- [`assets/build-template-ref`](https://github.com/Justype/cnt-scripts/blob/main/assets/build-template-ref) — data / reference (see [Custom Data Build Scripts](./custom_data.md))

## Share It

Open a PR against [cnt-scripts](https://github.com/Justype/cnt-scripts) so everyone gets the script via `condatainer avail`, or host your own repo for internal tools. See [Sharing Your Scripts](../deployment/share_scripts.md).

## Related

- [Build Script Manual](../manuals/build_script.md) — full header/variable reference
- [Sharing Your Scripts](../deployment/share_scripts.md) — upstreaming or hosting your own source
- [Custom Data Build Scripts](./custom_data.md) — packaging datasets, databases, and derived artifacts
- [Custom Helper Scripts](./custom_helper.md) — launching the app as a browser service
- [Custom OS Overlays](./custom_os.md) — when you need system packages (`apt`) instead
- [Custom Bundle Overlays](./custom_bundle.md) — read-only Conda environments
