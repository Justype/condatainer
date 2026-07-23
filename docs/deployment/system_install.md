# System-Wide Installation

For cluster and server maintainers: install CondaTainer once in a read-only location, so every user on the machine gets the binary, your curated build scripts, and any overlays you pre-build.

```{note}
Setting up for yourself or one lab instead?

- an individual user → [Installation](../user_guide/installation.md)
- a lab or group → [Shared Group Installation](./group_install.md)
```

A system install serves everyone on the machine, and leaves groups and individuals free to add their own on top. Users never write to your directory: CondaTainer searches every configured directory but writes to the first *writable* one, which for them is their own scratch or home. See [Data Layers](./data_layers.md).

## 1. Install and Initialize

Both paths lay down a read-only tree that every user reads but none can write. Pick the one that matches how your users get the binary: a versioned **module** tree, with the data pinned by `$CNT_ROOT`, or a plain **directory** the binary auto-detects as its app-root.

### With environment modules

The usual choice on a cluster: users opt in with `module load condatainer`, and the module file wires the environment (`$PATH`, `$CNT_ROOT`, Apptainer).

**The directory layout.** Module trees are versioned, so the binary lands at a path like `/opt/condatainer/1.4.0/bin/condatainer`. If that path defined the app-root, your build scripts and overlays would be tied to one release and stranded at the next upgrade. So keep the data in a version-independent directory and point `$CNT_ROOT` at it — then only the binary moves between versions:

```
/opt/condatainer/
├── 1.4.0/bin/condatainer     # each release in its own dir
├── 1.5.0/bin/condatainer
└── data/                     # CNT_ROOT — shared across all versions
    ├── build-scripts/
    ├── helper-scripts/
    └── images/
```

**Install a version.** Download the binary and lock the tree read-only:

```bash
mkdir -p /opt/condatainer/1.4.0/bin
curl -fsSL -o /opt/condatainer/1.4.0/bin/condatainer \
    https://github.com/Justype/condatainer/releases/download/v1.4.0/condatainer_linux_x86_64
chmod 755 /opt/condatainer/1.4.0/bin/condatainer
chmod -R go-w /opt/condatainer/data
```

**Write the module file.** It puts the binary on `$PATH`, points `$CNT_ROOT` at the data dir, and — because **CondaTainer needs Apptainer** — makes Apptainer available. Either pin an exact binary with `$CNT_APPTAINER_BIN`, or have the module load Apptainer's own module.

`/opt/modulefiles/condatainer/1.4.0.lua` (Lmod):

```lua
whatis("Conda environments, tools, and data as single-file overlays")

prepend_path("PATH", "/opt/condatainer/1.4.0/bin")
setenv("CNT_ROOT", "/opt/condatainer/data")

-- Apptainer: pin an exact binary...
setenv("CNT_APPTAINER_BIN", "/opt/apptainer/1.3.0/bin/apptainer")
-- ...or load its module instead (drop the setenv above):
-- depends_on("apptainer/1.3.0")
```

`/opt/modulefiles/condatainer/1.4.0.` (environment-modules/Tcl):

```tcl
#%Module1.0
prepend-path PATH /opt/condatainer/1.4.0/bin
setenv CNT_ROOT /opt/condatainer/data

# Apptainer: pin an exact binary...
setenv CNT_APPTAINER_BIN /opt/apptainer/1.3.0/bin/apptainer
# ...or load its module instead (drop the setenv above):
# depends-on apptainer/1.3.0     ;# environment-modules ≥ 4.4
# module load apptainer/1.3.0    ;# older environment-modules
```

**Site defaults.** Load the module you just wrote — it puts `condatainer` on `$PATH` and sets `$CNT_ROOT` — then set tunable defaults like `build.ncpus` or `channels` in the app-root `config.yaml`:

```bash
module load condatainer
condatainer config set build.ncpus 8 -l app-root
```

```{note}
`$CNT_ROOT` is shared across all versions, so this `config.yaml` applies to every installed version — a `config set` from any loaded version changes them all.
```

### Without modules

The binary lives in its own directory with the data beside it, so the app-root is auto-detected. A binary at `<cnt_root>/bin/condatainer` makes `<cnt_root>` the app-root, and `<cnt_root>/{build-scripts,helper-scripts,images}` are searched for every user:

```
/opt/condatainer/
├── bin/condatainer
├── config.yaml
├── build-scripts/
├── helper-scripts/
└── images/
```

Install, initialize, then lock it:

```bash
mkdir -p /opt/condatainer/bin
curl -fsSL -o /opt/condatainer/bin/condatainer \
    https://github.com/Justype/condatainer/releases/latest/download/condatainer_linux_x86_64
chmod 755 /opt/condatainer/bin/condatainer

/opt/condatainer/bin/condatainer config init -l app-root
chmod -R go-w /opt/condatainer
```

`config init` detects Apptainer and the scheduler and records them in the shared `config.yaml`, so no user has to. Put it on everyone's `$PATH` with a profile snippet in `/etc/profile.d/condatainer.sh`, or a symlink — symlinks are resolved first, so app-root detection still points back to the real install:

```bash
ln -s /opt/condatainer/bin/condatainer /usr/local/bin/condatainer
# Or add /opt/condatainer/bin to PATH
```

## 2. What to Put There

Anything you want every user to have without rebuilding:

```
<cnt_root>/           # beside the binary, or wherever CNT_ROOT points
├── config.yaml       # site defaults: apptainer path, scheduler, build limits
├── build-scripts/    # curated or site-specific recipes
├── helper-scripts/   # site-specific services
└── images/           # pre-built overlays: common tools, reference data
```

Pre-building `images/` is the highest-value part on a shared cluster: a genome index or a large tool overlay built once is one file that every user reads, instead of one copy per user.

```{important}
Re-apply `chmod -R go-w <cnt_root>` after adding content. Subdirectories keep the permissions they were created with, and a world-writable `images/` becomes the write target for every user — their builds would land in your tree instead of their own.
```

For scripts, you can also point at a remote source instead of copying files in — see [Sharing Your Scripts](./share_scripts.md).

## 3. Site Defaults in `config.yaml`

Settings here apply to everyone and can still be overridden per user:

```yaml
apptainer_bin: /usr/bin/apptainer
build:
  ncpus: 8
  mem: 32g
```

Scalars: the user's own config wins. Directory and source lists: merged across layers, so your entries are always searched even when a user has their own config.

```{note}
`/etc/condatainer/config.yaml` is also read, at the lowest priority. Use it when the binary is packaged system-wide and you have no app-root; otherwise keep site settings next to the install.
```

## 4. Verify as an Unprivileged User

```bash
su - someuser
module load condatainer     # if using modules

condatainer config paths    # your dirs listed; writable target is the user's own
condatainer avail           # your build-scripts appear
condatainer list            # your pre-built overlays appear
```

`config paths` marks which directory receives writes. On a correct setup that is the user's scratch or `~/.local/share/condatainer`, never yours.

## Group Use Under System Install

A system install is read-only, but it doesn't lock anyone out — groups and individuals layer their own directories on top of it. A lab that wants to share data among its members, without touching system install, points CondaTainer at a group directory with `$CNT_EXTRA_ROOT`:

```bash
export CNT_EXTRA_ROOT=/shared/labA/condatainer
```

That one variable adds the group's `build-scripts/`, `helper-scripts/`, and `images/` to the search path for everyone who sets it — most conveniently from a group `setup.sh` or an environment module. It's env-only, so put it in a module file or a shell profile.

The group root sits *above* the system install in the search order, so a group script can shadow a same-named one you ship, and writes go to the group directory. See [Shared Group Installation](./group_install.md) for setting up the group directory itself, and [Sharing Your Scripts](./share_scripts.md) for publishing into it.

## Related

- [Data Layers](./data_layers.md) — how the search order and write target work
- [Shared Group Installation](./group_install.md) — the lab-level setup
- [Sharing Your Scripts](./share_scripts.md) — publishing scripts to your users
- [Configuration Manual](../manuals/configuration.md#multi-tier-setup-system--group--user) — config layers and search paths
- [Installation](../user_guide/installation.md) — the standard single-user install
