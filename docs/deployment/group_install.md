# Shared Group Installation

To avoid duplicate installation for a group, CondaTainer can be installed into a shared directory. It detects its own install root from the binary's location, so everyone can read it and run it.

```{note}
Setting up for yourself or a whole cluster instead?

- an individual user → [Installation](../user_guide/installation.md)
- a cluster or server maintainer → [System-Wide Installation](./system_install.md)
```

## 1. Pick a Directory

The path must be:

- **owned by the group and group-writable** — members build overlays into it themselves;
- **visible from the compute nodes** — overlays are read from here at job runtime.

```bash
cd /shared/labA
mkdir condatainer
chgrp labA condatainer      # your group name
chmod 2775 condatainer      # group-writable + setgid
cd condatainer              # the rest of this guide runs from here
```

The **setgid** bit (the `2`) makes everything created inside inherit the `labA` group.

## 2. Install Into It

Download the binary into `bin/`:

```bash
mkdir -m 2775 bin
curl -fsSL -o bin/condatainer \
    https://github.com/Justype/condatainer/releases/latest/download/condatainer_linux_x86_64
chmod 775 bin/condatainer
```

## 3. Initialize the Shared Config

Write the config at the app-root layer so it applies to everyone using this install:

```bash
./bin/condatainer config init -l app-root
```

This detects Apptainer and the scheduler once, for everyone.

The directory is group-writable, so a member running `condatainer config set` would change this file for the whole group without meaning to. Drop the group write bit:

```bash
chmod 644 config.yaml
```

Their changes then land in their own config instead. See [Configuration Manual](../manuals/configuration.md#config-file-locations).

## 4. Give Everyone Access

Write one shell config at the group root, so settings can change later without anyone editing their rc file.

`/shared/labA/condatainer/setup.sh`:

```bash
export PATH="/shared/labA/condatainer/bin:$PATH"

source <(condatainer completion)   # shell is auto-detected

# umask 002   # uncomment if mostly work under the group dir
              # 002 makes hand-created files group-writable too
```

`umask` line is not requried: CondaTainer checks whether the parent directory is group-writable and adds the group-write bit to files/dirs it creates, regardless of umask.

Each member adds one line to their `~/.bashrc` or `~/.zshrc`:

```bash
source /shared/labA/condatainer/setup.sh
```

Verify from another member's account:

```bash
condatainer config paths   # shared dirs appear as app-root entries
condatainer avail
```

## What Gets Shared

The shared root is in everyone's search path, so:

- **Overlays** in `images/` : apps and data overlays. `condatainer exec/run`
- **Build scripts** in `build-scripts/` appear in everyone's `condatainer avail`.
- **Helper scripts** in `helper-scripts/` are launchable by everyone. `condatainer helper`

Writes go to the first *writable* directory in the search order ([Data Layers](./data_layers.md)). If a member's builds land in their own scratch instead of the shared root, check that the shared directory is still group-writable with the setgid bit set ([step 1](#1-pick-a-directory)).

## Related

- [Data Layers](./data_layers.md) — how the search order and write target work
- [Sharing Your Scripts](./share_scripts.md) — publishing build and helper scripts to the group
- [Installation](../user_guide/installation.md) — the standard single-user install
- [Configuration Manual](../manuals/configuration.md#multi-tier-setup-system--group--user) — config layers and search paths
