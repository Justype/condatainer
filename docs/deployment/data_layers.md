# Data Layers

CondaTainer keeps data under a small set of directories and searches them together. This page explains what those directories are, why the design is shaped this way, and how to work with more than one at once. For step-by-step setup, follow the guide for your situation:

| You are | Guide |
|---|---|
| an individual user | [Installation](../user_guide/installation.md) |
| setting up a lab or group | [Shared Group Installation](./group_install.md) |
| a cluster or server maintainer | [System-Wide Installation](./system_install.md) |

## app-root

`app-root` is an install's own home — the directory its `images/`, `build-scripts/`, and `helper-scripts/` live under:

```
<app-root>/
├── bin/
│   └── condatainer  # the binary; its location can define app-root
├── build-scripts/   # build recipes
├── helper-scripts/
└── images/          # overlays
```

You decide where it is: CondaTainer takes it from the binary's location (a binary at `<app-root>/bin/condatainer` makes `<app-root>` the app-root) or from `$CNT_ROOT`, which overrides that. For most installs, choosing where `app-root` points is the only decision that matters; the exact resolution rules are in [How app-root Is Found](#how-app-root-is-found).

## What It's Designed For

Packing environments into images is about reducing inodes, and the largest saving is not rebuilding what someone already built. CondaTainer is shaped around three ways people share — or deliberately don't:

- **A single user** wants one place to keep everything. Point `app-root` at a directory you own; every build and every read happens there.
- **A group** wants a shared place to read *and* write recipes and overlays together. Point `app-root` at a directory the group can write, and the whole lab contributes to and reads from one pool.
- **A cluster or system**, often alongside environment modules, wants a site-wide set of common overlays (RStudio Server, standard tools) built once by an admin and read by everyone. Here `app-root` is admin-owned and **read-only**, so no one's builds can land in it. Each consumer keeps their own writable place instead:
  - a **single user** writes to their personal space (`$SCRATCH`, or home) with no setup at all;
  - a **group** points `$CNT_EXTRA_ROOT` at a shared directory of their own, and their builds pool there.

The same binary and the same commands serve all three. What changes is only where `app-root` points, who may write it, and whether a writable layer sits beside it — the search rule below is what ties them together.

## The Layers

Highest priority first. Not all of them exist in every deployment — each is optional and appears only when something sets it, so a given install usually sees `app-root` and `user`.

| Layer | Set by | Typical path |
|---|---|---|
| extra-root | `$CNT_EXTRA_ROOT` env var | `/shared/labA/condatainer` |
| app-root | binary location, or `$CNT_ROOT` | `/anywhere/condatainer` |
| scratch | `$SCRATCH` | `$SCRATCH/condatainer` |
| user | XDG variables | `~/.local/share/condatainer` |

The stack is identical for everyone; only the layer that receives writes moves between the three cases from [What It's Designed For](#what-its-designed-for):

| Layer | single user | group | system-wide |
|---|---|---|---|
| extra-root | — | — | writes (group) |
| app-root | writes | writes (group) | read-only; writes (admin) |
| scratch | <span style="color:#888">fallback</span> | <span style="color:#888">fallback</span> | writes (solo) |
| user | <span style="color:#888">fallback</span> | <span style="color:#888">fallback</span> | <span style="color:#888">fallback</span> |

Reads always search every active layer; only the write target differs. *fallback* is used only when the higher writable layers are unset.

Each one holds the same three subdirectories:

```
<layer>/
├── build-scripts/
├── helper-scripts/
└── images/
```

## Reads Search Everything, Writes Take the First Writable

This is the rule that makes shared installs work:

- **Reading** — `condatainer avail`, `list`, `exec`, `helper` search *every* layer. A user sees the admin's build scripts, their lab's overlays, and their own, merged into one view.
- **Writing** — `condatainer create` walks the same order and uses the **first layer it can write to**.

So a read-only system install is simply skipped when a user builds something, and their overlay lands in their scratch or home instead.

The corollary: if a shared layer *is* writable by a user, their builds go **there** rather than into their own directory. That is the intended behaviour when you want members contributing shared builds. Permissions are what decide it — see [Shared Group Installation](./group_install.md#1-pick-a-directory).

`create` reports the layer it chose:

```
[CNT◇] Installing to /shared/labA/condatainer/images (extra-root)
```

An admin who expects `(app-root)` and sees `(user)` has just installed for themselves only — usually the wrong account.

```{tip}
`condatainer config paths` prints every directory in priority order, tagged with its layer and marked with which one receives writes:

    Images:
      1. /shared/labA/condatainer/images (extra-root) (writable, target)
      2. /opt/condatainer/images (app-root) (read-only)
      3. /scratch/me/condatainer/images (user) (not found)

It is the fastest way to check a deployment from a normal user's account.
```

## Working With a Specific Layer

The same overlay name can exist in more than one layer — a user rebuilding something the admin already ships, for example. `condatainer list` groups its output by directory and tags each with its layer, so both copies are visible:

```
════ /opt/condatainer/images (app-root) ════
Available app overlays:
 samtools/1.22

════ /scratch/me/condatainer/images (user) ════
Available app overlays:
 samtools/1.22
```

Commands that resolve a name use the highest-priority copy. To act on a different one, scope with `-l` — the same `u` / `r` / `e` vocabulary as `condatainer config -l`:

```bash
condatainer remove -l u samtools/1.22   # remove your copy, leave the site's alone
```

Without a scope, `remove` refuses to act on a name that exists in several layers and lists them, rather than guessing.

## How app-root Is Found

The app-root layer is the installation's own directory. It is resolved in one of two ways:

1. **From the binary** — a binary at `<cnt_root>/bin/condatainer` makes `<cnt_root>` the app-root. Symlinks are resolved first, so a link from `/usr/local/bin` still points back to the real install.
2. **From `$CNT_ROOT`** — set the variable and it wins outright, whatever the binary's location.

Shared binary directories (`$HOME/bin`, `~/.local/bin`, `/usr/bin`, `/usr/local/bin`) are deliberately excluded: they hold binaries for many programs, so their parent is nobody's install root. A binary installed there has no app-root unless `$CNT_ROOT` supplies one.

`$CNT_ROOT` also decouples the data from the binary, which matters for versioned module trees — see [System-Wide Installation](./system_install.md#with-environment-modules).

## Config Files Layer Too

Every layer can carry a `config.yaml`, and all of them are loaded together:

- **Scalar keys** (`apptainer_bin`, `submit_job`, `build.ncpus`) — the highest-priority layer that sets the key wins. A user's own config overrides the site default.
- **Directory and source keys** (`extra_*_dirs`, `extra_scripts_links`) — **merged** across layers, user entries first, then `extra-root`, `app-root`.
- **`channels`** — overwrite, not merged: Conda channel order decides which package wins, so merging two lists would silently change resolution.

The split is deliberate: users can override *settings* without being able to remove *sources* an admin published. See [Configuration Priority](../manuals/configuration.md#configuration-priority).

## Reference

- [Configuration Manual](../manuals/configuration.md#data-directory-search-paths) — the complete search-path and config-key tables
- [Multi-Tier Setup](../manuals/configuration.md#multi-tier-setup-system--group--user) — a worked system/group/user example
