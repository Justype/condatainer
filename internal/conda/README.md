# conda

Isolated, in-container management of the Conda environment mounted at `/cnt_env`.

## Isolation contract

- `CNT_CONDA_ROOT` identifies the mounted environment; mutations also require
  `CNT_CONDA_WRITABLE=1`.
- `CNT_CONDA_CHANNELS` carries CondaTainer's configured first-install channels.
- `CONDA_PREFIX` and `MAMBA_ROOT_PREFIX` are reset to the mounted environment for every command.
- Ambient Conda and Micromamba prefix and rc variables are removed.
- With `/cnt_env/.condarc`, Micromamba runs with `--rc-file` so no other rc files are loaded.
- Without `.condarc`, Micromamba runs with `--no-rc` during env creation.
- The first install uses explicit project channels, then saves them to `.condarc` after success.
- User Micromamba arguments are forwarded unchanged. Micromamba owns transaction locking.

The default first-install channels come from CondaTainer's `channels` configuration. Explicit CLI
channels are preserved ahead of those defaults; an environment file owns its channel list.

## Package layout

```text
environment.go  Validate the container, prefix, and writable state
runner.go       Isolate configuration and invoke Micromamba
config.go       Read and atomically update .condarc channels
pinned.go       Read and atomically update conda-meta/pinned
export.go       Canonical explicit.txt and environment.yml for an image to embed
packages.go     List installed packages, from a writable .img or a read-only .sqf
image.go        Read conda-meta/history for channels and explicitly-installed specs
```

`packages.go` and `image.go` are host-side reads (used by `info`, the dashboard, and
helper `#IMG_PACKAGES:` checks), unlike the rest of this package's in-container
contract above. `ListCondaPackages` reads a `.img`'s own `conda-meta` via `debugfs`;
`ListCondaPackagesSqf` is its read-only counterpart for a `.sqf`, via `unsquashfs -l`.
`ReadCondaInfo` reads one `conda-meta/history`; `ReadCondaInfoMerged` concatenates a
snapshot's history with a `.img`'s own and parses them as one continuous log, since a
`.img`'s history is genuinely the continuation of the snapshot's history it was frozen
from.

**This package itself never decides whether a pairing applies** — every function above
takes exactly the path(s) it is given and reads only those. A writable `.img` paired with
a frozen snapshot (`container.LookupSnapshot`) holds only the newest delta in its own
`conda-meta` and `conda-meta/history`, so reading it alone would report almost everything
the snapshot already provides as missing, or a nearly-empty install history — but knowing
*whether* a given path has a pair requires `container.LookupSnapshot`, which this package
cannot import (`container` depends on `internal/artifact/key`, which depends on `conda`).
`container.PairedPackages` and `container.PairedInfo` are where that decision is made:
they take a bare path, find its pair if one exists, and call the right function(s) here.
A caller that wants "packages installed at this path" or "history for this path" should
call those, not `ListCondaPackages`/`ReadCondaInfo` directly — see
[`internal/runtime/container/README.md`](../runtime/container/README.md), *Environment
Snapshots*.

## The two exports

A Conda app embeds `explicit.txt` and `environment.yml`, captured from the
environment that was **actually installed** rather than from a second solve.
Together they are the artifact's identity and equivalence, so `sha256sum` on
either reproduces a key by hand.

| file | from | pins |
|---|---|---|
| `explicit.txt` | `env export --explicit --no-md5` | exact package URLs — channel, subdir, name, version, build string |
| `environment.yml` | `env export --no-builds` | channels, package names and versions |

Both are **re-emitted, not stored as the tool printed them**. That is the whole
point: a Micromamba upgrade that reorders or re-spaces its output would otherwise
move the key of every artifact built after it, and every comparison across that
boundary would report *different* with an empty diff. Sorting is safe — Conda
does not depend on explicit-file order — and the result is still a working
`micromamba install --file` / `create -f` input.

`name:` and `prefix:` are dropped: they describe a temporary local environment
rather than its packages. A trailing `#<md5>` is dropped too — it checksums the
download, not what was installed, and the URL already distinguishes one package
from another.

**The channel set is the export's; only the order is corrected.** Those are the
channels that actually provided packages. `OrderChannels` reorders them into the
configured priority, because Micromamba alphabetizes on export and its sequence is
not the solve order; a channel used but not configured — a mirror, a
`bioconda::star` annotation — is appended after the known ones, sorted, rather
than dropped. So adding an unused channel to a site's config changes no
artifact's key, and a rebuild still resolves against the channels that mattered.

The split is deliberately asymmetric: changing a build string moves identity and
not equivalence, changing a version moves both. `numpy=1.26.4` built against MKL
and against OpenBLAS are equivalent — that is the intended answer, and identity
catches the difference for anyone who needs it.

A `pip:` sub-list is preserved and sorted as exported, and nothing chases it
further. `--explicit` omits pip packages entirely, so an environment with them has
an identity that ignores them; a centrally managed overlay does not carry pip
installs, which is why that is left alone rather than papered over.

The `condatainer env` commands are wired in `cmd/env.go`; `/usr/bin/mm` is only an in-container
shortcut to that command group.
