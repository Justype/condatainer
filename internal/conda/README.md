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
```

The `condatainer env` commands are wired in `cmd/env.go`; `/usr/bin/mm` is only an in-container
shortcut to that command group.
