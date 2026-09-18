# Self-Provisioned Toolchain Troubleshooting

This applies only when you are using the self-provisioned toolchain (`mksquashfs`, `squashfuse`,
`apptainer` installed by `condatainer update --libexec <package>`), not a system-installed or module-loaded
Apptainer.

## Mount or unmount errors mentioning `fusermount3`

If `exec`/`run`/`build` fails with apptainer output like:

```
Cleanup error: while stopping driver for .../mnt/session/rootfs: squashfuse_ll exited: Failed to call 'fusermount3': No such file or directory
Spawning fusermount3 to unmount failed: No such file or directory
```

the self-provisioned toolchain is out of date. Refresh it:

```bash
condatainer update --libexec
```

If any other condatainer session is currently running, stop it first — `update --libexec` refuses
to run while the toolchain is in use, rather than waiting.
