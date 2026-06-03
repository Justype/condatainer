# XFCE Desktop (VNC)

The `xfce4` helper launches an XFCE desktop session accessible in your browser via VNC. The `igv` and `cytoscape` helpers do the same with [IGV](#igv) or [Cytoscape](#cytoscape) pre-launched inside the desktop.

```bash
condatainer helper -u
condatainer helper xfce4
```

See [Helper Scripts](./helpers.md) for SSH port forwarding setup, resource flags, configuration, and reuse mode.

## VNC Provider

Choose between two VNC backends with `--vnc,-v`:

```bash
condatainer helper xfce4 --vnc kasm    # KasmVNC (default)
condatainer helper xfce4 --vnc turbo   # TurboVNC + noVNC
```

### TurboVNC

Uses TurboVNC with websockify as a noVNC WebSocket bridge. The session is password-protected with an auto-generated password embedded in the access URL.

### KasmVNC

Uses KasmVNC's built-in WebSocket server — no separate bridge needed. Browser-native protocol with better compression and additional features:

- **Clipboard sync**: bidirectional clipboard between browser and desktop.
- **Audio**: PulseAudio null sink streams desktop audio directly to your browser tab via kclient.
- **GPU acceleration**: Automatically detects NVIDIA EGL for hardware rendering. Falls back to software rendering on MIG-partitioned GPUs.

## GPU

Request a GPU when submitting the job:

```bash
condatainer helper xfce4 -g h100:1
```

### With KasmVNC (`--vnc kasm`)

GPU detection is automatic. On nodes where NVIDIA EGL is available (full GPU, non-MIG), hardware rendering is enabled and reported:

```
[CNT] GPU: using /dev/dri/renderD128 for hardware acceleration
```

On MIG-partitioned GPUs, EGL device enumeration is unavailable, so software rendering is used instead:

```
[CNT] GPU: DRI node found but NVIDIA EGL unavailable (MIG?), using software rendering
```

The DRI node is used for KasmVNC's own screen encoding pipeline. **CUDA/OptiX** workloads (e.g. Blender Cycles) work as normal via Apptainer's `--nv` injection. Apps using **EGL directly** also get GPU access. However, traditional **GLX-based** OpenGL apps still see software rendering (llvmpipe) inside the VNC session — `glxinfo` will show Mesa even with hw3d enabled. Use `vglrun -d egl <app>` for those.

### With TurboVNC (`--vnc turbo`)

TurboVNC does not use the GPU for the desktop session itself. To run a specific GPU-accelerated application inside the desktop, prefix it with VirtualGL:

```bash
vglrun -d egl <application>
```

Example — launch Blender with GPU rendering:

```bash
vglrun -d egl blender
```

## Flags

```bash
condatainer helper xfce4 --help
```

| Flag | Description | Default |
|------|-------------|---------|
| `--vnc,-v` | VNC provider: `kasm` or `turbo` | `kasm` |

See [Helper Scripts](./helpers.md#flags) for resource and overlay flags.

## IGV

The `igv` helper is identical to `xfce4` but automatically installs and launches IGV (Integrative Genomics Viewer) inside the desktop session.

```bash
condatainer helper igv --igv 2.19.7
```

| Flag | Description | Default |
|------|-------------|---------|
| `--igv,-i` | IGV version | lastest conda IGV version |
| `--vnc,-v` | VNC provider: `kasm` or `turbo` | `kasm` |

See [Helper Scripts](./helpers.md#flags) for resource and overlay flags.

## Cytoscape

The `cytoscape` helper is identical to `xfce4` but automatically installs and launches [Cytoscape](https://cytoscape.org/) (network biology visualization platform) inside the desktop session.

```bash
condatainer helper cytoscape --cytoscape 3.10.4 --openjdk 17.0.18
```

| Flag | Description | Default |
|------|-------------|---------|
| `--cytoscape,-C` | Cytoscape version | latest Cytoscape version |
| `--openjdk,-j` | OpenJDK version (must be 17.x) | latest 17.x |
| `--vnc,-v` | VNC provider: `kasm` or `turbo` | `kasm` |

See [Helper Scripts](./helpers.md#flags) for resource and overlay flags.

### Java version

Cytoscape 3.x requires **Java 17**. The `--openjdk` flag is locked to 17.x; using any other major version will prevent Cytoscape from starting.

### CytoscapeConfiguration

Cytoscape writes its configuration and bundle cache to `~/CytoscapeConfiguration/`. To keep it off your `$HOME` quota, symlink it to `$SCRATCH` before starting:

```bash
ln -s $SCRATCH/CytoscapeConfiguration ~/CytoscapeConfiguration
```

## Desktop Home Directory

XDG directories (Desktop, Documents, Downloads, etc.) are redirected to `$SCRATCH/desktop-home/` by default. This avoids filling your `$HOME` quota with desktop files across sessions.

## Common Issues

### Screen resolution

The default geometry is 1600×900. For TurboVNC, the noVNC URL includes `resize=remote` so the desktop scales to your browser window. For KasmVNC, scaling is handled natively.

### Audio not working (KasmVNC)

PulseAudio must be installed in the `xfce4` overlay. If the null sink fails to load, a warning is printed and audio is unavailable for that session.

