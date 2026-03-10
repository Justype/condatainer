# CondaTainer Tutorials

**CondaTainer** keeps your tools, data, and environments on the HPC cluster as overlay files — no manual downloads, no upload/download round-trips between your laptop and the cluster.

## Part 1: Managing Overlays

Choose the right overlay type for your workflow:

| Guide | Best For |
|-------|----------|
| [Module Overlays](../user_guide/module_overlays.md) | Installing individual tools (`samtools`, `salmon`) and reference datasets |
| [Workspace Overlays](../user_guide/workspace_overlays.md) | Writable Conda environments for interactive development |
| [Bundle Overlays](../user_guide/bundle_overlays.md) | Frozen, read-only Conda environments for sharing and reproducibility |

## Part 2: Upstream Analyses

Run HPC analysis pipelines with automatic dependency resolution, scheduler submission, and job chaining — all from a single `condatainer run` command.

Check out this repo: [cnt-tutorials](https://github.com/Justype/cnt-tutorials) for examples.

See the [CondaTainer Manual: run](../manuals/condatainer.md#runtime-check-run) for a full reference on flags, array jobs, job chaining, and scheduler integration.

## Part 3: Downstream Tools / IDEs

Run interactive analysis environments on HPC compute nodes.

### VS Code

| Tool | Port Forwarding | Best For |
|------|-----------------|----------|
| [VS Code Server](./vscode-server_on_HPC.md) | Required | Full VS Code with all extensions and Copilot |
| [VS Code Tunnel](./vscode-tunnel_on_HPC.md) | Not required | Full VS Code via Microsoft/GitHub relay — when port forwarding is difficult |

For non-HPC headless servers, see [Helper Scripts on Headless Servers](./helpers_on_server.md).

### RStudio Server

| Tool | R Installation | Best For |
|------|----------------|----------|
| [RStudio Server](./rstudio-server_on_HPC.md) | Posit R image overlays | Reproducible R versions, system-level R |
| [RStudio Server (Conda)](./rstudio-server-conda_on_HPC.md) | Conda (`mm-install r-base`) | Conda-managed R with packages |

For non-HPC remote servers, see [RStudio Server on Headless Server](./rstudio-server_on_server.md).

### Jupyter

Run Jupyter Lab or Notebook on HPC compute nodes:

```bash
condatainer helper jupyter
```

For non-HPC headless servers, see [Helper Scripts on Headless Servers](./helpers_on_server.md).

### Desktop / GUI Apps

| Tool | Description |
|------|-------------|
| [xfce4](./xfce4_on_HPC.md) | XFCE desktop session via VNC / noVNC, accessible through a browser |
| [igv](./xfce4_on_HPC.md) | XFCE desktop with IGV (Integrative Genomics Viewer) via noVNC |

For non-HPC headless servers, see [Helper Scripts on Headless Servers](./helpers_on_server.md).

---

See [Helper Scripts on HPC](./helpers_on_HPC.md) for how helper scripts work with overlays and SSH port forwarding.
