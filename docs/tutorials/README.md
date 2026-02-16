# CondaTainer Tutorials

This directory contains tutorials on how to use CondaTainer for various applications.

## Getting Started

- [Create a Workspace Environment](../user_guide/workspace_overlays.md)
- [Helper Scripts on HPC](./helpers_on_HPC.md) - How helper scripts work with overlays and SSH port forwarding

## IDE & Editor Tutorials

### VS Code Options

| Tool | Port Forwarding | Features | Best For |
|------|-----------------|----------|----------|
| [VS Code Server](./vscode-server_on_HPC.md) | Required | Full VS Code, all extensions, Copilot | Full-featured development |
| [VS Code Tunnel](./vscode-tunnel_on_HPC.md) | Not required | Full VS Code via Microsoft/GitHub relay | When port forwarding is difficult |

Now [code-server](./code-server_on_HPC.md) is replaced by VS Code Server. (Pylance, Copilot, and other Microsoft extensions are supported only on VS Code.)

### RStudio Server

| Tool | R Installation | Best For |
|------|----------------|----------|
| [RStudio Server](./rstudio-server_on_HPC.md) | Posit R image overlays | Reproducible R versions, system-level R |
| [RStudio Server (Conda)](./rstudio-server-conda_on_HPC.md) | Conda (`mm-install r-base`) | Conda-managed R with packages |

For non-HPC remote servers, see [RStudio Server on Headless Server](./rstudio-server_on_server.md).

### Jupyter

Run Jupyter Lab or Notebook on HPC. Use `-j notebook` for classic Notebook mode.

```bash
condatainer helper jupyter
```

### Desktop / GUI Apps

| Tool | Description |
|------|-------------|
| [xfce4](./xfce4_on_HPC.md) | XFCE desktop session via VNC / noVNC, accessible through a browser |
| [igv](./xfce4_on_HPC.md) | XFCE desktop with IGV (Integrative Genomics Viewer) via noVNC |
