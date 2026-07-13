# VS Code

Three variants are available:

| Variant | Port Forwarding | Extensions | Best For |
|---------|-----------------|------------|----------|
| `code-server` | Required | Open VSX only | Open-source alternative, no Microsoft account |
| `vscode-server` | Required | All (Copilot, Pylance) | Full VS Code, direct connection |
| `vscode-tunnel` | Not required | All (Copilot, Pylance) | When port forwarding is difficult |

```bash
condatainer helper -u
condatainer helper code-server      # open-source
condatainer helper vscode-server    # full VS Code
condatainer helper vscode-tunnel    # no port forwarding
```

See [Helper Scripts](./helpers.md) for SSH port forwarding setup, resource flags, configuration, and reuse mode.

## code-server

Open-source VS Code served in the browser. Does **not** support Microsoft-proprietary extensions (Copilot, Pylance). Use `vscode-server` if you need those.

### Flags

```bash
condatainer helper code-server --help
```

| Flag | Description | Default |
|------|-------------|---------|
| `--password,-p` | Browser password (empty = no auth) | (empty) |

See [Helper Scripts](./helpers.md#flags) for resource and overlay flags.

### Issues

**Some extensions may not render correctly in Firefox** — if an extension's webview Content Security Policy omits the webview resource origin from `style-src`, Firefox blocks external stylesheets while Chromium silently exempts `vscode-cdn.net` URLs from CSP. Use Chromium/Chrome for affected extensions.

## vscode-server

Uses the official VS Code CLI (`code`) to serve a full VS Code instance in your browser. All extensions — including Microsoft-proprietary ones (Copilot, Pylance, Remote, etc.) — are available.

### Flags

```bash
condatainer helper vscode-server --help
```

| Flag | Description | Default |
|------|-------------|---------|
| `--token,-a` | Connection token for the web UI | (auto-generated) |

See [Helper Scripts](./helpers.md#flags) for resource and overlay flags.

## vscode-tunnel

Uses Microsoft's relay service — no SSH port forwarding needed. Authenticate once with GitHub or Microsoft, then connect your local VS Code by tunnel name.

```
Local VS Code ←→ Microsoft relay ←→ Compute node
```

```{note}
`vscode-tunnel` is a singleton helper — only one instance can run at a time (tunnel names must be unique per account).
```

### Flags

```bash
condatainer helper vscode-tunnel --help
```

| Flag | Description | Default |
|------|-------------|---------|
| `--auth,-a` | Auth provider: `github` or `microsoft` | `github` |
| `--name,-n` | Tunnel machine name | `<username>-<hostname>` |

See [Helper Scripts](./helpers.md#flags) for resource and overlay flags.

### Authentication

On first run (or after token expiry):

```
To sign in, use a web browser to open the page https://microsoft.com/devicelogin
and enter the code FE2G6WJQK to authenticate.
```

Authentication times out after 10 minutes — run again if it expires.

The access URL is `https://vscode.dev/tunnel/<name>/...` — connect via the [VS Code Remote - Tunnels extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode.remote-server) or open it in a browser.

## Common Issues

### Inode quota issues

VS Code stores data under `$HOME` by default, which can exhaust the HOME inode quota on HPC systems. Redirect to high-inode storage (SCRATCH, VAST, or similar) with symlinks — set this up once on the login node:

```bash
# Replace /path/to/high-inode-storage with your SCRATCH, VAST, etc.
HIGH_INODE=/path/to/high-inode-storage

# code-server: extensions
mkdir -p "$HIGH_INODE/.local/share/code-server"
ln -sfn "$HIGH_INODE/.local/share/code-server" "$HOME/.local/share/code-server"

# vscode-server and vscode-tunnel: CLI runtime
mkdir -p "$HIGH_INODE/.vscode"
ln -sfn "$HIGH_INODE/.vscode" "$HOME/.vscode"

# vscode-server: extensions
mkdir -p "$HIGH_INODE/.vscode-server"
ln -sfn "$HIGH_INODE/.vscode-server" "$HOME/.vscode-server"
```

### Connection token issues (vscode-server)

The token is embedded in the URL as `tkn=<token>`. Copy the full URL from the output.

### Machine name already taken (vscode-tunnel)

```bash
condatainer helper vscode-tunnel config set MACHINE_NAME my-unique-name
```
