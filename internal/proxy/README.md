# proxy

SSH tunnel + dual-protocol proxy management for HPC compute nodes.

## Two Modes

| | Shared | Per-job |
|---|---|---|
| Started on | Login node | Compute node (inside job) |
| Bind | `0.0.0.0:PORT` | `127.0.0.1:PORT` |
| PID file | `config.GetUserDataDir()/proxy.pid` (NFS) | `utils.GetTmpDir()/proxy.pid` (local) |
| Lifetime | Until stopped | Until job ends |

## Daemon Architecture

```
job/login node
  RunDaemon()
    ├── tryEstablishTunnel()   ← SSH tunnel to login node
    │     ├── 1. DialGoSSH()        Go SSH library (no subprocess)
    │     ├── 2. ssh -D unix-sock   system ssh, Unix socket SOCKS5
    │     └── 3. ssh -D 127.0.0.1  system ssh, TCP SOCKS5 (OpenSSH <6.7)
    └── RunProxy(ctx, listenAddr, dial)
          ├── SOCKS5 (RFC 1928)      detected by first byte 0x05
          └── HTTP CONNECT           everything else
```

`RunProxy` serves both protocols on the same port. The `DialFunc` from whichever tunnel method succeeded is passed directly — no intermediate SOCKS5 layer for the go-ssh path.

## SSH Tunnel Fallback Chain (`daemon.go`, `ssh_dial.go`)

`tryEstablishTunnel` tries three methods in order, returning on the first success:

| Priority | Method | Condition |
|---|---|---|
| 1 | `DialGoSSH` — Go SSH library | Any non-interactive auth available |
| 2 | `ssh -D /path/to.sock` | system `ssh` binary in PATH |
| 3 | `ssh -D 127.0.0.1:PORT` | fallback for OpenSSH < 6.7 (no Unix socket `-D`) |

## Go SSH Auth Methods (`ssh_dial.go`)

`buildAuthMethods` assembles all available non-interactive auth methods — never prompts:

| Method | Condition |
|---|---|
| `ssh.Hostbased` | `ssh-keysign` setuid binary found + host key in `/etc/ssh/` |
| `ssh.PublicKeysCallback` via agent | `$SSH_AUTH_SOCK` is set and connectable |
| `ssh.PublicKeys` via key files | `~/.ssh/id_{ed25519,ecdsa,ecdsa_sk,ed25519_sk,rsa}` readable |

RSA keys are wrapped in `rsaSHA2Signer` to advertise and sign with `rsa-sha2-256` consistently (fixes a `golang.org/x/crypto/ssh` mismatch where `Type()` returns `ssh-rsa` but `Sign()` uses SHA-256).

The Go SSH client sends a `keepalive@openssh.com` request every 60 s to prevent firewalls from silently dropping idle connections.

## Key Functions

```go
proxy.GetJobProxy() (string, bool)                   // cached; returns "", false on login nodes
proxy.FindActiveProxy(autostart bool) (string, bool) // per-job → shared → auto-start
proxy.ResolveViaHost(flagVia string) (string, error) // --via > CNT_PROXY_VIA > error
proxy.StartOnNode(host, via string, port int) error  // SSH delegation; waits for NFS PID file
proxy.ProxyAlive(host string, port int) bool         // TCP dial, 500ms timeout
proxy.RunDaemon(sshDest string, port int, localOnly bool, reportFd int) error
proxy.DialGoSSH(sshDest string) (DialFunc, func(), <-chan struct{}, error)
```

## Usage

`proxy.GetJobProxy()` is the call site for injection — resolves once per process via `sync.Once`,
no-ops on login nodes. Used by `internal/exec/run.go` (container env vars) and
`internal/build/fetch.go` (HTTP transport).

## Daemon Readiness Protocol

`proxy start` and `autoStartLocalProxy` use an `os.Pipe` to synchronise with the daemon:

1. Parent creates `(pr, pw)`, passes `pw` as `ExtraFiles[0]` (fd 3 in daemon) and `--report-fd 3`
2. Daemon closes `pw` (EOF) once tunnel is up and PID file is written → success
3. Daemon writes an error message before closing → startup failed, surfaced to user
4. Parent blocks on `io.ReadAll(pr)` — no polling needed

## CNT_PROXY_VIA

Baked at job submission time by `condatainer build`/`run` using `os.Hostname()` on the login node.
Used by `ResolveViaHost` so compute nodes know which login node to tunnel through.
