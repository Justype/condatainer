# proxy

SSH SOCKS5 proxy tunnel management for HPC compute nodes.

## Two Modes

| | Shared | Per-job |
|---|---|---|
| Started on | Login node | Compute node (inside job) |
| Bind | `0.0.0.0:PORT` | `127.0.0.1:PORT` |
| PID file | `config.GetUserDataDir()/proxy.pid` (NFS) | `utils.GetTmpDir()/proxy.pid` (local) |
| Lifetime | Until stopped | Until job ends |

## Key Functions

```go
proxy.GetJobProxy() (string, bool)               // cached; returns "", false on login nodes
proxy.FindActiveProxy(autostart bool) (string, bool) // per-job → shared → auto-start
proxy.ResolveViaHost(flagVia string) (string, error) // --via > CNT_PROXY_VIA > error
proxy.StartOnNode(host, via string, port int) error  // SSH delegation; waits for NFS PID file
proxy.ProxyAlive(host string, port int) bool         // TCP dial, 500ms timeout
```

## Usage

`proxy.GetJobProxy()` is the call site for injection — resolves once per process via `sync.Once`,
no-ops on login nodes. Used by `internal/exec/run.go` (container env vars) and
`internal/build/fetch.go` (HTTP transport).

## CNT_PROXY_VIA

Baked at job submission time by `condatainer build`/`run` using `os.Hostname()` on the login node.
Used by `ResolveViaHost` so compute nodes know which login node to tunnel through.
