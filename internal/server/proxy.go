package server

import (
	"fmt"
	"log/slog"
	"net/http"
	"net/http/httputil"
	"net/url"
	"os"
	"strings"
	"sync"

	internalproxy "github.com/Justype/condatainer/internal/proxy"
)

// proxyEntry holds a live reverse-proxy for one helper run.
type proxyEntry struct {
	rp   *httputil.ReverseProxy
	name string // helper script name (e.g. "code-server")
	node string
	port int
}

// proxyRegistry caches one reverse-proxy per helper ID.
type proxyRegistry struct {
	mu         sync.RWMutex
	entries    map[string]*proxyEntry
	pending    map[string]struct{} // IDs with an in-flight Open() dial
	lastErrors map[string]string   // last Open() failure reason per ID
	log        *slog.Logger
}

func newProxyRegistry(log *slog.Logger) *proxyRegistry {
	return &proxyRegistry{
		entries:    make(map[string]*proxyEntry),
		pending:    make(map[string]struct{}),
		lastErrors: make(map[string]string),
		log:        log,
	}
}

// isLocalHost returns true when node refers to the current machine, meaning no
// SSH tunnel is needed (headless / same-node jobs).
func isLocalHost(node string) bool {
	if node == "localhost" || node == "127.0.0.1" || node == "::1" {
		return true
	}
	if h, err := os.Hostname(); err == nil && h == node {
		return true
	}
	return false
}

// Open creates (or returns existing) a reverse-proxy for the given helper using
// a three-level fallback: direct TCP (same host) → Go SSH → system ssh -W.
// The SSH dial happens outside the registry lock so HTTP Get() calls are never
// blocked by a slow or failing dial. Concurrent Open() calls for the same ID
// are deduplicated via the pending set.
func (r *proxyRegistry) Open(id, name, node string, port int) {
	r.mu.Lock()
	if _, ok := r.entries[id]; ok {
		r.mu.Unlock()
		return
	}
	if _, ok := r.pending[id]; ok {
		r.mu.Unlock()
		return // another goroutine is already dialing
	}
	r.pending[id] = struct{}{}
	r.mu.Unlock()

	// Level 1: same-host — direct TCP, no SSH needed (headless jobs).
	if isLocalHost(node) {
		target, _ := url.Parse(fmt.Sprintf("http://127.0.0.1:%d", port))
		rp := httputil.NewSingleHostReverseProxy(target)
		rp.Director = hostRewriteDirector(target)
		r.mu.Lock()
		delete(r.pending, id)
		if _, ok := r.entries[id]; !ok {
			r.entries[id] = &proxyEntry{rp: rp, name: name, node: node, port: port}
			delete(r.lastErrors, id)
			r.log.Debug("server: proxy direct TCP", "id", id, "port", port)
		}
		r.mu.Unlock()
		return
	}

	// SSH dial happens outside the lock — may take up to 25 s on unreachable nodes.
	// Level 2: Go SSH.
	dial, stop, done, err := internalproxy.DialGoSSH(node)
	var goSSHErr error
	if err != nil {
		goSSHErr = err
		r.log.Debug("server: go-ssh failed", "node", node, "err", err)
		r.log.Info("server: falling back to system ssh", "node", node)
		// Level 3: system ssh -W fallback.
		dial, stop, done, err = internalproxy.DialSystemSSH(node)
	}

	r.mu.Lock()
	defer r.mu.Unlock()
	delete(r.pending, id)

	if err != nil {
		r.log.Warn("server: cannot open tunnel", "node", node, "id", id, "err", err)
		if goSSHErr != nil {
			r.lastErrors[id] = fmt.Sprintf("SSH to %s failed: %v", node, goSSHErr)
		} else {
			r.lastErrors[id] = fmt.Sprintf("SSH to %s failed: %v", node, err)
		}
		return
	}
	// Double-check: another goroutine may have opened it while we were dialing.
	if _, already := r.entries[id]; already {
		stop()
		return
	}

	target, _ := url.Parse(fmt.Sprintf("http://127.0.0.1:%d", port))
	transport := &http.Transport{
		DialContext: dial,
	}
	// Watch for unexpected tunnel closure and clean up.
	go func() {
		<-done
		r.Close(id)
		stop()
	}()
	rp := httputil.NewSingleHostReverseProxy(target)
	rp.Director = hostRewriteDirector(target)
	rp.Transport = transport

	r.entries[id] = &proxyEntry{rp: rp, name: name, node: node, port: port}
	delete(r.lastErrors, id)
	r.log.Debug("server: proxy tunnel opened", "id", id, "node", node, "port", port)
}

// Close removes and discards the proxy entry for a helper.
func (r *proxyRegistry) Close(id string) {
	r.mu.Lock()
	delete(r.entries, id)
	delete(r.pending, id)
	delete(r.lastErrors, id)
	r.mu.Unlock()
}

// Get returns the proxy entry for a helper, or nil if not open.
func (r *proxyRegistry) Get(id string) *proxyEntry {
	r.mu.RLock()
	defer r.mu.RUnlock()
	return r.entries[id]
}

// GetByName returns the oldest active proxy entry whose helper name matches,
// or nil if none is open. Used for stable-name subdomain routing
// (e.g. code-server.localhost → oldest running code-server instance).
// IDs have the form "{name}-YYYYMMDD-HHMMSS", so the lexicographically smallest
// ID is the oldest one; this ensures the name URL stays stable when a newer
// instance starts alongside an existing one.
func (r *proxyRegistry) GetByName(name string) *proxyEntry {
	r.mu.RLock()
	defer r.mu.RUnlock()
	var oldestID string
	var oldest *proxyEntry
	for id, entry := range r.entries {
		if entry.name == name && (oldestID == "" || id < oldestID) {
			oldestID = id
			oldest = entry
		}
	}
	return oldest
}

// LastError returns the most recent Open() failure reason for a helper, or "".
func (r *proxyRegistry) LastError(id string) string {
	r.mu.RLock()
	defer r.mu.RUnlock()
	return r.lastErrors[id]
}

// serveSubdomainProxy handles requests to {id}.localhost:{port} by forwarding
// to the helper's tunnel. The full request path is forwarded as-is — no URL
// rewriting needed since the service is at the root of its own subdomain.
// Lookup order: exact ID match first, then first active instance by name
// (e.g. "code-server.localhost" → first running code-server).
func (s *srv) serveSubdomainProxy(w http.ResponseWriter, r *http.Request, id string) {
	entry := s.proxies.Get(id)
	if entry == nil {
		entry = s.proxies.GetByName(id)
	}
	if entry == nil {
		msg := "no active tunnel for " + id
		if reason := s.proxies.LastError(id); reason != "" {
			msg += ": " + reason
		}
		http.Error(w, msg, http.StatusServiceUnavailable)
		return
	}
	entry.rp.ServeHTTP(w, r)
}

// hostRewriteDirector returns a director func that rewrites req.Host and the
// Origin header to the backend address. This prevents host-header and
// WebSocket origin checks in apps like JupyterLab and code-server (which
// reject requests whose Host/Origin doesn't match their bound address) from
// blocking proxied requests that arrive with a subdomain Host header.
func hostRewriteDirector(target *url.URL) func(*http.Request) {
	std := httputil.NewSingleHostReverseProxy(target).Director
	origin := "http://" + target.Host
	return func(req *http.Request) {
		std(req)
		req.Host = target.Host
		if req.Header.Get("Origin") != "" {
			req.Header.Set("Origin", origin)
		}
	}
}

// hostDispatch wraps the dashboard mux with subdomain-based proxy routing.
// Requests to {id}.localhost:{port} are forwarded to the matching helper tunnel;
// all other requests (localhost:{port}) are handled by the dashboard mux.
func (s *srv) hostDispatch(next http.Handler) http.Handler {
	return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		host := r.Host
		// Strip port suffix.
		if i := strings.LastIndex(host, ":"); i >= 0 {
			host = host[:i]
		}
		// {id}.localhost → proxy to helper.
		const suffix = ".localhost"
		if strings.HasSuffix(host, suffix) {
			id := host[:len(host)-len(suffix)]
			if id != "" {
				s.serveSubdomainProxy(w, r, id)
				return
			}
		}
		next.ServeHTTP(w, r)
	})
}
