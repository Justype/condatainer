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
func (r *proxyRegistry) Open(id, node string, port int) {
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
		r.mu.Lock()
		delete(r.pending, id)
		if _, ok := r.entries[id]; !ok {
			r.entries[id] = &proxyEntry{rp: rp, node: node, port: port}
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
	rp.Transport = transport

	r.entries[id] = &proxyEntry{rp: rp, node: node, port: port}
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

// LastError returns the most recent Open() failure reason for a helper, or "".
func (r *proxyRegistry) LastError(id string) string {
	r.mu.RLock()
	defer r.mu.RUnlock()
	return r.lastErrors[id]
}

// ServeProxy handles /proxy/{id}/{path...} by forwarding to the helper's tunnel.
func (s *srv) ServeProxy(w http.ResponseWriter, r *http.Request) {
	// Path: /proxy/{id}/{rest...}
	parts := strings.SplitN(strings.TrimPrefix(r.URL.Path, "/proxy/"), "/", 2)
	if len(parts) == 0 || parts[0] == "" {
		http.Error(w, "missing helper ID", http.StatusBadRequest)
		return
	}
	id := parts[0]

	entry := s.proxies.Get(id)
	if entry == nil {
		msg := "no active tunnel for " + id
		if reason := s.proxies.LastError(id); reason != "" {
			msg += ": " + reason
		}
		http.Error(w, msg, http.StatusServiceUnavailable)
		return
	}

	// Rewrite the request path to strip the /proxy/{id} prefix.
	rest := "/"
	if len(parts) > 1 {
		rest = "/" + parts[1]
	}
	r2 := r.Clone(r.Context())
	r2.URL.Path = rest
	r2.URL.RawPath = rest
	if r.URL.RawQuery != "" {
		r2.URL.RawQuery = r.URL.RawQuery
	}
	entry.rp.ServeHTTP(w, r2)
}
