package registry

import (
	"fmt"
	"net"
	"strings"

	"oras.land/oras-go/v2/registry/remote"
)

// newRepository builds an authenticated client for "<base>/<repo>", with no tag —
// a caller supplies that per operation.
//
// Plain HTTP is enabled for a loopback host so a local test registry needs no
// certificate. The rule is deliberately narrow: anything that is not loopback
// keeps TLS, because a silent downgrade would send a token in the clear.
func newRepository(base, repo string) (*remote.Repository, error) {
	ref := TrimBaseScheme(base) + "/" + repo
	r, err := remote.NewRepository(ref)
	if err != nil {
		return nil, fmt.Errorf("invalid registry reference %q: %w", ref, err)
	}
	r.Client = newAuthClient()
	r.PlainHTTP = isLoopback(r.Reference.Registry)
	return r, nil
}

// isLoopback reports whether host is this machine, with or without a port.
//
// Parsed rather than string-matched: an IPv6 host is bracketed and full of
// colons, so splitting on the first one turns "[::1]:5000" into nothing, and
// "127.0.0.1.example.test" must not pass a prefix test.
func isLoopback(host string) bool {
	if h, _, err := net.SplitHostPort(host); err == nil {
		host = h
	}
	host = strings.Trim(host, "[]")
	if host == "localhost" {
		return true
	}
	ip := net.ParseIP(host)
	return ip != nil && ip.IsLoopback()
}
