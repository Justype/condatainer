package registry

import (
	"cmp"
	"context"
	"net/http"
	"os"
	"strings"

	"oras.land/oras-go/v2/registry/remote/auth"
	"oras.land/oras-go/v2/registry/remote/retry"
)

// Environment credential, for a host with no usable Docker credential store.
const (
	// EnvGitHubToken is GitHub's own name and applies to ghcr.io alone. It is not
	// ours to rename, and it is already set in every GitHub Actions job.
	EnvGitHubToken = "GITHUB_TOKEN"

	// ghcrHost is the one registry EnvGitHubToken speaks for.
	ghcrHost = "ghcr.io"
	// defaultTokenUser is what a registry expects beside a token when the token
	// itself carries the identity. GHCR accepts any username with a valid PAT.
	defaultTokenUser = "x-access-token"
)

// credentialFunc resolves credentials for the repository at scope
// ("host/owner/repo", no tag), in order: GITHUB_TOKEN for ghcr.io, then the
// stored credentials, then anonymous. Anonymous is the normal case for public
// artifacts, so an unreadable store is not an error.
//
// A stored credential is the first hit going from the most specific key (the
// repository) to the general one (the host), trying the layers nearest first
// for each key.
func credentialFunc(scope string) auth.CredentialFunc {
	var files []authFile
	loaded := false
	return func(ctx context.Context, host string) (auth.Credential, error) {
		if cred, ok := envCredential(host); ok {
			return cred, nil
		}
		if !loaded {
			files, loaded = readLayers(), true
		}
		for _, key := range storeKeys(scope, host) {
			for _, file := range files {
				if entry, ok := file.Auths[key]; ok {
					if cred := entry.credential(); cred != auth.EmptyCredential {
						return cred, nil
					}
				}
			}
		}
		return auth.EmptyCredential, nil
	}
}

// readLayers reads every layer's credential file, nearest first. A layer that is
// unavailable or unreadable holds nothing.
func readLayers() []authFile {
	var files []authFile
	for _, layer := range credentialLayerNames {
		path, err := credentialFilePath(layer)
		if err != nil {
			continue
		}
		if file, err := readAuthFile(path); err == nil {
			files = append(files, file)
		}
	}
	return files
}

// storeKeys lists the keys that may hold a credential for host, most specific
// first: scope itself, each parent path, then the bare host. A scope on another
// host contributes only that host's key.
func storeKeys(scope, host string) []string {
	scopeHost, _, _ := strings.Cut(scope, "/")
	if !strings.EqualFold(scopeHost, host) {
		return []string{host}
	}
	var keys []string
	for path := scope; strings.Contains(path, "/"); path = path[:strings.LastIndex(path, "/")] {
		keys = append(keys, path)
	}
	return append(keys, scopeHost)
}

// HasCredential reports whether any credential is available for the registry
// named by ref, which may be a bare host or a host/prefix coordinate.
//
// It separates the two refusals that arrive as one status: a credential that
// does not open an artifact is the expired-token case ErrUnauthorized exists to
// report, while a refusal with nothing in hand is a closed door — an endpoint
// this user cannot pull from at all, which for a caller is the same outcome as
// the artifact not being published.
func HasCredential(ctx context.Context, ref string) bool {
	host, _, _ := strings.Cut(TrimBaseScheme(ref), "/")
	if host == "" {
		return false
	}
	cred, err := credentialFunc(TrimBaseScheme(ref))(ctx, host)
	if err != nil {
		return false
	}
	return cred != auth.EmptyCredential
}

// envCredential builds a credential from GITHUB_TOKEN for ghcr.io, reporting
// false for any other host so a job's token is never sent elsewhere.
func envCredential(host string) (auth.Credential, bool) {
	if !strings.EqualFold(host, ghcrHost) {
		return auth.Credential{}, false
	}
	token := strings.TrimSpace(os.Getenv(EnvGitHubToken))
	if token == "" {
		return auth.Credential{}, false
	}
	return auth.Credential{Username: defaultTokenUser, Password: token}, true
}

// newAuthClient is the client every registry operation on the repository at scope
// uses: retrying transport, per-host token cache, and the credential chain above.
//
// The retrying transport stays for the small requests — token exchange, HEAD,
// manifest reads — and is irrelevant to a rate-limited blob: its waits are
// milliseconds where a secondary limit needs minutes, it does not recognize
// GitHub's 403, and it cannot replay a body without GetBody. That case belongs
// to [retryPolicy]. It sets no client-level timeout, which must stay true: one
// would kill a long upload outright.
//
// [inspectTransport] wraps it to keep the response metadata ORAS discards and to
// ask permission before sending a large body.
func newAuthClient(scope string) *auth.Client {
	inner := *retry.DefaultClient
	inner.Transport = &inspectTransport{next: cmp.Or(inner.Transport, http.DefaultTransport)}

	client := &auth.Client{
		Client:     &inner,
		Cache:      auth.NewCache(),
		Credential: credentialFunc(scope),
	}
	client.SetUserAgent(userAgent())
	return client
}
