package registry

import (
	"cmp"
	"context"
	"net/http"
	"os"
	"strings"

	"oras.land/oras-go/v2/registry/remote/auth"
	"oras.land/oras-go/v2/registry/remote/credentials"
	"oras.land/oras-go/v2/registry/remote/retry"
)

// Environment credentials, for a host with no usable Docker credential store.
const (
	// EnvToken and EnvUser are checked for every host.
	EnvToken = "CNT_REGISTRY_TOKEN"
	EnvUser  = "CNT_REGISTRY_USER"
	// EnvGitHubToken is GitHub's own name and applies to ghcr.io alone. It is not
	// ours to rename, and it is already set in every GitHub Actions job.
	EnvGitHubToken = "GITHUB_TOKEN"

	// ghcrHost is the one registry EnvGitHubToken speaks for.
	ghcrHost = "ghcr.io"
	// defaultTokenUser is what a registry expects beside a token when the token
	// itself carries the identity. GHCR accepts any username with a valid PAT.
	defaultTokenUser = "x-access-token"
)

// credentialFunc resolves credentials per registry host, in order: an explicit
// environment token, the Docker/OCI credential store, then anonymous.
//
// The environment path exists for compute nodes, where there is no interactive
// login and $HOME may not be the one holding ~/.docker/config.json. Anonymous is
// the normal case for public artifacts, so a missing store is not an error.
func credentialFunc() auth.CredentialFunc {
	store, storeErr := credentials.NewStoreFromDocker(credentials.StoreOptions{})
	return func(ctx context.Context, host string) (auth.Credential, error) {
		if cred, ok := envCredential(host); ok {
			return cred, nil
		}
		if storeErr != nil {
			return auth.EmptyCredential, nil
		}
		return credentials.Credential(store)(ctx, host)
	}
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
	cred, err := credentialFunc()(ctx, host)
	if err != nil {
		return false
	}
	return cred != auth.EmptyCredential
}

// envCredential builds a credential from the environment for host, reporting
// false when nothing applies. EnvToken wins for any host; EnvGitHubToken is
// consulted only for ghcr.io, so a job's GitHub token is never sent elsewhere.
func envCredential(host string) (auth.Credential, bool) {
	token := strings.TrimSpace(os.Getenv(EnvToken))
	if token == "" && strings.EqualFold(host, ghcrHost) {
		token = strings.TrimSpace(os.Getenv(EnvGitHubToken))
	}
	if token == "" {
		return auth.Credential{}, false
	}
	user := strings.TrimSpace(os.Getenv(EnvUser))
	if user == "" {
		user = defaultTokenUser
	}
	return auth.Credential{Username: user, Password: token}, true
}

// newAuthClient is the client every registry operation uses: retrying transport,
// per-host token cache, and the credential chain above.
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
func newAuthClient() *auth.Client {
	inner := *retry.DefaultClient
	inner.Transport = &inspectTransport{next: cmp.Or(inner.Transport, http.DefaultTransport)}

	client := &auth.Client{
		Client:     &inner,
		Cache:      auth.NewCache(),
		Credential: credentialFunc(),
	}
	client.SetUserAgent(userAgent())
	return client
}
