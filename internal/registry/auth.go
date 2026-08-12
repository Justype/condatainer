package registry

import (
	"context"
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
func newAuthClient() *auth.Client {
	return &auth.Client{
		Client:     retry.DefaultClient,
		Cache:      auth.NewCache(),
		Credential: credentialFunc(),
	}
}
