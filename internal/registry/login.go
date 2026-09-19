package registry

import (
	"context"
	"encoding/base64"
	"encoding/json"
	"errors"
	"fmt"
	"io/fs"
	"os"
	"sort"
	"strings"

	"oras.land/oras-go/v2/registry/remote"
	"oras.land/oras-go/v2/registry/remote/auth"
	"oras.land/oras-go/v2/registry/remote/credentials"
	"oras.land/oras-go/v2/registry/remote/retry"
)

// Login verifies a credential against registry and stores it in the Docker/OCI
// credential store used by every other operation in this package.
func Login(ctx context.Context, registry, username, password string) error {
	registry = TrimBaseScheme(registry)
	if registry == "" || strings.Contains(registry, "/") {
		return fmt.Errorf("registry login needs a host, got %q", registry)
	}
	if password == "" {
		return fmt.Errorf("registry password or token is empty")
	}

	store, err := credentials.NewStoreFromDocker(credentials.StoreOptions{AllowPlaintextPut: true})
	if err != nil {
		return fmt.Errorf("cannot open registry credential store: %w", err)
	}
	reg, err := remote.NewRegistry(registry)
	if err != nil {
		return fmt.Errorf("invalid registry %q: %w", registry, err)
	}
	reg.Client = &auth.Client{Client: retry.DefaultClient, Cache: auth.NewCache()}
	reg.PlainHTTP = isLoopback(reg.Reference.Registry)

	cred := auth.Credential{Username: username, Password: password}
	if username == "" {
		cred = auth.Credential{Password: password, AccessToken: password}
	}
	if err := credentials.Login(ctx, store, reg, cred); err != nil {
		return fmt.Errorf("login to %s failed: %w", registry, classify(err))
	}
	return nil
}

// Logout removes the stored credential for registry.
func Logout(ctx context.Context, registry string) error {
	registry = TrimBaseScheme(registry)
	if registry == "" || strings.Contains(registry, "/") {
		return fmt.Errorf("registry logout needs a host, got %q", registry)
	}
	store, err := credentials.NewStoreFromDocker(credentials.StoreOptions{})
	if err != nil {
		return fmt.Errorf("cannot open registry credential store: %w", err)
	}
	if err := credentials.Logout(ctx, store, registry); err != nil {
		return fmt.Errorf("logout from %s failed: %w", registry, err)
	}
	return nil
}

// StoredCredential is one host with a saved credential. The secret is never
// read into it.
type StoredCredential struct {
	Host     string
	Username string // empty when a helper holds the credential or it is a bare token
	Helper   string // credential helper serving this host, empty for the config file
}

// StoredCredentials lists the hosts the Docker/OCI credential store has
// credentials for, sorted by host. A host held only by a global credsStore
// helper is not listed: the helper cannot be asked which hosts it knows.
func StoredCredentials(ctx context.Context) ([]StoredCredential, error) {
	store, err := credentials.NewStoreFromDocker(credentials.StoreOptions{})
	if err != nil {
		return nil, fmt.Errorf("cannot open registry credential store: %w", err)
	}
	data, err := os.ReadFile(store.ConfigPath())
	if errors.Is(err, fs.ErrNotExist) {
		return nil, nil
	}
	if err != nil {
		return nil, fmt.Errorf("cannot read %s: %w", store.ConfigPath(), err)
	}
	var cfg struct {
		Auths map[string]struct {
			Auth     string `json:"auth"`
			Username string `json:"username"`
		} `json:"auths"`
		CredHelpers map[string]string `json:"credHelpers"`
	}
	if err := json.Unmarshal(data, &cfg); err != nil {
		return nil, fmt.Errorf("cannot parse %s: %w", store.ConfigPath(), err)
	}

	byHost := map[string]StoredCredential{}
	for host, entry := range cfg.Auths {
		user := entry.Username
		if raw, err := base64.StdEncoding.DecodeString(entry.Auth); err == nil && user == "" {
			user, _, _ = strings.Cut(string(raw), ":")
		}
		byHost[host] = StoredCredential{Host: host, Username: user}
	}
	for host, helper := range cfg.CredHelpers {
		byHost[host] = StoredCredential{Host: host, Helper: helper}
	}
	out := make([]StoredCredential, 0, len(byHost))
	for _, c := range byHost {
		out = append(out, c)
	}
	sort.Slice(out, func(i, j int) bool { return out[i].Host < out[j].Host })
	return out, nil
}
