package registry

import (
	"context"
	"fmt"
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
