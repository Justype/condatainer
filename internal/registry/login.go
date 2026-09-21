package registry

import (
	"context"
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"oras.land/oras-go/v2/registry/remote"
	"oras.land/oras-go/v2/registry/remote/auth"
	"oras.land/oras-go/v2/registry/remote/retry"
)

// Login verifies a credential against the registry host in target and stores it
// in one config layer, under target itself: a bare host, or host/owner/repo to
// cover one repository. The layer is "user", "extra-root", "app-root" or
// "system". The layer's directory must already exist unless it is the user's.
func Login(ctx context.Context, target, layer, username, password string) error {
	key, err := splitTarget(target)
	if err != nil {
		return err
	}
	if password == "" {
		return fmt.Errorf("registry password or token is empty")
	}
	path, err := layerFile(layer)
	if err != nil {
		return err
	}

	host, _, _ := strings.Cut(key, "/")
	reg, err := remote.NewRegistry(host)
	if err != nil {
		return fmt.Errorf("invalid registry %q: %w", host, err)
	}
	cred := auth.Credential{Username: username, Password: password}
	reg.Client = &auth.Client{
		Client:     retry.DefaultClient,
		Cache:      auth.NewCache(),
		Credential: auth.StaticCredential(reg.Reference.Registry, cred),
	}
	reg.PlainHTTP = isLoopback(reg.Reference.Registry)
	if err := reg.Ping(ctx); err != nil {
		return fmt.Errorf("login to %s failed: %w", key, classify(err))
	}

	file, err := readAuthFile(path)
	if err != nil {
		return err
	}
	if file.Auths == nil {
		file.Auths = map[string]authEntry{}
	}
	file.Auths[key] = entryFor(cred)
	return writeAuthFile(path, file)
}

// Logout removes the credential stored under target in one config layer.
func Logout(ctx context.Context, target, layer string) error {
	key, err := splitTarget(target)
	if err != nil {
		return err
	}
	path, err := layerFile(layer)
	if err != nil {
		return err
	}
	file, err := readAuthFile(path)
	if err != nil {
		return err
	}
	if _, ok := file.Auths[key]; !ok {
		return fmt.Errorf("no stored credential for %s in the %s layer", key, layer)
	}
	delete(file.Auths, key)
	return writeAuthFile(path, file)
}

// layerFile is the credential file for a layer, refusing a directory that does
// not exist unless it is the user's: a shared layer's directory is not ours to
// create.
func layerFile(layer string) (string, error) {
	path, err := credentialFilePath(layer)
	if err != nil {
		return "", err
	}
	dir := filepath.Dir(path)
	if layer == "user" {
		return path, os.MkdirAll(dir, 0o700)
	}
	if info, err := os.Stat(dir); err != nil || !info.IsDir() {
		return "", fmt.Errorf("%s does not exist", dir)
	}
	return path, nil
}

// StoredCredential is one saved credential. The secret is never read into it.
type StoredCredential struct {
	Key      string // a registry host, or host/owner/repo
	Layer    string
	Username string
	Path     string // the file that holds it

	// ReadableBy is who can read that file: "you", "group <name>" or "everyone".
	ReadableBy string
}

// StoredCredentials lists every saved credential across the layers, by key and
// then nearest layer first.
func StoredCredentials(ctx context.Context) ([]StoredCredential, error) {
	var out []StoredCredential
	for _, layer := range credentialLayerNames {
		path, err := credentialFilePath(layer)
		if err != nil {
			continue
		}
		file, err := readAuthFile(path)
		if err != nil {
			return nil, err
		}
		for key, entry := range file.Auths {
			out = append(out, StoredCredential{
				Key: key, Layer: layer, Username: entry.credential().Username, Path: path,
				ReadableBy: ReadableBy(path),
			})
		}
	}
	rank := map[string]int{}
	for i, layer := range credentialLayerNames {
		rank[layer] = i
	}
	sort.Slice(out, func(i, j int) bool {
		if out[i].Key != out[j].Key {
			return out[i].Key < out[j].Key
		}
		return rank[out[i].Layer] < rank[out[j].Layer]
	})
	return out, nil
}
