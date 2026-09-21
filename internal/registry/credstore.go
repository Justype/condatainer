package registry

import (
	"encoding/base64"
	"encoding/json"
	"errors"
	"fmt"
	"io/fs"
	"os"
	"path/filepath"
	"strings"

	"oras.land/oras-go/v2/registry/remote/auth"

	"github.com/Justype/condatainer/internal/config"
)

// credentialFileName is the credential file in each config layer's directory.
const credentialFileName = "registry-auth.json"

// credentialLayerNames are the layers that can hold credentials, nearest first.
var credentialLayerNames = []string{"user", "extra-root", "app-root", "system"}

// credentialFilePath is where a config layer keeps its credentials. Tests
// replace it.
var credentialFilePath = func(layer string) (string, error) {
	path, err := config.GetConfigPathByLayer(layer)
	if err != nil {
		return "", err
	}
	return filepath.Join(filepath.Dir(path), credentialFileName), nil
}

// authFile is one layer's credentials, keyed by a registry host or by
// "host/owner/repo".
type authFile struct {
	Auths map[string]authEntry `json:"auths"`
}

// authEntry is Docker's config shape: "auth" is base64 of "user:password".
type authEntry struct {
	Auth string `json:"auth"`
}

func (e authEntry) credential() auth.Credential {
	raw, err := base64.StdEncoding.DecodeString(e.Auth)
	if err != nil {
		return auth.EmptyCredential
	}
	user, password, _ := strings.Cut(string(raw), ":")
	return auth.Credential{Username: user, Password: password}
}

func entryFor(cred auth.Credential) authEntry {
	return authEntry{Auth: base64.StdEncoding.EncodeToString([]byte(cred.Username + ":" + cred.Password))}
}

// readAuthFile reads a credential file; a missing file holds nothing.
func readAuthFile(path string) (authFile, error) {
	data, err := os.ReadFile(path)
	if errors.Is(err, fs.ErrNotExist) {
		return authFile{}, nil
	}
	if err != nil {
		return authFile{}, fmt.Errorf("cannot read %s: %w", path, err)
	}
	var file authFile
	if err := json.Unmarshal(data, &file); err != nil {
		return authFile{}, fmt.Errorf("cannot parse %s: %w", path, err)
	}
	return file, nil
}

// writeAuthFile replaces path atomically with mode 0600. The directory is not
// created: a shared layer's directory belongs to whoever set it up.
func writeAuthFile(path string, file authFile) error {
	data, err := json.MarshalIndent(file, "", "  ")
	if err != nil {
		return err
	}
	tmp, err := os.CreateTemp(filepath.Dir(path), ".registry-auth-*")
	if err != nil {
		return fmt.Errorf("cannot write %s: %w", path, err)
	}
	defer os.Remove(tmp.Name())
	if _, err := tmp.Write(append(data, '\n')); err != nil {
		tmp.Close()
		return err
	}
	if err := tmp.Close(); err != nil {
		return err
	}
	return os.Rename(tmp.Name(), path)
}

// splitTarget separates a login target into its host and "host/owner/repo" key.
// It accepts a bare host or a host with a repository path, nothing else.
func splitTarget(target string) (key string, err error) {
	key = TrimBaseScheme(target)
	host, _, _ := strings.Cut(key, "/")
	if host == "" || strings.Contains(key, "//") || strings.ContainsAny(key, " \t?#@") {
		return "", fmt.Errorf("%q is not a registry host or host/owner/repository", target)
	}
	return key, nil
}
