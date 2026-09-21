package registry

import (
	"errors"
	"path/filepath"
	"testing"

	"oras.land/oras-go/v2/registry/remote/auth"
)

func TestEnvCredential(t *testing.T) {
	tests := []struct {
		name  string
		token string
		host  string
		want  bool
	}{
		{"no token is not a credential", "", "ghcr.io", false},
		{"GITHUB_TOKEN speaks for ghcr.io", "gh", "ghcr.io", true},
		{"the host matches case-insensitively", "gh", "GHCR.IO", true},
		// A CI job's GitHub token must never reach an unrelated registry.
		{"GITHUB_TOKEN goes nowhere else", "gh", "registry.example.test", false},
		{"whitespace is not a token", "   ", "ghcr.io", false},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			t.Setenv(EnvGitHubToken, tt.token)
			cred, ok := envCredential(tt.host)
			if ok != tt.want {
				t.Fatalf("ok = %v, want %v", ok, tt.want)
			}
			if ok && (cred.Username != defaultTokenUser || cred.Password != tt.token) {
				t.Errorf("credential = %+v", cred)
			}
		})
	}
}

// useLayers points each layer's credential file at a temp directory holding the
// given entries, and leaves the other layers unavailable.
func useLayers(t *testing.T, layers map[string]authFile) map[string]string {
	t.Helper()
	paths := map[string]string{}
	for layer, file := range layers {
		dir := t.TempDir()
		paths[layer] = filepath.Join(dir, credentialFileName)
		if err := writeAuthFile(paths[layer], file); err != nil {
			t.Fatal(err)
		}
	}
	prev := credentialFilePath
	credentialFilePath = func(layer string) (string, error) {
		if path, ok := paths[layer]; ok {
			return path, nil
		}
		return "", errors.New("layer not available")
	}
	t.Cleanup(func() { credentialFilePath = prev })
	return paths
}

func entry(user string) authEntry {
	return entryFor(auth.Credential{Username: user, Password: "secret"})
}

// A missing store means anonymous access, which is the normal case for a public
// artifact — never an error.
func TestCredentialFuncFallsBackToAnonymous(t *testing.T) {
	t.Setenv(EnvGitHubToken, "")
	useLayers(t, nil)

	cred, err := credentialFunc("registry.example.test/lab/repo")(t.Context(), "registry.example.test")
	if err != nil {
		t.Fatalf("credentialFunc: %v", err)
	}
	if cred != auth.EmptyCredential {
		t.Errorf("credential = %+v, want empty", cred)
	}
}

// The most specific key wins over a nearer layer, and among layers holding the
// same key the nearest wins.
func TestCredentialFuncPrefersTheMostSpecificKeyThenTheNearestLayer(t *testing.T) {
	t.Setenv(EnvGitHubToken, "")
	useLayers(t, map[string]authFile{
		"user": {Auths: map[string]authEntry{"ghcr.io": entry("personal-host")}},
		"extra-root": {Auths: map[string]authEntry{
			"ghcr.io":               entry("group-host"),
			"ghcr.io/my-lab/rnaseq": entry("group-repo"),
			"ghcr.io/my-lab":        entry("group-owner"),
		}},
	})
	tests := []struct{ scope, want string }{
		{"ghcr.io/my-lab/rnaseq/cnt", "group-repo"},
		{"ghcr.io/my-lab/other/cnt", "group-owner"},
		{"ghcr.io/someone/else", "personal-host"},
		{"ghcr.io", "personal-host"},
	}
	for _, tt := range tests {
		cred, err := credentialFunc(tt.scope)(t.Context(), "ghcr.io")
		if err != nil {
			t.Fatal(err)
		}
		if cred.Username != tt.want {
			t.Errorf("scope %s: user = %q, want %q", tt.scope, cred.Username, tt.want)
		}
	}
}

// A repository's credential is never sent to another host, and GITHUB_TOKEN
// still comes first.
func TestCredentialFuncKeepsHostsApartAndEnvironmentFirst(t *testing.T) {
	useLayers(t, map[string]authFile{
		"user": {Auths: map[string]authEntry{"ghcr.io/my-lab/rnaseq": entry("repo")}},
	})
	t.Setenv(EnvGitHubToken, "")
	cred, _ := credentialFunc("ghcr.io/my-lab/rnaseq/cnt")(t.Context(), "registry.example.test")
	if cred != auth.EmptyCredential {
		t.Errorf("another host received %+v", cred)
	}
	t.Setenv(EnvGitHubToken, "env-token")
	cred, _ = credentialFunc("ghcr.io/my-lab/rnaseq/cnt")(t.Context(), "ghcr.io")
	if cred.Password != "env-token" {
		t.Errorf("GITHUB_TOKEN did not win: %+v", cred)
	}
}
