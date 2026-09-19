package registry

import (
	"encoding/base64"
	"os"
	"path/filepath"
	"reflect"
	"testing"
)

func TestLoginAndLogoutRequireARegistryHost(t *testing.T) {
	for _, registry := range []string{"", "ghcr.io/owner"} {
		if err := Login(t.Context(), registry, "user", "token"); err == nil {
			t.Errorf("Login(%q) accepted a non-host", registry)
		}
		if err := Logout(t.Context(), registry); err == nil {
			t.Errorf("Logout(%q) accepted a non-host", registry)
		}
	}
}

func TestStoredCredentialsListsHostsWithoutSecrets(t *testing.T) {
	dir := t.TempDir()
	t.Setenv("DOCKER_CONFIG", dir)
	config := `{"auths":{"ghcr.io":{"auth":"` + base64.StdEncoding.EncodeToString([]byte("alice:s3cret")) + `"}},
"credHelpers":{"registry.lab.example":"pass"}}`
	if err := os.WriteFile(filepath.Join(dir, "config.json"), []byte(config), 0o600); err != nil {
		t.Fatal(err)
	}
	got, err := StoredCredentials(t.Context())
	if err != nil {
		t.Fatal(err)
	}
	want := []StoredCredential{
		{Host: "ghcr.io", Username: "alice"},
		{Host: "registry.lab.example", Helper: "pass"},
	}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("got %+v, want %+v", got, want)
	}
}
