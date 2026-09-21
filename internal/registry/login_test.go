package registry

import (
	"os"
	"path/filepath"
	"reflect"
	"strings"
	"testing"
)

func TestLoginAndLogoutRejectATargetThatIsNotAHostOrRepository(t *testing.T) {
	useLayers(t, map[string]authFile{"user": {}})
	for _, target := range []string{"", "ghcr.io//lab", "ghcr.io/my lab", "https://ghcr.io/lab@x"} {
		if err := Login(t.Context(), target, "user", "user", "token"); err == nil {
			t.Errorf("Login(%q) accepted it", target)
		}
		if err := Logout(t.Context(), target, "user"); err == nil {
			t.Errorf("Logout(%q) accepted it", target)
		}
	}
}

func TestStoredCredentialsListsKeysAndLayersWithoutSecrets(t *testing.T) {
	paths := useLayers(t, map[string]authFile{
		"user":       {Auths: map[string]authEntry{"ghcr.io/my-lab/rnaseq": entry("bob")}},
		"extra-root": {Auths: map[string]authEntry{"ghcr.io": entry("group")}},
	})
	got, err := StoredCredentials(t.Context())
	if err != nil {
		t.Fatal(err)
	}
	want := []StoredCredential{
		{Key: "ghcr.io", Layer: "extra-root", Username: "group", Path: paths["extra-root"], ReadableBy: "you"},
		{Key: "ghcr.io/my-lab/rnaseq", Layer: "user", Username: "bob", Path: paths["user"], ReadableBy: "you"},
	}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("got %+v, want %+v", got, want)
	}
}

func TestLogoutRemovesOnlyTheNamedKey(t *testing.T) {
	paths := useLayers(t, map[string]authFile{"user": {Auths: map[string]authEntry{
		"ghcr.io": entry("host"), "ghcr.io/my-lab/rnaseq": entry("repo"),
	}}})
	if err := Logout(t.Context(), "ghcr.io/my-lab/rnaseq", "user"); err != nil {
		t.Fatal(err)
	}
	file, err := readAuthFile(paths["user"])
	if err != nil {
		t.Fatal(err)
	}
	if _, ok := file.Auths["ghcr.io"]; !ok || len(file.Auths) != 1 {
		t.Errorf("auths = %v, want only the host entry", file.Auths)
	}
	if err := Logout(t.Context(), "ghcr.io/my-lab/rnaseq", "user"); err == nil {
		t.Error("removing a missing credential succeeded")
	}
}

// A shared layer's directory is not ours to create.
func TestLayerFileRefusesAMissingSharedDirectory(t *testing.T) {
	missing := filepath.Join(t.TempDir(), "nope", credentialFileName)
	prev := credentialFilePath
	credentialFilePath = func(string) (string, error) { return missing, nil }
	t.Cleanup(func() { credentialFilePath = prev })

	if _, err := layerFile("extra-root"); err == nil || !strings.Contains(err.Error(), "does not exist") {
		t.Fatalf("extra-root: %v", err)
	}
	if _, err := layerFile("user"); err != nil {
		t.Fatalf("user: %v", err)
	}
	if info, err := os.Stat(filepath.Dir(missing)); err != nil || info.Mode().Perm() != 0o700 {
		t.Errorf("the user's directory was not created private: %v", err)
	}
}

func TestWriteAuthFileIsPrivate(t *testing.T) {
	path := filepath.Join(t.TempDir(), credentialFileName)
	if err := writeAuthFile(path, authFile{Auths: map[string]authEntry{"ghcr.io": entry("u")}}); err != nil {
		t.Fatal(err)
	}
	if info, err := os.Stat(path); err != nil || info.Mode().Perm() != 0o600 {
		t.Fatalf("mode = %v, %v", info.Mode(), err)
	}
}
