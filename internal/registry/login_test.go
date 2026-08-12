package registry

import "testing"

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
