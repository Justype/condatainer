package registry

import "testing"

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

// A missing or unreadable Docker credential store means anonymous access, which
// is the normal case for a public artifact — never an error.
func TestCredentialFuncFallsBackToAnonymous(t *testing.T) {
	t.Setenv(EnvGitHubToken, "")
	t.Setenv("DOCKER_CONFIG", t.TempDir()) // exists, holds no config.json

	cred, err := credentialFunc()(t.Context(), "registry.example.test")
	if err != nil {
		t.Fatalf("credentialFunc: %v", err)
	}
	if cred.Username != "" || cred.Password != "" || cred.RefreshToken != "" || cred.AccessToken != "" {
		t.Errorf("credential = %+v, want empty", cred)
	}
}
