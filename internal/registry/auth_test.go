package registry

import "testing"

// The environment path exists for hosts with no usable credential store, so its
// precedence and its blast radius are both load-bearing.
func TestEnvCredential(t *testing.T) {
	tests := []struct {
		name       string
		token      string
		user       string
		githubTok  string
		host       string
		wantOK     bool
		wantUser   string
		wantSecret string
	}{
		{
			name: "no environment token is not a credential",
			host: "ghcr.io",
		},
		{
			name:  "the explicit token applies to any host",
			token: "t0k", host: "registry.example.test",
			wantOK: true, wantUser: defaultTokenUser, wantSecret: "t0k",
		},
		{
			name:  "an explicit user overrides the default",
			token: "t0k", user: "alice", host: "ghcr.io",
			wantOK: true, wantUser: "alice", wantSecret: "t0k",
		},
		{
			name:      "GITHUB_TOKEN speaks for ghcr.io",
			githubTok: "gh", host: "ghcr.io",
			wantOK: true, wantUser: defaultTokenUser, wantSecret: "gh",
		},
		{
			name:      "GITHUB_TOKEN is matched case-insensitively",
			githubTok: "gh", host: "GHCR.IO",
			wantOK: true, wantUser: defaultTokenUser, wantSecret: "gh",
		},
		{
			// A CI job's GitHub token must never reach an unrelated registry.
			name:      "GITHUB_TOKEN goes nowhere else",
			githubTok: "gh", host: "registry.example.test",
		},
		{
			name:  "the explicit token wins over GITHUB_TOKEN on ghcr.io",
			token: "t0k", githubTok: "gh", host: "ghcr.io",
			wantOK: true, wantUser: defaultTokenUser, wantSecret: "t0k",
		},
		{
			name:  "whitespace is not a token",
			token: "   ", host: "ghcr.io",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			t.Setenv(EnvToken, tt.token)
			t.Setenv(EnvUser, tt.user)
			t.Setenv(EnvGitHubToken, tt.githubTok)

			cred, ok := envCredential(tt.host)
			if ok != tt.wantOK {
				t.Fatalf("ok = %v, want %v (cred = %+v)", ok, tt.wantOK, cred)
			}
			if !ok {
				return
			}
			if cred.Username != tt.wantUser || cred.Password != tt.wantSecret {
				t.Errorf("credential = %s/%s, want %s/%s",
					cred.Username, cred.Password, tt.wantUser, tt.wantSecret)
			}
		})
	}
}

// A missing or unreadable Docker credential store means anonymous access, which
// is the normal case for a public artifact — never an error.
func TestCredentialFuncFallsBackToAnonymous(t *testing.T) {
	t.Setenv(EnvToken, "")
	t.Setenv(EnvUser, "")
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

// An environment token must be used even when a credential store is readable:
// on a compute node the store is usually the stale half of the pair.
func TestCredentialFuncPrefersTheEnvironmentToken(t *testing.T) {
	t.Setenv(EnvToken, "t0k")
	t.Setenv(EnvUser, "")
	t.Setenv("DOCKER_CONFIG", t.TempDir())

	cred, err := credentialFunc()(t.Context(), "registry.example.test")
	if err != nil {
		t.Fatalf("credentialFunc: %v", err)
	}
	if cred.Password != "t0k" || cred.Username != defaultTokenUser {
		t.Errorf("credential = %+v, want the environment token", cred)
	}
}
