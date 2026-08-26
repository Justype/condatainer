package publish

import "testing"

func TestNormalizeRemoteURL(t *testing.T) {
	// The point is that these are not cosmetic variants: SSH and HTTPS clones of
	// one repository leave different strings in different collaborators'
	// checkouts, and publishing whichever the pusher happened to have would put
	// two annotations on one project's packages.
	want := "https://github.com/my-lab/rnaseq-2026"
	for _, spelling := range []string{
		"git@github.com:my-lab/rnaseq-2026.git",
		"git@github.com:my-lab/rnaseq-2026",
		"https://github.com/my-lab/rnaseq-2026.git",
		"https://github.com/my-lab/rnaseq-2026",
		"https://github.com/my-lab/rnaseq-2026/",
		"ssh://git@github.com/my-lab/rnaseq-2026.git",
		"  https://github.com/my-lab/rnaseq-2026\n",
	} {
		if got := NormalizeRemoteURL(spelling); got != want {
			t.Errorf("NormalizeRemoteURL(%q) = %q, want %q", spelling, got, want)
		}
	}

	// Any host, and a nested group is a repository path rather than a page.
	for _, tc := range []struct{ raw, want string }{
		{"git@gitlab.com:my-lab/p.git", "https://gitlab.com/my-lab/p"},
		{"https://gitlab.com/my-lab/p", "https://gitlab.com/my-lab/p"},
		{"git@gitlab.company.example:group/sub/p.git", "https://gitlab.company.example/group/sub/p"},
		{"ssh://git@git.institute.example/lab/p.git", "https://git.institute.example/lab/p"},
		{"http://gitea.lab.example/team/p.git", "https://gitea.lab.example/team/p"},
	} {
		if got := NormalizeRemoteURL(tc.raw); got != tc.want {
			t.Errorf("NormalizeRemoteURL(%q) = %q, want %q", tc.raw, got, tc.want)
		}
	}

	// Names no browsable repository, so nothing is recorded rather than a URL
	// nothing serves.
	for _, other := range []string{
		"",
		"https://github.com/my-lab", // no repository
		"https://github.com/",       // no path
		"https://github.com/a/b/c",  // a page inside a repository
		"git@server:/srv/git/p.git", // absolute path: a bare repository
		"/srv/repos/p.git",          // a local path, not a remote
		"file:///srv/repos/p.git",   // ditto
		"https://example.invalid",   // host only
	} {
		if got := NormalizeRemoteURL(other); got != "" {
			t.Errorf("NormalizeRemoteURL(%q) = %q, want \"\"", other, got)
		}
	}
}
