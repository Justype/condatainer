package conda

import (
	"strings"
	"testing"
)

// The property that matters: the file written from a solve and the file written
// from the installed environment's export are the same bytes, because the
// identity hashes them and one artifact may not have two identities.
func TestExplicitFromMatchesTheExportPath(t *testing.T) {
	packages := []Package{
		{Name: "python", Version: "3.10.14", URL: "https://conda.anaconda.org/conda-forge/linux-64/python-3.10.14-hd12c33a_0.tar.bz2", Channel: "conda-forge"},
		{Name: "cutadapt", Version: "5.0", URL: "https://conda.anaconda.org/bioconda/linux-64/cutadapt-5.0-py310h1fe012e_0.tar.bz2", Channel: "bioconda"},
	}
	fromSolve, err := ExplicitFrom(packages)
	if err != nil {
		t.Fatal(err)
	}
	// What `micromamba env export --explicit` emits: commentary, the marker, then
	// the URLs in whatever order it likes.
	export := "# This file may be used to create an environment using:\n" +
		"# platform: linux-64\n@EXPLICIT\n" +
		packages[1].URL + "\n" + packages[0].URL + "\n"
	fromExport, err := CanonicalExplicit([]byte(export))
	if err != nil {
		t.Fatal(err)
	}
	if string(fromSolve) != string(fromExport) {
		t.Fatalf("solve and export disagree:\nsolve:\n%s\nexport:\n%s", fromSolve, fromExport)
	}
}

func TestEnvironmentFromMatchesTheExportPath(t *testing.T) {
	packages := []Package{
		{Name: "python", Version: "3.10.14", Channel: "conda-forge", URL: "https://conda.anaconda.org/conda-forge/linux-64/python-3.10.14-hd12c33a_0.tar.bz2"},
		{Name: "cutadapt", Version: "5.0", Channel: "bioconda", URL: "https://conda.anaconda.org/bioconda/linux-64/cutadapt-5.0-py310h1fe012e_0.tar.bz2"},
	}
	priority := []string{"conda-forge", "bioconda"}
	fromSolve, err := EnvironmentFrom(packages, priority)
	if err != nil {
		t.Fatal(err)
	}
	export := "name: scratch\nchannels:\n  - bioconda\n  - conda-forge\ndependencies:\n" +
		"  - cutadapt=5.0\n  - python=3.10.14\nprefix: /cnt/tmp\n"
	fromExport, err := CanonicalEnvironment([]byte(export), priority)
	if err != nil {
		t.Fatal(err)
	}
	if string(fromSolve) != string(fromExport) {
		t.Fatalf("solve and export disagree:\nsolve:\n%s\nexport:\n%s", fromSolve, fromExport)
	}
}

// A solver reports a channel as a URL, with the subdir attached; an export names
// the channel alone. They have to reduce to the same thing or the key moves.
func TestChannelNameReducesSolverSpellings(t *testing.T) {
	for _, tc := range []struct{ in, want string }{
		{"bioconda", "bioconda"},
		{"https://conda.anaconda.org/bioconda/linux-64", "bioconda"},
		{"https://conda.anaconda.org/conda-forge/noarch", "conda-forge"},
		{"conda-forge/linux-64", "conda-forge"},
		{"", ""},
	} {
		if got := channelName(tc.in); got != tc.want {
			t.Errorf("channelName(%q) = %q, want %q", tc.in, got, tc.want)
		}
	}
}

// FETCH lists only what is missing from the package cache, so a warm cache would
// silently shrink the package set — and with it the identity.
func TestResolvedPrefersLinkOverFetch(t *testing.T) {
	var d DryRun
	d.Actions.Fetch = []Package{{Name: "only-uncached"}}
	d.Actions.Link = []Package{{Name: "a"}, {Name: "b"}}
	if got := d.Resolved(); len(got) != 2 || got[0].Name != "a" {
		t.Fatalf("Resolved() = %#v, want the LINK set", got)
	}
	d.Actions.Link = nil
	if got := d.Resolved(); len(got) != 1 || got[0].Name != "only-uncached" {
		t.Fatalf("Resolved() = %#v, want the FETCH fallback", got)
	}
}

// A package with no URL cannot go in an explicit file, and guessing one would
// invent an identity.
func TestExplicitFromRefusesAPackageWithoutAURL(t *testing.T) {
	_, err := ExplicitFrom([]Package{{Name: "mystery", Version: "1.0"}})
	if err == nil || !strings.Contains(err.Error(), "mystery") {
		t.Fatalf("error = %v, want one naming the package", err)
	}
}
