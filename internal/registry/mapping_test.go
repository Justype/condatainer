package registry

import (
	"slices"
	"strings"
	"testing"
	"time"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
)

var buildTime = time.Date(2026, 7, 21, 13, 45, 0, 0, time.UTC)

func manifestFor(typ catalog.Type, name string) meta.Manifest {
	return meta.Manifest{
		Name:  name,
		Type:  typ,
		Build: meta.Build{Created: buildTime},
	}
}

// Version-less is a property of type *and* name shape, not of either alone:
// ubuntu24/build-essential and hello/1.0 are the same shape and map differently.
func TestIsVersionLess(t *testing.T) {
	tests := []struct {
		typ  catalog.Type
		name string
		want bool
	}{
		{catalog.TypeBase, "ubuntu24/base", true},
		{catalog.TypeOS, "ubuntu24/build-essential", true},
		{catalog.TypeOS, "ubuntu24/r/4.4.4", false}, // a real version segment
		{catalog.TypeApp, "hello/1.0", false},       // same shape as the OS above
		{catalog.TypeApp, "cellranger/9.0.1", false},
		{catalog.TypeData, "grch38/star/2.7.11b/gencode49-101", false},
		// No slash at all: nothing in the name is a version, whatever the type.
		{catalog.TypeApp, "myenv", true},
		{catalog.TypeData, "scratchset", true},
		{catalog.TypeOS, "ubuntu24", true},
	}
	for _, tt := range tests {
		if got := isVersionLess(tt.typ, tt.name); got != tt.want {
			t.Errorf("isVersionLess(%s, %q) = %v, want %v", tt.typ, tt.name, got, tt.want)
		}
	}
}

// The versioned branch splits at the last segment; the version-less branch keeps
// the whole name as the repository and dates the tag.
func TestPushReference(t *testing.T) {
	tests := []struct {
		typ      catalog.Type
		name     string
		wantRepo string
		wantTags []string
	}{
		{catalog.TypeApp, "cellranger/9.0.1", "cellranger", []string{"9.0.1"}},
		{catalog.TypeData, "grch38/star/2.7.11b/gencode49-101",
			"grch38/star/2.7.11b", []string{"gencode49-101"}},
		{catalog.TypeOS, "ubuntu24/r/4.4.4", "ubuntu24/r", []string{"4.4.4"}},
		{catalog.TypeBase, "ubuntu24/base", "ubuntu24/base",
			[]string{"20260721", RollingTag}},
		{catalog.TypeOS, "ubuntu24/build-essential", "ubuntu24/build-essential",
			[]string{"20260721", RollingTag}},
		// A conda environment built from a file needs no version, so it rolls.
		{catalog.TypeApp, "myenv", "myenv", []string{"20260721", RollingTag}},
	}
	for _, tt := range tests {
		repo, tags, err := PushReference(manifestFor(tt.typ, tt.name))
		if err != nil {
			t.Errorf("PushReference(%q): %v", tt.name, err)
			continue
		}
		if repo != tt.wantRepo || strings.Join(tags, ",") != strings.Join(tt.wantTags, ",") {
			t.Errorf("PushReference(%q) = (%q, %v), want (%q, %v)",
				tt.name, repo, tags, tt.wantRepo, tt.wantTags)
		}
	}
}

// A date tag is an address. One recording the upload rather than the build is a
// lie that only surfaces as a wrong artifact months later, so refuse it.
func TestPushReferenceRefusesAVersionLessArtifactWithNoBuildTime(t *testing.T) {
	m := manifestFor(catalog.TypeBase, "ubuntu24/base")
	m.Build.Created = time.Time{}
	if _, _, err := PushReference(m); err == nil {
		t.Fatal("a base with no build time was given a date tag anyway")
	}

	// A versioned artifact takes its tag from the name, so it needs no build time.
	versioned := manifestFor(catalog.TypeApp, "cellranger/9.0.1")
	versioned.Build.Created = time.Time{}
	if _, _, err := PushReference(versioned); err != nil {
		t.Errorf("a versioned artifact should not need a build time: %v", err)
	}
}

// A date tag names a day in UTC wherever it was built, or two sites publishing
// the same artifact would disagree about which day it was.
func TestPushReferenceDateTagIsUTC(t *testing.T) {
	m := manifestFor(catalog.TypeBase, "ubuntu24/base")
	// 21:00 UTC on the 21st is already the 22nd in Auckland.
	m.Build.Created = time.Date(2026, 7, 21, 21, 0, 0, 0,
		time.FixedZone("NZST", 13*60*60))

	_, tags, err := PushReference(m)
	if err != nil {
		t.Fatal(err)
	}
	if tags[0] != "20260721" {
		t.Errorf("date tag = %q, want 20260721 — the UTC day", tags[0])
	}
}

func TestPullReference(t *testing.T) {
	tests := []struct {
		typ            catalog.Type
		name           string
		wantRepo, want string
	}{
		{catalog.TypeApp, "cellranger/9.0.1", "cellranger", "9.0.1"},
		{catalog.TypeData, "grch38/star/2.7.11b/gencode47-101", "grch38/star/2.7.11b", "gencode47-101"},
		{catalog.TypeBase, "ubuntu24/base", "ubuntu24/base", RollingTag},
		{catalog.TypeOS, "ubuntu24/build-essential", "ubuntu24/build-essential", RollingTag},
		{catalog.TypeOS, "ubuntu24/r/4.4.4", "ubuntu24/r", "4.4.4"},
	}
	for _, tt := range tests {
		repo, tag, err := PullReference(tt.typ, tt.name)
		if err != nil {
			t.Errorf("PullReference(%s, %q): %v", tt.typ, tt.name, err)
			continue
		}
		if repo != tt.wantRepo || tag != tt.want {
			t.Errorf("PullReference(%s, %q) = (%q, %q), want (%q, %q)",
				tt.typ, tt.name, repo, tag, tt.wantRepo, tt.want)
		}
	}
}

// Push and pull must agree on the repository for the same artifact, or a pull
// looks in a place nothing was ever published to.
func TestPushAndPullAgreeOnTheRepository(t *testing.T) {
	for _, tt := range []struct {
		typ  catalog.Type
		name string
	}{
		{catalog.TypeApp, "cellranger/9.0.1"},
		{catalog.TypeData, "grch38/star/2.7.11b/gencode49-101"},
		{catalog.TypeOS, "ubuntu24/build-essential"},
		{catalog.TypeOS, "ubuntu24/r/4.4.4"},
		{catalog.TypeBase, "ubuntu24/base"},
	} {
		pushRepo, pushTags, err := PushReference(manifestFor(tt.typ, tt.name))
		if err != nil {
			t.Fatalf("PushReference(%q): %v", tt.name, err)
		}
		pullRepo, pullTag, err := PullReference(tt.typ, tt.name)
		if err != nil {
			t.Fatalf("PullReference(%q): %v", tt.name, err)
		}
		if pushRepo != pullRepo {
			t.Errorf("%s: push repo %q, pull repo %q", tt.name, pushRepo, pullRepo)
		}
		// A pull resolves the tag it can know: the pushed one when versioned,
		// the rolling one otherwise.
		if !slices.Contains(pushTags, pullTag) {
			t.Errorf("%s: pull tag %q is not among pushed tags %v", tt.name, pullTag, pushTags)
		}
	}
}

// A slashless name is not an error — it is version-less, so it rolls.
func TestPullReferenceSlashlessNameRolls(t *testing.T) {
	repo, tag, err := PullReference(catalog.TypeApp, "myenv")
	if err != nil {
		t.Fatalf("PullReference: %v", err)
	}
	if repo != "myenv" || tag != RollingTag {
		t.Errorf("= (%q, %q), want (myenv, %s)", repo, tag, RollingTag)
	}
}

func TestPullReferenceRejectsUnmappableNames(t *testing.T) {
	if _, _, err := PullReference(catalog.TypeData, "GRCh38/STAR/2.7.11b"); err == nil {
		t.Error("an uppercase repository was accepted")
	}
	if _, _, err := PullReference(catalog.TypeApp, "cell ranger/1.0"); err == nil {
		t.Error("a name with a space was accepted")
	}
	if _, _, err := PullReference(catalog.TypeApp, ""); err == nil {
		t.Error("an empty name was accepted")
	}
}

func TestValidateIdentity(t *testing.T) {
	valid := []struct {
		name        string
		versionLess bool
	}{
		{"cellranger/9.0.1", false},
		{"grch38/star/2.7.11b/gencode47-101", false},
		{"ubuntu24/base", true},
		{"ubuntu24/build-essential", true},
	}
	for _, tc := range valid {
		if err := ValidateIdentity(tc.name, tc.versionLess); err != nil {
			t.Errorf("ValidateIdentity(%q, %v): %v", tc.name, tc.versionLess, err)
		}
	}

	invalid := []struct {
		why         string
		name        string
		versionLess bool
	}{
		{"uppercase repository", "CellRanger/9.0.1", false},
		{"traversal as a tag", "cellranger/..", false},
		{"a digest selector inside the name", "cellranger/1.0@sha256:x", false},
		{"a space", "cell ranger/1.0", false},
		{"uppercase in a version-less name", "ubuntu24/Base_Image", true},
		{"empty", "", false},
		{"no version segment", "cellranger", false},
	}
	for _, tc := range invalid {
		if err := ValidateIdentity(tc.name, tc.versionLess); err == nil {
			t.Errorf("ValidateIdentity(%q, %v) accepted %s", tc.name, tc.versionLess, tc.why)
		}
	}
}

func TestSplitPullSpec(t *testing.T) {
	sha := digestPrefix + strings.Repeat("a", 64)
	tests := []struct{ in, name, selector string }{
		{"ubuntu24/base", "ubuntu24/base", ""},
		{"ubuntu24/base:20260721", "ubuntu24/base", "20260721"},
		{"ubuntu24/base:latest", "ubuntu24/base", "latest"},
		{"ubuntu24/base@" + sha, "ubuntu24/base", sha},
		{"cellranger/9.0.1", "cellranger/9.0.1", ""},
		// A separator the catalog accepts normalizes before the split.
		{"cellranger--9.0.1", "cellranger/9.0.1", ""},
	}
	for _, tt := range tests {
		name, selector, err := SplitPullSpec(tt.in)
		if err != nil {
			t.Errorf("SplitPullSpec(%q): %v", tt.in, err)
			continue
		}
		if name != tt.name || selector != tt.selector {
			t.Errorf("SplitPullSpec(%q) = (%q, %q), want (%q, %q)",
				tt.in, name, selector, tt.name, tt.selector)
		}
	}
}

func TestSplitPullSpecRejectsAMalformedDigest(t *testing.T) {
	if _, _, err := SplitPullSpec("ubuntu24/base@sha256:nope"); err == nil {
		t.Fatal("a malformed digest was accepted")
	}
}
