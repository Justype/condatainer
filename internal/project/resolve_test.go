package project

import (
	"context"
	"encoding/json"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/capsule"
	"github.com/Justype/condatainer/internal/artifact/key"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/project/lock"
	"github.com/Justype/condatainer/internal/store"
)

func projectRoot(t *testing.T) string {
	t.Helper()
	root := t.TempDir()
	if err := os.MkdirAll(filepath.Join(root, lock.DirName), 0o775); err != nil {
		t.Fatal(err)
	}
	return root
}

func vendor(t *testing.T, root, name, recipe string) (string, meta.Manifest) {
	t.Helper()
	manifest := meta.Manifest{
		SchemaVersion: meta.SchemaVersion,
		Name:          name,
		Type:          catalog.TypeApp,
		BuildType:     "script",
		Platform:      meta.NativePlatform(),
		Source:        meta.Source{Files: []string{meta.RecipeFileName}},
	}
	derived, err := key.Generate(manifest, key.Sources{meta.RecipeFileName: []byte(recipe)})
	if err != nil {
		t.Fatal(err)
	}
	manifest.Keys = derived.Keys()
	data, err := json.MarshalIndent(manifest, "", "  ")
	if err != nil {
		t.Fatal(err)
	}
	relative, err := lock.StageEntry(root, capsule.EntryName(name, manifest.Keys.Identity.Digest()),
		map[string][]byte{meta.FileName: data, meta.RecipeFileName: []byte(recipe)})
	if err != nil {
		t.Fatal(err)
	}
	return relative, manifest
}

func present(name, path string) LookupFunc {
	return func(got string, keys meta.Keys, _ Match, _ []string) (store.Candidate, bool) {
		if got != name {
			return store.Candidate{}, false
		}
		return store.Candidate{Name: got, Path: path, Identity: keys.Identity, Equiv: keys.Equiv}, true
	}
}

func absent(string, meta.Keys, Match, []string) (store.Candidate, bool) {
	return store.Candidate{}, false
}

func absentAt(string, string, meta.Keys, Match) (store.Candidate, bool) {
	return store.Candidate{}, false
}

func TestResolveMountsALockedName(t *testing.T) {
	root := projectRoot(t)
	relative, _ := vendor(t, root, "star/2.7.11b", "echo star\n")
	l := lock.New()
	l.Pins["star/2.7.11b"] = lock.PinEntry{Artifact: relative}
	requests := []lock.Request{{Key: "star/2.7.11b", Kind: lock.KindName}}

	got, err := Resolve(context.Background(), root, l, requests, ResolveOptions{lookup: present("star/2.7.11b", "/images/star.sqf")})
	if err != nil {
		t.Fatal(err)
	}
	if !got.Complete() {
		t.Fatalf("unresolved: %#v", got.Unresolved)
	}
	if len(got.Mounts) != 1 || got.Mounts[0].Path != "/images/star.sqf" {
		t.Fatalf("mounts = %#v", got.Mounts)
	}
	if got.Mounts[0].Identity == "" {
		t.Error("the locked identity was not reported")
	}
}

// Inside a project an unsatisfiable declaration is always an error. Falling back
// to whatever answers to the name is the failure a lock exists to prevent.
func TestResolveRefusesRatherThanFallingBackToTheName(t *testing.T) {
	root := projectRoot(t)
	relative, _ := vendor(t, root, "star/2.7.11b", "echo star\n")
	l := lock.New()
	l.Pins["star/2.7.11b"] = lock.PinEntry{Artifact: relative}
	requests := []lock.Request{{Key: "star/2.7.11b", Kind: lock.KindName}}

	got, err := Resolve(context.Background(), root, l, requests, ResolveOptions{lookup: absent})
	if err != nil {
		t.Fatal(err)
	}
	if got.Complete() {
		t.Fatal("resolution succeeded with nothing available")
	}
	if len(got.Mounts) != 0 {
		t.Errorf("mounts = %#v, want none", got.Mounts)
	}
	if !strings.Contains(got.Unresolved[0].Reason, "star/2.7.11b") {
		t.Errorf("reason does not name the artifact: %q", got.Unresolved[0].Reason)
	}
}

// A declaration nothing selected cannot be mounted, and the remedy is to pin it.
func TestResolveReportsAnUnpinnedDeclaration(t *testing.T) {
	root := projectRoot(t)
	l := lock.New()
	requests := []lock.Request{{Key: "star/2.7.11b", Kind: lock.KindName}}

	got, err := Resolve(context.Background(), root, l, requests, ResolveOptions{lookup: absent})
	if err != nil {
		t.Fatal(err)
	}
	if got.Complete() {
		t.Fatal("an unpinned declaration resolved")
	}
	if !strings.Contains(got.Unresolved[0].Reason, "project pin") {
		t.Errorf("reason does not name the remedy: %q", got.Unresolved[0].Reason)
	}
}

// A bare name no pin answers resolves live when LiveResolve is set, the same
// as it would outside a project — never claimed to be pinned, so nothing was
// promised and broken.
func TestResolveLiveResolvesALooseNameWhenEnabled(t *testing.T) {
	root := projectRoot(t)
	imagesDir := t.TempDir()
	if err := os.WriteFile(filepath.Join(imagesDir, "samtools--1.22.1.sqf"), nil, 0o644); err != nil {
		t.Fatal(err)
	}
	prevSources := config.Global.Sources
	config.Global.Sources = nil
	config.ResetCatalog()
	t.Cleanup(func() {
		config.Global.Sources = prevSources
		config.ResetCatalog()
	})

	l := lock.New()
	requests := []lock.Request{{Key: "samtools", Kind: lock.KindName, Dep: catalog.Dep{Name: "samtools"}}}

	got, err := Resolve(context.Background(), root, l, requests,
		ResolveOptions{lookup: absent, LiveResolve: true, SearchDirs: []string{imagesDir}})
	if err != nil {
		t.Fatal(err)
	}
	if !got.Complete() {
		t.Fatalf("unresolved: %#v", got.Unresolved)
	}
	want := filepath.Join(imagesDir, "samtools--1.22.1.sqf")
	if len(got.Mounts) != 1 || got.Mounts[0].Path != want || got.Mounts[0].Name != "samtools/1.22.1" {
		t.Fatalf("mounts = %#v, want a live mount at %s", got.Mounts, want)
	}
	if !got.Mounts[0].Live {
		t.Error("Live was not set on a live-resolved mount")
	}
}

// The same unpinned bare name still refuses without LiveResolve — the fixed
// requirements a project resolves for (a helper's required overlays, the
// project's own root) never opt into it.
func TestResolveRefusesALooseNameWhenLiveResolveDisabled(t *testing.T) {
	root := projectRoot(t)
	imagesDir := t.TempDir()
	if err := os.WriteFile(filepath.Join(imagesDir, "samtools--1.22.1.sqf"), nil, 0o644); err != nil {
		t.Fatal(err)
	}

	l := lock.New()
	requests := []lock.Request{{Key: "samtools", Kind: lock.KindName, Dep: catalog.Dep{Name: "samtools"}}}

	got, err := Resolve(context.Background(), root, l, requests,
		ResolveOptions{lookup: absent, SearchDirs: []string{imagesDir}})
	if err != nil {
		t.Fatal(err)
	}
	if got.Complete() {
		t.Fatal("resolved live without LiveResolve set")
	}
}

// A path pin is verified at its declared project path, not in the store.
func TestResolveAnchorsAPathSelectionOnTheProjectRoot(t *testing.T) {
	root := projectRoot(t)
	relative, _ := vendor(t, root, "tool/1.0", "echo tool\n")
	l := lock.New()
	l.Pins[lock.PathPrefix+"overlays/tool.sqf"] = lock.PinEntry{Artifact: relative}
	requests := []lock.Request{{
		Key: lock.PathPrefix + "overlays/tool.sqf", Kind: lock.KindPath, Path: "overlays/tool.sqf",
	}}

	want := filepath.Join(root, "overlays", "tool.sqf")
	var asked string
	lookupAt := func(path, name string, keys meta.Keys, _ Match) (store.Candidate, bool) {
		asked = path
		return store.Candidate{Name: name, Path: path, Identity: keys.Identity, Equiv: keys.Equiv}, true
	}

	got, err := Resolve(context.Background(), root, l, requests, ResolveOptions{lookup: absent, lookupAt: lookupAt})
	if err != nil {
		t.Fatal(err)
	}
	if asked != want {
		t.Fatalf("looked at %q, want the project-root anchored %q", asked, want)
	}
	if !got.Complete() || got.Mounts[0].Path != want {
		t.Fatalf("mounts = %#v", got.Mounts)
	}
}

// A path pinned under a declaring script's own directory — because the
// root-relative file did not exist when it was locked — is still verified at
// the resolved location, not at the literal text as scanned.
func TestResolveMountsAPathPinnedUnderAFallbackKey(t *testing.T) {
	root := projectRoot(t)
	relative, _ := vendor(t, root, "tool/1.0", "echo tool\n")
	l := lock.New()
	l.Pins[lock.PathPrefix+"steps1/xxx.sqf"] = lock.PinEntry{Artifact: relative}
	requests := []lock.Request{{
		Key: lock.PathPrefix + "xxx.sqf", Kind: lock.KindPath, Path: "xxx.sqf",
		Scripts: []string{"steps1/run.sh"},
	}}

	want := filepath.Join(root, "steps1", "xxx.sqf")
	var asked string
	lookupAt := func(path, name string, keys meta.Keys, _ Match) (store.Candidate, bool) {
		asked = path
		return store.Candidate{Name: name, Path: path, Identity: keys.Identity, Equiv: keys.Equiv}, true
	}

	got, err := Resolve(context.Background(), root, l, requests, ResolveOptions{lookup: absent, lookupAt: lookupAt})
	if err != nil {
		t.Fatal(err)
	}
	if asked != want {
		t.Fatalf("looked at %q, want the fallback location %q", asked, want)
	}
	if !got.Complete() || got.Mounts[0].Path != want {
		t.Fatalf("mounts = %#v", got.Mounts)
	}
}

// An unpinnable declaration is mounted as the literal path it names, and a
// project-relative one is still anchored on the root.
func TestResolveMountsUnpinnableDeclarationsLiterally(t *testing.T) {
	root := projectRoot(t)
	l := lock.New()
	requests := []lock.Request{
		{Key: lock.PathPrefix + "env.img", Kind: lock.KindWritable, Path: "env.img"},
		{Key: lock.PathPrefix + "/shared/genome.sqf", Kind: lock.KindExternal, Path: "/shared/genome.sqf"},
	}

	got, err := Resolve(context.Background(), root, l, requests, ResolveOptions{lookup: absent, lookupAt: absentAt})
	if err != nil {
		t.Fatal(err)
	}
	if !got.Complete() {
		t.Fatalf("unresolved: %#v", got.Unresolved)
	}
	if len(got.Mounts) != 2 {
		t.Fatalf("mounts = %#v", got.Mounts)
	}
	if got.Mounts[0].Path != filepath.Join(root, "env.img") || !got.Mounts[0].Unpinned {
		t.Errorf("writable mount = %#v, want it anchored on the root", got.Mounts[0])
	}
	if got.Mounts[1].Path != "/shared/genome.sqf" {
		t.Errorf("external mount = %#v, want the path as written", got.Mounts[1])
	}
}

// Every returned path is absolute: the runtime resolves a relative overlay
// against the process working directory, which a scheduler chooses.
func TestResolveReturnsAbsolutePaths(t *testing.T) {
	root := projectRoot(t)
	l := lock.New()
	requests := []lock.Request{
		{Key: lock.PathPrefix + "env.img", Kind: lock.KindWritable, Path: "env.img"},
	}
	got, err := Resolve(context.Background(), root, l, requests, ResolveOptions{lookup: absent, lookupAt: absentAt})
	if err != nil {
		t.Fatal(err)
	}
	for _, mount := range got.Mounts {
		if !filepath.IsAbs(mount.Path) {
			t.Errorf("mount %q is relative", mount.Path)
		}
	}
}

// An unsound lock is reported, never resolved around.
func TestResolveReportsAnInvalidLock(t *testing.T) {
	root := projectRoot(t)
	l := lock.New()
	l.Pins["star/2.7.11b"] = lock.PinEntry{Artifact: "provenance/star--2.7.11b@000000000000"}
	requests := []lock.Request{{Key: "star/2.7.11b", Kind: lock.KindName}}

	got, err := Resolve(context.Background(), root, l, requests, ResolveOptions{lookup: absent})
	if err != nil {
		t.Fatal(err)
	}
	if got.Complete() {
		t.Fatal("an invalid lock resolved")
	}
}

// An equivalent substitution is reported, never silent.
func TestResolveReportsASubstitution(t *testing.T) {
	root := projectRoot(t)
	relative, manifest := vendor(t, root, "star/2.7.11b", "echo star\n")
	l := lock.New()
	l.Pins["star/2.7.11b"] = lock.PinEntry{Artifact: relative}
	requests := []lock.Request{{Key: "star/2.7.11b", Kind: lock.KindName}}

	other := meta.KeyRef{Scheme: manifest.Keys.Identity.Scheme, SHA256: strings.Repeat("e", 64)}
	lookup := func(name string, keys meta.Keys, _ Match, _ []string) (store.Candidate, bool) {
		return store.Candidate{Name: name, Path: "/images/star.sqf", Identity: other, Equiv: keys.Equiv}, true
	}

	got, err := Resolve(context.Background(), root, l, requests, ResolveOptions{lookup: lookup})
	if err != nil {
		t.Fatal(err)
	}
	if !got.Complete() {
		t.Fatalf("unresolved: %#v", got.Unresolved)
	}
	if got.Mounts[0].Found != other.Digest() {
		t.Errorf("substitution was not reported: %#v", got.Mounts[0])
	}
}

// A pinnable kind always needs a pin. Nothing in a script opts out of one, so a
// project path with no pin does not resolve.
func TestResolveRequiresAPinForAPinnablePath(t *testing.T) {
	root := projectRoot(t)
	l := lock.New()
	requests := []lock.Request{{
		Key: lock.PathPrefix + "overlays/tool.sqf", Kind: lock.KindPath, Path: "overlays/tool.sqf",
	}}

	got, err := Resolve(context.Background(), root, l, requests, ResolveOptions{lookup: absent, lookupAt: absentAt})
	if err != nil {
		t.Fatal(err)
	}
	if got.Complete() {
		t.Fatal("a project path resolved with no pin")
	}
}

func TestHigherComparesVersions(t *testing.T) {
	for _, tc := range []struct {
		a, b string
		want bool
	}{
		{"1.25", "1.21", true},
		{"1.21", "1.25", false},
		{"1.21", "1.21", false},
		{"2.7.11b", "2.7.11a", true},
		{"2.7.9a", "2.7.11b", false},
	} {
		if got := higher(tc.a, tc.b); got != tc.want {
			t.Errorf("higher(%q, %q) = %v, want %v", tc.a, tc.b, got, tc.want)
		}
	}
}

// The tool name is the dependency name without its version, which is what a
// history edge is allowed to match at any version.
func TestToolNameDropsTheVersion(t *testing.T) {
	for input, want := range map[string]string{
		"samtools/1.21":         "samtools",
		"grch38/genome/gencode": "grch38/genome",
		"ubuntu24/pytorch/2.9":  "ubuntu24/pytorch",
		// A dep recorded without a version is already just the tool.
		"samtools": "samtools",
	} {
		if got := toolName(input); got != want {
			t.Errorf("toolName(%q) = %q, want %q", input, got, want)
		}
	}
}

// Under an identity project nothing but the recorded identity is accepted, whatever
// the edge's role says: ScriptIdentityV1 hashes dependency identities, so a
// substitution would produce a dependent the mode rejects anyway.
func TestLookupInputRefusesASubstituteUnderMatchIdentity(t *testing.T) {
	dep := meta.Dependency{
		Name: "samtools/1.21", Role: meta.RoleHistory,
		Identity: meta.KeyRef{Scheme: "script-identity-v1", SHA256: strings.Repeat("a", 64)},
		Equiv:    meta.KeyRef{Scheme: "script-equiv-v1", SHA256: strings.Repeat("b", 64)},
	}
	// An empty images root, so only the exact-identity probe can answer.
	if _, ok := LookupInput(dep, MatchIdentity, []string{t.TempDir()}); ok {
		t.Fatal("a substitute was accepted under an identity project")
	}
}

// A caller that sets no mode is held to the lock's, so run and exec refuse what
// restore would have refused.
func TestResolveUsesTheLocksMatchMode(t *testing.T) {
	root := projectRoot(t)
	relative, _ := vendor(t, root, "star/2.7.11b", "echo star\n")
	l := lock.New()
	l.Pins["star/2.7.11b"] = lock.PinEntry{Artifact: relative}
	l.Match = lock.MatchIdentity
	requests := []lock.Request{{Key: "star/2.7.11b", Kind: lock.KindName}}

	var seen Match
	spy := func(_ string, _ meta.Keys, match Match, _ []string) (store.Candidate, bool) {
		seen = match
		return store.Candidate{}, false
	}
	if _, err := Resolve(context.Background(), root, l, requests, ResolveOptions{lookup: spy}); err != nil {
		t.Fatal(err)
	}
	if seen != MatchIdentity {
		t.Errorf("lookup ran under %q, want the lock's identity", seen)
	}
	if _, err := Resolve(context.Background(), root, l, requests, ResolveOptions{Match: MatchEquivalence, lookup: spy}); err != nil {
		t.Fatal(err)
	}
	if seen != MatchEquivalence {
		t.Errorf("an explicit mode did not win: %q", seen)
	}
}
