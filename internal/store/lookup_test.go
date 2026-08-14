package store

import (
	"errors"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/internal/artifact/compare"
	"github.com/Justype/condatainer/internal/artifactcache"
)

func TestScanVerifiesAddressAndResolveDetectsAmbiguity(t *testing.T) {
	root := t.TempDir()
	storeDir := filepath.Join(root, DirName)
	if err := os.Mkdir(storeDir, 0o755); err != nil {
		t.Fatal(err)
	}
	a := strings.Repeat("a", 64)
	b := "aaaaaaaaaaaa" + strings.Repeat("b", 52)
	paths := map[string]compare.Artifact{}
	add := func(path, name, identity string) {
		t.Helper()
		if err := os.WriteFile(path, []byte("fixture"), 0o644); err != nil {
			t.Fatal(err)
		}
		paths[path] = compare.Artifact{
			Name: name, IdentityScheme: "script-identity-v1", Identity: "sha256:" + identity,
			EquivScheme: "script-equiv-v1", Equiv: "sha256:" + strings.Repeat("e", 64),
		}
	}
	flat := filepath.Join(root, "star--2.7.sqf")
	storedA := filepath.Join(storeDir, "star--2.7@aaaaaaaaaaaa.sqf")
	storedB := filepath.Join(storeDir, "star--2.7@aaaaaaaaaaaab.sqf")
	bad := filepath.Join(storeDir, "star--2.7@cccccccccccc.sqf")
	add(flat, "star/2.7", a)
	add(storedA, "star/2.7", a)
	add(storedB, "star/2.7", b)
	add(bad, "star/2.7", a)

	read := func(path string) (compare.Artifact, error) {
		artifact, ok := paths[path]
		if !ok {
			return compare.Artifact{}, errors.New("missing fixture")
		}
		return artifact, nil
	}
	report := scan(ScanOptions{Dirs: []string{root}, Name: "star/2.7"}, read, nil)
	if len(report.Candidates) != 3 || len(report.Issues) != 1 {
		t.Fatalf("scan got %d candidates and %d issues: %#v", len(report.Candidates), len(report.Issues), report)
	}
	if report.Candidates[0].Layout != LayoutFlat || report.Candidates[0].Path != flat {
		t.Fatalf("flat candidate did not retain priority: %#v", report.Candidates[0])
	}

	candidate, err := resolveReport("star/2.7", IdentityQuery{SHA256: strings.Repeat("a", 64)}, report)
	if err != nil || candidate.Path != flat {
		t.Fatalf("exact duplicate resolution = %#v, %v", candidate, err)
	}
	if _, err := resolveReport("star/2.7", IdentityQuery{SHA256: "aaaaaaaaaaaa"}, report); !errors.Is(err, ErrAmbiguous) {
		t.Fatalf("prefix resolution error = %v", err)
	}
}

func TestScanUsesVerifiedFingerprintCache(t *testing.T) {
	root := t.TempDir()
	path := filepath.Join(root, "star--2.7.sqf")
	if err := os.WriteFile(path, []byte("fixture"), 0o644); err != nil {
		t.Fatal(err)
	}
	cache := artifactcache.New(func() string { return "" })
	reads := 0
	read := func(string) (compare.Artifact, error) {
		reads++
		return compare.Artifact{
			Name: "star/2.7", IdentityScheme: "script-identity-v1", Identity: "sha256:" + strings.Repeat("a", 64),
			EquivScheme: "script-equiv-v1", Equiv: "sha256:" + strings.Repeat("e", 64),
		}, nil
	}
	opts := ScanOptions{Dirs: []string{root}, Name: "star/2.7", Flat: true}
	if got := scan(opts, read, cache); len(got.Candidates) != 1 || reads != 1 {
		t.Fatalf("cold scan = %#v, reads = %d", got, reads)
	}
	if got := scan(opts, read, cache); len(got.Candidates) != 1 || reads != 1 {
		t.Fatalf("warm scan = %#v, reads = %d; wanted cached keys", got, reads)
	}
	if got := scan(opts, read, nil); len(got.Candidates) != 1 || reads != 2 {
		t.Fatalf("uncached validation scan = %#v, reads = %d; wanted regenerated keys", got, reads)
	}
}
