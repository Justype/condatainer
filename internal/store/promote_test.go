package store

import (
	"errors"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/internal/artifact/compare"
	"github.com/Justype/condatainer/internal/image"
)

// fixture builds a reader over a map of content → artifact, so promotion can be
// tested without real SquashFS. Keyed by content rather than path because
// promotion renames: an artifact is the same artifact wherever it lands.
type fixture struct {
	t         *testing.T
	artifacts map[string]compare.Artifact
}

func newFixture(t *testing.T) *fixture {
	return &fixture{t: t, artifacts: map[string]compare.Artifact{}}
}

func (f *fixture) write(path, name, identity string) string {
	f.t.Helper()
	if err := os.MkdirAll(filepath.Dir(path), 0o775); err != nil {
		f.t.Fatal(err)
	}
	if err := os.WriteFile(path, []byte(identity), 0o664); err != nil {
		f.t.Fatal(err)
	}
	f.artifacts[identity] = transactionArtifact(name, identity)
	return path
}

// flat writes a bare-name artifact in root.
func (f *fixture) flat(root, name, identity string) string {
	return f.write(filepath.Join(root, image.EncodeArtifactName(name)+".sqf"), name, identity)
}

// stored writes an identity-addressed artifact in root's store/.
func (f *fixture) stored(root, name, identity string) string {
	filename, err := Filename(name, keyRef("identity-v1", "sha256:"+identity), DefaultPrefixChars)
	if err != nil {
		f.t.Fatal(err)
	}
	return f.write(filepath.Join(root, DirName, filename), name, identity)
}

func (f *fixture) read(path string) (compare.Artifact, error) {
	data, err := os.ReadFile(path)
	if err != nil {
		return compare.Artifact{}, err
	}
	if artifact, ok := f.artifacts[string(data)]; ok {
		return artifact, nil
	}
	return compare.Artifact{}, errors.New("not an artifact")
}

func query(identity string) IdentityQuery {
	return IdentityQuery{Scheme: "identity-v1", SHA256: identity}
}

var (
	idA = strings.Repeat("a", 64)
	idB = strings.Repeat("b", 64)
)

// The ordinary swap: the store entry takes the bare name and the incumbent goes
// to the same root's store/, where it is still resolvable by identity.
func TestPromoteSwapsInsideOneRoot(t *testing.T) {
	root := t.TempDir()
	f := newFixture(t)
	flat := f.flat(root, "star/2.7", idA)
	entry := f.stored(root, "star/2.7", idB)

	result, err := promote("star/2.7", query(idB), []string{root}, f.read)
	if err != nil {
		t.Fatal(err)
	}
	if result.Promoted.Path != flat {
		t.Errorf("promoted to %q, want %q", result.Promoted.Path, flat)
	}
	if result.DemotedTo == "" || filepath.Dir(result.DemotedTo) != filepath.Join(root, DirName) {
		t.Errorf("demoted to %q, want a path in the root's store/", result.DemotedTo)
	}
	if _, err := os.Stat(entry); !os.IsNotExist(err) {
		t.Errorf("the promoted entry is still in the store: %v", err)
	}
	// Both identities still resolve, which is what keeps other people's locks working.
	for _, id := range []string{idA, idB} {
		if _, _, err := resolveWith(t, f, root, id); err != nil {
			t.Errorf("identity %s… no longer resolves: %v", id[:6], err)
		}
	}
}

func resolveWith(t *testing.T, f *fixture, root, identity string) (Candidate, Report, error) {
	t.Helper()
	report := scan(ScanOptions{Dirs: []string{root}, Name: "star/2.7"}, f.read, nil)
	candidate, err := resolveReport("star/2.7", query(identity), report)
	return candidate, report, err
}

// A root with no flat copy needs no demotion: one rename, and the farther copy
// is left in place and reported as shadowed.
func TestPromoteInNearerRootShadowsRatherThanMoves(t *testing.T) {
	near, far := t.TempDir(), t.TempDir()
	f := newFixture(t)
	farFlat := f.flat(far, "star/2.7", idA)
	f.stored(near, "star/2.7", idB)

	result, err := promote("star/2.7", query(idB), []string{near, far}, f.read)
	if err != nil {
		t.Fatal(err)
	}
	if result.DemotedTo != "" {
		t.Errorf("demoted %q, but the nearer root held no flat copy", result.DemotedTo)
	}
	if result.Shadows != farFlat {
		t.Errorf("shadows = %q, want %q", result.Shadows, farFlat)
	}
	if _, err := os.Stat(farFlat); err != nil {
		t.Errorf("the shadowed copy was disturbed: %v", err)
	}
	if want := filepath.Join(near, "star--2.7.sqf"); result.Promoted.Path != want {
		t.Errorf("promoted to %q, want %q", result.Promoted.Path, want)
	}
}

// Promoting where a nearer root already answers to the name would change
// nothing a read resolves, so it is refused instead of silently doing nothing.
func TestPromoteRefusesWhenANearerRootShadows(t *testing.T) {
	near, far := t.TempDir(), t.TempDir()
	f := newFixture(t)
	nearFlat := f.flat(near, "star/2.7", idA)
	entry := f.stored(far, "star/2.7", idB)

	_, err := promote("star/2.7", query(idB), []string{near, far}, f.read)
	if !errors.Is(err, ErrShadowed) {
		t.Fatalf("error = %v, want ErrShadowed", err)
	}
	if !strings.Contains(err.Error(), nearFlat) {
		t.Errorf("error does not name the shadowing file: %v", err)
	}
	if _, err := os.Stat(entry); err != nil {
		t.Errorf("a refused promotion moved the entry: %v", err)
	}
}

// The requested identity already answers to the name: report it, do not rename.
func TestPromoteIsANoOpWhenAlreadyFlat(t *testing.T) {
	root := t.TempDir()
	f := newFixture(t)
	flat := f.flat(root, "star/2.7", idA)

	result, err := promote("star/2.7", query(idA), []string{root}, f.read)
	if err != nil {
		t.Fatal(err)
	}
	if !result.AlreadyFlat || result.Promoted.Path != flat {
		t.Errorf("result = %#v, want a no-op on %q", result, flat)
	}
	if _, err := os.Stat(filepath.Join(root, DirName)); !os.IsNotExist(err) {
		t.Errorf("a no-op created a store directory: %v", err)
	}
}

// A cleared write bit is how an artifact is pinned, and it protects the
// incumbent from being demoted out from under whoever pinned it.
func TestPromoteRefusesAProtectedIncumbent(t *testing.T) {
	root := t.TempDir()
	f := newFixture(t)
	flat := f.flat(root, "star/2.7", idA)
	entry := f.stored(root, "star/2.7", idB)
	if err := os.Chmod(flat, 0o444); err != nil {
		t.Fatal(err)
	}
	defer os.Chmod(flat, 0o664) //nolint:errcheck

	_, err := promote("star/2.7", query(idB), []string{root}, f.read)
	if !errors.Is(err, image.ErrProtected) {
		t.Fatalf("error = %v, want ErrProtected", err)
	}
	if _, err := os.Stat(entry); err != nil {
		t.Errorf("a refused promotion moved the entry: %v", err)
	}
	if _, err := os.Stat(flat); err != nil {
		t.Errorf("a refused promotion moved the incumbent: %v", err)
	}
}

// The interrupted state — both artifacts in store/, the bare name free — is
// what a crash between the two renames leaves, and re-running repairs it.
func TestPromoteRepairsAnInterruptedSwap(t *testing.T) {
	root := t.TempDir()
	f := newFixture(t)
	f.stored(root, "star/2.7", idA)
	f.stored(root, "star/2.7", idB)

	result, err := promote("star/2.7", query(idB), []string{root}, f.read)
	if err != nil {
		t.Fatal(err)
	}
	if want := filepath.Join(root, "star--2.7.sqf"); result.Promoted.Path != want {
		t.Errorf("promoted to %q, want %q", result.Promoted.Path, want)
	}
	if result.DemotedTo != "" {
		t.Errorf("demoted %q, but nothing held the bare name", result.DemotedTo)
	}
}

// A demoted incumbent whose 12-character prefix is taken gets a longer one
// rather than colliding.
func TestDemoteLengthensACollidingPrefix(t *testing.T) {
	root := t.TempDir()
	f := newFixture(t)
	flat := f.flat(root, "star/2.7", idA)
	f.stored(root, "star/2.7", idB)

	// Squat the incumbent's natural store filename with something else.
	occupied, err := Filename("star/2.7", keyRef("identity-v1", "sha256:"+idA), DefaultPrefixChars)
	if err != nil {
		t.Fatal(err)
	}
	f.write(filepath.Join(root, DirName, occupied), "star/2.7", strings.Repeat("c", 64))

	result, err := promote("star/2.7", query(idB), []string{root}, f.read)
	if err != nil {
		t.Fatal(err)
	}
	if filepath.Base(result.DemotedTo) == occupied {
		t.Fatalf("demotion overwrote the occupied filename %q", occupied)
	}
	if len(result.DemotedTo) == 0 {
		t.Fatal("nothing was demoted")
	}
	if _, err := os.Stat(flat); err == nil {
		if data, _ := os.ReadFile(flat); string(data) == idA {
			t.Error("the incumbent still holds the bare name")
		}
	}
}
