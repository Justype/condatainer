package registry

import (
	"bytes"
	"context"
	"encoding/json"
	"errors"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"sync"
	"testing"

	"github.com/opencontainers/go-digest"
	specs "github.com/opencontainers/image-spec/specs-go"
	ocispec "github.com/opencontainers/image-spec/specs-go/v1"
	"oras.land/oras-go/v2"
	"oras.land/oras-go/v2/registry"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/key"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/image"
)

func requireSquashfsTools(t *testing.T) {
	t.Helper()
	for _, bin := range []string{"mksquashfs", "unsquashfs"} {
		if _, err := exec.LookPath(bin); err != nil {
			t.Skipf("%s not available", bin)
		}
	}
}

// imageSpec describes a test artifact to build.
type imageSpec struct {
	name   string
	typ    catalog.Type
	arch   string // "" is the native architecture
	recipe string
	ext    string // "" is .sqf
	// noKeys omits the recorded keys, as an imported artifact has none.
	noKeys bool
	// tamper edits the manifest after its keys are derived, so what it records
	// and what its files reproduce disagree.
	tamper func(*meta.Manifest)
}

// packImage builds a real image carrying real metadata and real keys, because
// both gates regenerate those keys from the payload's own files — a hand-written
// fixture would prove nothing about either.
func packImage(t *testing.T, s imageSpec) (path string, m meta.Manifest) {
	t.Helper()
	if s.typ == "" {
		s.typ = catalog.TypeApp
	}
	if s.arch == "" {
		s.arch = meta.NativeArch()
	}
	if s.ext == "" {
		s.ext = ".sqf"
	}
	platform := meta.Platform{OS: "linux", Arch: s.arch}

	root := t.TempDir()
	dir := filepath.Join(root, meta.DirName)
	if err := meta.StageRuntime(dir, meta.Runtime{
		SchemaVersion: meta.SchemaVersion,
		Name:          s.name,
		Type:          s.typ,
		Platform:      platform,
		Prefix:        meta.Prefix(s.name, s.typ),
	}); err != nil {
		t.Fatal(err)
	}

	// The build type is not free to choose: an os or base is produced from an
	// Apptainer definition, an app or data from a script recipe, and the key
	// schemes refuse the other pairing outright.
	buildType := "script"
	if s.typ == catalog.TypeOS || s.typ == catalog.TypeBase {
		buildType = "def"
	}
	m = meta.Manifest{
		SchemaVersion: meta.SchemaVersion,
		Name:          s.name,
		Type:          s.typ,
		BuildType:     buildType,
		Platform:      platform,
		// Stamped as a real build stamps it. A version-less artifact has no date
		// tag without this, and no key scheme hashes it.
		Build: meta.Build{Created: buildTime},
	}
	if err := meta.StageBytes(dir, meta.RecipeFileName, []byte(s.recipe)); err != nil {
		t.Fatal(err)
	}
	m.Source.Files = []string{meta.RecipeFileName}
	if !s.noKeys {
		derived, err := key.Generate(m, key.Sources{meta.RecipeFileName: []byte(s.recipe)})
		if err != nil {
			t.Fatal(err)
		}
		m.Keys = derived.Keys()
	}
	if s.tamper != nil {
		s.tamper(&m)
	}
	if err := meta.StageManifest(dir, m); err != nil {
		t.Fatal(err)
	}

	payload := filepath.Join(root, "cnt", s.name)
	if err := os.MkdirAll(payload, 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(payload, "data"), []byte(s.recipe), 0o644); err != nil {
		t.Fatal(err)
	}

	path = filepath.Join(t.TempDir(), "image"+s.ext)
	cmd := exec.Command("mksquashfs", root, path, "-no-progress", "-noappend", "-quiet", "-no-xattrs")
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("mksquashfs: %v\n%s", err, out)
	}
	return path, m
}

// publishImage pushes an artifact into the fake registry the way push will,
// through the real chunked path, and returns what a puller resolves.
func publishImage(t *testing.T, f *fakeRegistry, repo, artifactPath string, ann map[string]string, artifactType, layerType string, tag string) ocispec.Descriptor {
	t.Helper()
	// A real derived size would leave these fixtures in one layer; 64 KiB keeps
	// the multi-layer reassembly path under test without a multi-gigabyte image.
	layers, err := pushArtifactLayers(context.Background(), fakeBlobs{f}, artifactPath, layerType, 64<<10)
	if err != nil {
		t.Fatal(err)
	}
	// oras fetches the config blob like any other.
	f.blobs[ocispec.DescriptorEmptyJSON.Digest] = ocispec.DescriptorEmptyJSON.Data

	child := f.publish(t, repo, ocispec.MediaTypeImageManifest, ocispec.Manifest{
		Versioned:    specs.Versioned{SchemaVersion: 2},
		MediaType:    ocispec.MediaTypeImageManifest,
		ArtifactType: artifactType,
		Config:       ocispec.DescriptorEmptyJSON,
		Layers:       layers,
		Annotations:  ann,
	})
	child.Platform = ptr(nativePlatform(t))
	f.publish(t, repo, ocispec.MediaTypeImageIndex, ocispec.Index{
		Versioned:   specs.Versioned{SchemaVersion: 2},
		MediaType:   ocispec.MediaTypeImageIndex,
		Manifests:   []ocispec.Descriptor{child},
		Annotations: indexAnnotations(ann),
	}, tag)
	return child
}

func ptr[T any](v T) *T { return &v }

// pullFixture is a published artifact and a place to install it.
type pullFixture struct {
	registry *fakeRegistry
	desc     ocispec.Descriptor
	ann      map[string]string
	source   string
	destDir  string
	destPath string
}

func newPullFixture(t *testing.T) *pullFixture {
	t.Helper()
	requireSquashfsTools(t)

	source, m := packImage(t, imageSpec{name: "hello/1.0", recipe: "#!/bin/bash\necho hello\n"})
	f := newFakeRegistry(t)
	ann := Annotations(m, "zstd")
	desc := publishImage(t, f, "hello", source, ann, ArtifactTypeOverlay, MediaTypeOverlayBlob, "1.0")

	destDir := t.TempDir()
	return &pullFixture{
		registry: f, desc: desc, ann: ann, source: source,
		destDir: destDir, destPath: filepath.Join(destDir, "hello--1.0.sqf"),
	}
}

func (p *pullFixture) pull(ctx context.Context) error {
	return Pull(ctx, p.registry.base(), "hello", p.desc, p.ann, p.destPath)
}

func sameBytes(t *testing.T, a, b string) bool {
	t.Helper()
	x, err := os.ReadFile(a)
	if err != nil {
		t.Fatal(err)
	}
	y, err := os.ReadFile(b)
	if err != nil {
		t.Fatal(err)
	}
	return string(x) == string(y)
}

func TestPullInstallsTheArtifact(t *testing.T) {
	p := newPullFixture(t)
	if err := p.pull(context.Background()); err != nil {
		t.Fatalf("Pull: %v", err)
	}
	if !sameBytes(t, p.source, p.destPath) {
		t.Error("the installed artifact is not the published one")
	}
	// Staging is cleaned up: nothing but the artifact is left behind.
	entries, err := os.ReadDir(p.destDir)
	if err != nil {
		t.Fatal(err)
	}
	if len(entries) != 1 || entries[0].Name() != filepath.Base(p.destPath) {
		t.Errorf("destination holds %v, want only the installed artifact", entries)
	}
}

// A chunked artifact must reassemble to exactly what was pushed, through the
// real push and pull paths on both ends.
func TestPullReassemblesAChunkedArtifact(t *testing.T) {
	requireSquashfsTools(t)

	p := newPullFixture(t)
	if len(p.registry.blobs) < 3 {
		t.Fatalf("only %d blobs published; the artifact did not chunk", len(p.registry.blobs))
	}
	if err := p.pull(context.Background()); err != nil {
		t.Fatalf("Pull: %v", err)
	}
	if !sameBytes(t, p.source, p.destPath) {
		t.Error("the reassembled artifact differs from the published one")
	}
}

// A pinned artifact is never replaced, and the reason given is that it is pinned
// — not a transport failure the reader would go looking for a network problem in.
func TestPullRefusesAProtectedDestination(t *testing.T) {
	p := newPullFixture(t)
	if err := os.WriteFile(p.destPath, []byte("pinned"), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.Chmod(p.destPath, 0o444); err != nil {
		t.Fatal(err)
	}
	t.Cleanup(func() { os.Chmod(p.destPath, 0o644) }) //nolint:errcheck

	err := p.pull(context.Background())
	if !errors.Is(err, image.ErrProtected) {
		t.Fatalf("err = %v, want image.ErrProtected", err)
	}
	got, readErr := os.ReadFile(p.destPath)
	if readErr != nil {
		t.Fatal(readErr)
	}
	if string(got) != "pinned" {
		t.Error("a pinned artifact was replaced anyway")
	}
}

// The availability probe belongs after the download: a lock held when the pull
// starts but released while it runs must not refuse an install that is fine by
// the time it matters.
func TestPullProbesAfterTheDownload(t *testing.T) {
	p := newPullFixture(t)
	if err := os.WriteFile(p.destPath, []byte("old"), 0o644); err != nil {
		t.Fatal(err)
	}
	lock, err := image.AcquireLock(p.destPath, false)
	if err != nil {
		t.Fatal(err)
	}
	// Guarded: the release happens on the server's goroutine and is read back here.
	var mu sync.Mutex
	released := false
	p.registry.onBlob = func() {
		mu.Lock()
		defer mu.Unlock()
		if !released {
			released = true
			lock.Close() //nolint:errcheck
		}
	}

	if err := p.pull(context.Background()); err != nil {
		t.Fatalf("Pull refused an install that was locked only before the download: %v", err)
	}
	mu.Lock()
	defer mu.Unlock()
	if !released {
		t.Fatal("no blob was served, so the test proved nothing")
	}
	if !sameBytes(t, p.source, p.destPath) {
		t.Error("the artifact was not installed")
	}
}

// The mirror image: a reader that arrives *during* the download must still be
// caught, which is only possible if the probe runs after it.
func TestPullRefusesADestinationLockedDuringTheDownload(t *testing.T) {
	p := newPullFixture(t)
	if err := os.WriteFile(p.destPath, []byte("in use"), 0o644); err != nil {
		t.Fatal(err)
	}
	// The lock is taken on the server's goroutine and released on this one, so
	// the handoff is guarded: an unsynchronized pointer can read back nil here
	// and leave the file open, which on NFS leaves an undeletable silly-rename.
	var mu sync.Mutex
	var lock *image.Lock
	p.registry.onBlob = func() {
		mu.Lock()
		defer mu.Unlock()
		if lock == nil {
			lock, _ = image.AcquireLock(p.destPath, false)
		}
	}
	release := func() {
		mu.Lock()
		defer mu.Unlock()
		if lock != nil {
			lock.Close() //nolint:errcheck
			lock = nil
		}
	}
	t.Cleanup(release)

	err := p.pull(context.Background())
	if !errors.Is(err, image.ErrInUse) {
		t.Fatalf("err = %v, want image.ErrInUse", err)
	}
	release()
	got, readErr := os.ReadFile(p.destPath)
	if readErr != nil {
		t.Fatal(readErr)
	}
	if string(got) != "in use" {
		t.Error("an artifact in use was replaced anyway")
	}
}

// A base served where an overlay was asked for fails at the transport, not at
// mount. The destination's extension is what states the expectation.
func TestPullRefusesTheWrongArtifactType(t *testing.T) {
	requireSquashfsTools(t)
	source, m := packImage(t, imageSpec{name: "hello/1.0", recipe: "#!/bin/bash\n"})
	f := newFakeRegistry(t)
	ann := Annotations(m, "")
	desc := publishImage(t, f, "hello", source, ann, ArtifactTypeBase, MediaTypeBaseBlob, "1.0")

	dest := filepath.Join(t.TempDir(), "hello--1.0.sqf")
	err := Pull(context.Background(), f.base(), "hello", desc, ann, dest)
	if !errors.Is(err, ErrInvalidArtifact) {
		t.Fatalf("err = %v, want ErrInvalidArtifact", err)
	}
	if _, statErr := os.Stat(dest); statErr == nil {
		t.Error("a wrong-typed artifact was installed anyway")
	}
}

// The payload's own files regenerate the keys; annotations are only a claim
// about them. A publisher whose two disagree is refused after the transfer.
func TestPullRefusesAPayloadThatContradictsItsAnnotations(t *testing.T) {
	p := newPullFixture(t)
	p.ann[AnnIdentitySHA] = strings.Repeat("f", 64)

	err := p.pull(context.Background())
	if !errors.Is(err, ErrInvalidArtifact) {
		t.Fatalf("err = %v, want ErrInvalidArtifact", err)
	}
	if _, statErr := os.Stat(p.destPath); statErr == nil {
		t.Error("an incoherent artifact was installed anyway")
	}
}

// An artifact installed into a shared directory must be group-writable, so the
// next person in the lab can replace it. Staging is 0700; the install must not
// inherit that.
func TestPullInstallsGroupWritableIntoASharedDirectory(t *testing.T) {
	p := newPullFixture(t)
	if err := os.Chmod(p.destDir, 0o2775); err != nil {
		t.Skipf("cannot set the shared mode on %s: %v", p.destDir, err)
	}
	if err := p.pull(context.Background()); err != nil {
		t.Fatalf("Pull: %v", err)
	}
	info, err := os.Stat(p.destPath)
	if err != nil {
		t.Fatal(err)
	}
	if info.Mode().Perm()&0o060 != 0o060 {
		t.Errorf("installed mode is %v, want group read-write in a shared directory", info.Mode().Perm())
	}
}

func TestImageTypes(t *testing.T) {
	tests := []struct {
		path         string
		artifactType string
		wantErr      bool
	}{
		{"/images/hello--1.0.sqf", ArtifactTypeOverlay, false},
		{"/images/ubuntu24.sif", ArtifactTypeBase, false},
		{"/images/dev.img", "", true},
		{"/images/hello", "", true},
	}
	for _, tt := range tests {
		got, _, err := imageTypes(tt.path)
		if (err != nil) != tt.wantErr {
			t.Errorf("imageTypes(%q) err = %v, wantErr %v", tt.path, err, tt.wantErr)
		}
		if got != tt.artifactType {
			t.Errorf("imageTypes(%q) = %q, want %q", tt.path, got, tt.artifactType)
		}
	}
}

func TestDownloadSize(t *testing.T) {
	desc := ocispec.Descriptor{Size: 100}
	single := ocispec.Manifest{
		Config: ocispec.Descriptor{Size: 10},
		Layers: []ocispec.Descriptor{{Size: 1000}},
	}
	if got, want := downloadSize(desc, single), int64(1110); got != want {
		t.Errorf("downloadSize = %d, want %d", got, want)
	}

	// A chunked artifact costs no more than the same bytes in one layer: each
	// one streams into its own range of the finished file, so the payload never
	// exists twice. This used to be doubled, and for 40 GB that was 40 GB of
	// scratch nobody needed.
	chunked := ocispec.Manifest{
		Config: ocispec.Descriptor{Size: 10},
		Layers: []ocispec.Descriptor{{Size: 600}, {Size: 400}},
	}
	if got, want := downloadSize(desc, chunked), int64(1110); got != want {
		t.Errorf("downloadSize = %d, want %d — chunking must not cost a second copy", got, want)
	}

	// A descriptor with no size recorded must not subtract from the total.
	negative := ocispec.Manifest{Layers: []ocispec.Descriptor{{Size: -1}, {Size: 500}}}
	if got := downloadSize(ocispec.Descriptor{Size: -1}, negative); got != 500 {
		t.Errorf("downloadSize = %d, want 500", got)
	}
}

func TestRequireFreeSpace(t *testing.T) {
	dir := t.TempDir()
	if err := requireFreeSpace(dir, 0); err != nil {
		t.Errorf("a zero requirement was refused: %v", err)
	}
	if err := requireFreeSpace(dir, 1); err != nil {
		t.Errorf("one byte was refused: %v", err)
	}
	// An exabyte fits nowhere.
	err := requireFreeSpace(dir, 1<<60)
	if err == nil {
		t.Fatal("an impossible requirement was accepted")
	}
	if !strings.Contains(err.Error(), dir) {
		t.Errorf("the refusal does not name the directory: %v", err)
	}
	if err := requireFreeSpace(filepath.Join(dir, "absent"), 1); err == nil {
		t.Error("a missing directory was accepted")
	}
}

// A pull runs once per node per artifact, so its request rate is the larger of
// the two directions across a cluster. Whatever the number is, it is chosen.
func TestPullConcurrencyIsStated(t *testing.T) {
	// ORAS leaves Concurrency zero and fills in its own default at copy time, so
	// a positive value here is the difference between choosing the number and
	// inheriting whichever one the library happens to ship.
	if oras.DefaultCopyOptions.Concurrency != 0 {
		t.Fatalf("ORAS now defaults Concurrency to %d; the reasoning below needs rechecking",
			oras.DefaultCopyOptions.Concurrency)
	}
	if pullConcurrency <= 0 {
		t.Errorf("pullConcurrency = %d, which leaves the choice to ORAS", pullConcurrency)
	}
}

// The whole point of writing at offsets: a 40 GB artifact used to need 80 GB of
// scratch, because every layer was staged as its own file and then joined into a
// third. Nothing but the finished artifact may exist while a pull runs.
func TestPullStagesNoSecondCopy(t *testing.T) {
	p := newPullFixture(t)

	var peak int
	done := make(chan struct{})
	var wg sync.WaitGroup
	wg.Add(1)
	go func() {
		defer wg.Done()
		for {
			select {
			case <-done:
				return
			default:
			}
			if entries, err := os.ReadDir(p.destDir); err == nil {
				for _, e := range entries {
					if !e.IsDir() {
						continue
					}
					staged, err := os.ReadDir(filepath.Join(p.destDir, e.Name()))
					if err == nil && len(staged) > peak {
						peak = len(staged)
					}
				}
			}
		}
	}()

	if err := p.pull(context.Background()); err != nil {
		t.Fatalf("Pull: %v", err)
	}
	close(done)
	wg.Wait()

	if peak > 1 {
		t.Errorf("staging held %d files at once; a layer was written somewhere other than "+
			"its own range of the artifact", peak)
	}
}

// Layers arrive concurrently and land at offsets, so the finished file has to be
// the exact concatenation whatever order they completed in.
func TestPullReassemblesUnderConcurrency(t *testing.T) {
	if pullConcurrency < 2 {
		t.Skip("pull is serial; there is no ordering to disturb")
	}
	p := newPullFixture(t)
	if len(p.registry.blobs) < 3 {
		t.Fatalf("only %d blobs published; the artifact did not chunk", len(p.registry.blobs))
	}
	if err := p.pull(context.Background()); err != nil {
		t.Fatalf("Pull: %v", err)
	}
	if !sameBytes(t, p.source, p.destPath) {
		t.Error("concurrent offset writes did not reproduce the artifact")
	}
}

// A layer whose bytes are not what the manifest promised must never reach the
// installed file. ORAS checks Content-Length and the Docker-Content-Digest
// header on a fetch but never hashes the body, so this path verifies it itself.
func TestDownloadRejectsACorruptedLayer(t *testing.T) {
	p := newPullFixture(t)

	// Specifically one of this manifest's layers. Picking any large blob drew
	// from map order, and the index blob is never fetched — Pull is handed the
	// child manifest — so corrupting that one left the pull free to succeed.
	var manifest ocispec.Manifest
	if err := json.Unmarshal(p.registry.blobs[p.desc.Digest], &manifest); err != nil {
		t.Fatalf("read the published manifest: %v", err)
	}
	if len(manifest.Layers) == 0 {
		t.Fatal("the published manifest carries no layers")
	}
	target := manifest.Layers[0].Digest
	corrupted := append([]byte(nil), p.registry.blobs[target]...)
	if len(corrupted) == 0 {
		t.Fatalf("layer %s was not published", target)
	}
	corrupted[len(corrupted)/2] ^= 0xff
	p.registry.blobs[target] = corrupted

	if err := p.pull(context.Background()); err == nil {
		t.Fatal("a layer whose bytes contradict its digest was installed")
	}
	if _, err := os.Stat(p.destPath); !os.IsNotExist(err) {
		t.Error("a corrupted download left a file at the destination")
	}
}

// flakyBlobStore serves a blob, refusing the first refusals attempts the way a
// throttled registry does. Only Fetch is implemented; the download path calls
// nothing else.
type flakyBlobStore struct {
	registry.BlobStore
	content  []byte
	refusals int
	fetches  int
	// truncate cuts the body short on the first attempt, standing in for a
	// connection lost partway through a layer.
	truncate bool
}

func (f *flakyBlobStore) Fetch(_ context.Context, _ ocispec.Descriptor) (io.ReadCloser, error) {
	f.fetches++
	if f.fetches <= f.refusals {
		if !f.truncate {
			return nil, codedErr(429, "TOOMANYREQUESTS", "slow down")
		}
		return io.NopCloser(bytes.NewReader(f.content[:len(f.content)/2])), nil
	}
	return io.NopCloser(bytes.NewReader(f.content)), nil
}

func layerOf(content []byte) ocispec.Descriptor {
	return ocispec.Descriptor{
		MediaType: MediaTypeOverlayBlob,
		Digest:    digest.FromBytes(content),
		Size:      int64(len(content)),
	}
}

// A rate limit on the way down is waited out, not failed. It cannot fall back to
// a local build either — ErrRateLimited is deliberately not ErrUnavailable — so
// waiting is the only thing that turns it into a working pull.
func TestDownloadWaitsOutARateLimit(t *testing.T) {
	content := []byte(strings.Repeat("payload", 64))
	blobs := &flakyBlobStore{content: content, refusals: 2}

	path := filepath.Join(t.TempDir(), "hello--1.0.sqf")
	f, err := os.Create(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close() //nolint:errcheck

	policy := newTestPolicy()
	ctx := withRetryPolicy(context.Background(), policy.retryPolicy)
	progress := newDownloadProgress(ctx, int64(len(content)))

	if err := fetchLayerAt(ctx, blobs, layerOf(content), f, 0, progress); err != nil {
		t.Fatalf("fetchLayerAt: %v", err)
	}
	if blobs.fetches != 3 {
		t.Errorf("fetches = %d, want the two refusals retried", blobs.fetches)
	}
	if len(*policy.waits) != 2 {
		t.Errorf("waits = %v, want one per refusal", *policy.waits)
	}
	got, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	if string(got) != string(content) {
		t.Error("the retried layer did not land intact")
	}
}

// What writing at a fixed offset buys that a push cannot have: a layer lost
// partway through is simply written again over the same range, so a partial
// transfer needs no resume protocol to recover from.
func TestDownloadRewritesAPartialLayer(t *testing.T) {
	content := []byte(strings.Repeat("payload", 64))
	blobs := &flakyBlobStore{content: content, refusals: 1, truncate: true}

	path := filepath.Join(t.TempDir(), "hello--1.0.sqf")
	f, err := os.Create(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close() //nolint:errcheck

	// A truncated body fails its digest, which is what makes the attempt
	// retryable rather than silently short.
	progress := newDownloadProgress(context.Background(), int64(len(content)))
	err = fetchLayerAt(context.Background(), blobs, layerOf(content), f, 0, progress)
	if err == nil {
		t.Fatal("a truncated layer was accepted")
	}
	if progress.done != 0 {
		t.Errorf("done = %d after a failed attempt, want the bytes withdrawn", progress.done)
	}
}
