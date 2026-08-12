package registry

import (
	"context"
	"encoding/json"
	"errors"
	"fmt"
	"io"
	"net/http"
	"net/http/httptest"
	"slices"
	"strings"
	"sync"
	"testing"

	"github.com/opencontainers/go-digest"
	specs "github.com/opencontainers/image-spec/specs-go"
	ocispec "github.com/opencontainers/image-spec/specs-go/v1"
)

// fakeRegistry serves just enough of the distribution API to resolve a manifest
// and list tags. It exists so the classification every fallback decision depends
// on is exercised through the real oras client rather than asserted on a mock.
type fakeRegistry struct {
	server *httptest.Server
	// mu guards the maps: push and pull both run requests on server goroutines.
	mu sync.Mutex
	// content is keyed "<repo>/<tag-or-digest>".
	content map[string]ocispec.Descriptor
	blobs   map[digest.Digest][]byte
	tags    map[string][]string
	// status, when set, is returned for every request, standing in for a
	// registry that is down, private, or empty.
	status int
	// onBlob, when set, runs just before a payload blob is served. It is the
	// seam for a test that needs something to happen *during* a transfer.
	onBlob func()
}

func newFakeRegistry(t *testing.T) *fakeRegistry {
	t.Helper()
	f := &fakeRegistry{
		content: map[string]ocispec.Descriptor{},
		blobs:   map[digest.Digest][]byte{},
		tags:    map[string][]string{},
	}
	f.server = httptest.NewServer(http.HandlerFunc(f.serve))
	t.Cleanup(f.server.Close)
	return f
}

// base is the registry host, which is loopback, so newRepository speaks plain
// HTTP to it and no certificate is involved.
func (f *fakeRegistry) base() string { return strings.TrimPrefix(f.server.URL, "http://") }

func (f *fakeRegistry) serve(w http.ResponseWriter, r *http.Request) {
	f.mu.Lock()
	defer f.mu.Unlock()

	if f.status != 0 {
		// No Www-Authenticate challenge: there is no token endpoint to retry against.
		w.WriteHeader(f.status)
		return
	}
	path := strings.TrimPrefix(r.URL.Path, "/v2/")

	if r.Method == http.MethodPost || r.Method == http.MethodPut || r.Method == http.MethodPatch {
		f.accept(w, r, path)
		return
	}

	if repo, ok := strings.CutSuffix(path, "/tags/list"); ok {
		tags, known := f.tags[repo]
		if !known {
			w.WriteHeader(http.StatusNotFound)
			return
		}
		w.Header().Set("Content-Type", "application/json")
		json.NewEncoder(w).Encode(map[string]any{"name": repo, "tags": tags}) //nolint:errcheck
		return
	}

	if _, ref, ok := cutLast(path, "/blobs/"); ok {
		body, found := f.blobs[digest.Digest(ref)]
		if !found {
			w.WriteHeader(http.StatusNotFound)
			return
		}
		if f.onBlob != nil {
			f.onBlob()
		}
		f.write(w, r, "application/octet-stream", ref, body)
		return
	}

	repo, ref, ok := cutLast(path, "/manifests/")
	if !ok {
		w.WriteHeader(http.StatusNotFound)
		return
	}
	desc, found := f.content[repo+"/"+ref]
	if !found {
		w.WriteHeader(http.StatusNotFound)
		return
	}
	f.write(w, r, desc.MediaType, desc.Digest.String(), f.blobs[desc.Digest])
}

// accept handles the write half of the distribution API: a blob upload session
// and a manifest PUT. Enough for a real push to complete against it, so a push
// test exercises the protocol rather than a mock of it.
func (f *fakeRegistry) accept(w http.ResponseWriter, r *http.Request, path string) {
	if repo, _, ok := cutLast(path, "/blobs/uploads/"); ok {
		// POST opens a session; the PUT that follows carries the bytes.
		if r.Method == http.MethodPost && r.URL.Query().Get("digest") == "" {
			w.Header().Set("Location", "/v2/"+repo+"/blobs/uploads/session")
			w.WriteHeader(http.StatusAccepted)
			return
		}
		body, err := io.ReadAll(r.Body)
		if err != nil {
			http.Error(w, err.Error(), http.StatusInternalServerError)
			return
		}
		dgst := digest.FromBytes(body)
		if want := r.URL.Query().Get("digest"); want != "" && want != dgst.String() {
			http.Error(w, "digest mismatch: uploaded "+dgst.String(), http.StatusBadRequest)
			return
		}
		f.blobs[dgst] = body
		w.Header().Set("Docker-Content-Digest", dgst.String())
		w.WriteHeader(http.StatusCreated)
		return
	}

	if repo, ref, ok := cutLast(path, "/manifests/"); ok {
		body, err := io.ReadAll(r.Body)
		if err != nil {
			http.Error(w, err.Error(), http.StatusInternalServerError)
			return
		}
		desc := ocispec.Descriptor{
			MediaType: r.Header.Get("Content-Type"),
			Digest:    digest.FromBytes(body),
			Size:      int64(len(body)),
		}
		f.blobs[desc.Digest] = body
		f.content[repo+"/"+ref] = desc
		f.content[repo+"/"+desc.Digest.String()] = desc
		if !strings.HasPrefix(ref, digestPrefix) && !slices.Contains(f.tags[repo], ref) {
			f.tags[repo] = append(f.tags[repo], ref)
		}
		w.Header().Set("Docker-Content-Digest", desc.Digest.String())
		w.WriteHeader(http.StatusCreated)
		return
	}
	w.WriteHeader(http.StatusNotFound)
}

func (f *fakeRegistry) write(w http.ResponseWriter, r *http.Request, mediaType, dgst string, body []byte) {
	w.Header().Set("Content-Type", mediaType)
	w.Header().Set("Docker-Content-Digest", dgst)
	w.Header().Set("Content-Length", fmt.Sprint(len(body)))
	if r.Method == http.MethodHead {
		w.WriteHeader(http.StatusOK)
		return
	}
	w.Write(body) //nolint:errcheck
}

// blobStore lets a test publish payloads through the real chunked push, so a
// pull test is a round trip rather than an assertion about a hand-built layout.
type fakeBlobs struct{ f *fakeRegistry }

func (b fakeBlobs) Exists(_ context.Context, target ocispec.Descriptor) (bool, error) {
	_, ok := b.f.blobs[target.Digest]
	return ok, nil
}

func (b fakeBlobs) Push(_ context.Context, expected ocispec.Descriptor, content io.Reader) error {
	data, err := io.ReadAll(content)
	if err != nil {
		return err
	}
	b.f.blobs[expected.Digest] = data
	return nil
}

// cutLast splits around the final occurrence of sep, so a repository path with
// slashes in it stays intact.
func cutLast(s, sep string) (before, after string, found bool) {
	i := strings.LastIndex(s, sep)
	if i < 0 {
		return s, "", false
	}
	return s[:i], s[i+len(sep):], true
}

// publish stores content under its digest and every given reference.
func (f *fakeRegistry) publish(t *testing.T, repo, mediaType string, payload any, refs ...string) ocispec.Descriptor {
	t.Helper()
	body, err := json.Marshal(payload)
	if err != nil {
		t.Fatal(err)
	}
	desc := ocispec.Descriptor{
		MediaType: mediaType,
		Digest:    digest.FromBytes(body),
		Size:      int64(len(body)),
	}
	f.blobs[desc.Digest] = body
	for _, ref := range append(refs, desc.Digest.String()) {
		f.content[repo+"/"+ref] = desc
	}
	return desc
}

// publishArtifact stores a per-platform manifest under an index, the shape a
// native push produces.
func publishArtifact(t *testing.T, f *fakeRegistry, repo string, plat ocispec.Platform, ann map[string]string, tags ...string) ocispec.Descriptor {
	t.Helper()
	child := f.publish(t, repo, ocispec.MediaTypeImageManifest, ocispec.Manifest{
		Versioned:    specs.Versioned{SchemaVersion: 2},
		MediaType:    ocispec.MediaTypeImageManifest,
		ArtifactType: ArtifactTypeOverlay,
		Config:       ocispec.DescriptorEmptyJSON,
		Layers:       []ocispec.Descriptor{{MediaType: MediaTypeOverlayBlob, Size: 10}},
		Annotations:  ann,
	})
	child.Platform = &plat
	f.publish(t, repo, ocispec.MediaTypeImageIndex, ocispec.Index{
		Versioned:   specs.Versioned{SchemaVersion: 2},
		MediaType:   ocispec.MediaTypeImageIndex,
		Manifests:   []ocispec.Descriptor{child},
		Annotations: indexAnnotations(ann),
	}, tags...)
	return child
}

func nativePlatform(t *testing.T) ocispec.Platform {
	t.Helper()
	plat, ok := platform()
	if !ok {
		t.Skip("not a distribution architecture")
	}
	return plat
}

func TestResolveArtifactDescendsTheIndex(t *testing.T) {
	f := newFakeRegistry(t)
	want := Annotations(fullManifest(), "zstd")
	child := publishArtifact(t, f, "grch38/star/2.7.11b", nativePlatform(t), want, "gencode49-101")

	desc, ann, err := ResolveArtifact(context.Background(), f.base(), "grch38/star/2.7.11b", "gencode49-101")
	if err != nil {
		t.Fatalf("ResolveArtifact: %v", err)
	}
	if desc.Digest != child.Digest {
		t.Errorf("resolved %s, want this platform's child %s", desc.Digest, child.Digest)
	}
	if desc.ArtifactType != ArtifactTypeOverlay {
		t.Errorf("artifact type = %q, want %q", desc.ArtifactType, ArtifactTypeOverlay)
	}
	if ann[AnnIdentitySHA] != want[AnnIdentitySHA] || ann[AnnCompression] != "zstd" {
		t.Errorf("annotations = %v, want the child's own", ann)
	}
	// The resolved descriptor must be usable as an address on its own.
	if _, _, err := ResolveArtifact(context.Background(), f.base(), "grch38/star/2.7.11b", desc.Digest.String()); err != nil {
		t.Errorf("resolving the returned digest: %v", err)
	}
}

// A noarch artifact is tagged as a bare manifest with no index above it.
func TestResolveArtifactAcceptsABareManifest(t *testing.T) {
	f := newFakeRegistry(t)
	m := fullManifest()
	m.Platform.Arch = "noarch"
	want := Annotations(m, "zstd")
	f.publish(t, "hello", ocispec.MediaTypeImageManifest, ocispec.Manifest{
		Versioned:    specs.Versioned{SchemaVersion: 2},
		MediaType:    ocispec.MediaTypeImageManifest,
		ArtifactType: ArtifactTypeOverlay,
		Config:       ocispec.DescriptorEmptyJSON,
		Annotations:  want,
	}, "1.0")

	desc, ann, err := ResolveArtifact(context.Background(), f.base(), "hello", "1.0")
	if err != nil {
		t.Fatalf("ResolveArtifact: %v", err)
	}
	if !IsNoarch(ann) {
		t.Errorf("annotations = %v, want the noarch mark", ann)
	}
	if desc.ArtifactType != ArtifactTypeOverlay {
		t.Errorf("artifact type = %q, want %q", desc.ArtifactType, ArtifactTypeOverlay)
	}
}

// The artifact exists but nobody published this architecture. That is a
// different answer from "no such artifact", and the message must name what was
// published so the reader knows a push is missing rather than a name wrong.
func TestResolveArtifactReportsAMissingPlatform(t *testing.T) {
	f := newFakeRegistry(t)
	native := nativePlatform(t)
	foreign := ocispec.Platform{OS: "linux", Architecture: "s390x"}
	if native.Architecture == foreign.Architecture {
		t.Skip("running on the architecture used as the foreign one")
	}
	publishArtifact(t, f, "cellranger", foreign, Annotations(fullManifest(), ""), "9.0.1")

	_, _, err := ResolveArtifact(context.Background(), f.base(), "cellranger", "9.0.1")
	if !errors.Is(err, ErrUnsupportedPlatform) {
		t.Fatalf("err = %v, want ErrUnsupportedPlatform", err)
	}
	if !strings.Contains(err.Error(), "linux/s390x") {
		t.Errorf("the refusal does not name what was published: %v", err)
	}
}

func TestResolveArtifactNotFound(t *testing.T) {
	f := newFakeRegistry(t)
	_, _, err := ResolveArtifact(context.Background(), f.base(), "cellranger", "9.0.1")
	if !errors.Is(err, ErrNotFound) {
		t.Fatalf("err = %v, want ErrNotFound", err)
	}
}

// An expired token must never read as "not published": that is the difference
// between an error and a silent forty-minute rebuild.
func TestResolveArtifactUnauthorized(t *testing.T) {
	f := newFakeRegistry(t)
	f.status = http.StatusUnauthorized

	_, _, err := ResolveArtifact(context.Background(), f.base(), "cellranger", "9.0.1")
	if !errors.Is(err, ErrUnauthorized) {
		t.Fatalf("err = %v, want ErrUnauthorized", err)
	}
	if errors.Is(err, ErrNotFound) {
		t.Error("an authorization failure was also reported as not published")
	}
}

func TestExists(t *testing.T) {
	f := newFakeRegistry(t)
	publishArtifact(t, f, "cellranger", nativePlatform(t), Annotations(fullManifest(), ""), "9.0.1")
	ctx := context.Background()

	if ok, err := Exists(ctx, f.base(), "cellranger", "9.0.1"); err != nil || !ok {
		t.Errorf("Exists = (%v, %v), want (true, nil)", ok, err)
	}
	if ok, err := Exists(ctx, f.base(), "cellranger", "8.0.0"); err != nil || ok {
		t.Errorf("Exists = (%v, %v), want (false, nil)", ok, err)
	}
}

// "Nothing is published here" and "the registry would not say" are different
// answers, and push overwrites nothing only on the first.
func TestExistsReportsAnUnreadableRegistryRatherThanFalse(t *testing.T) {
	f := newFakeRegistry(t)
	f.status = http.StatusUnauthorized

	ok, err := Exists(context.Background(), f.base(), "cellranger", "9.0.1")
	if ok {
		t.Error("an unreadable registry reported the tag as present")
	}
	if !errors.Is(err, ErrUnauthorized) {
		t.Errorf("err = %v, want ErrUnauthorized rather than a bare false", err)
	}
}

func TestListTags(t *testing.T) {
	f := newFakeRegistry(t)
	f.tags["ubuntu24/base"] = []string{"20260721", "20260801", RollingTag}

	tags, err := ListTags(context.Background(), f.base(), "ubuntu24/base")
	if err != nil {
		t.Fatalf("ListTags: %v", err)
	}
	for _, want := range []string{"20260721", "20260801", RollingTag} {
		if !slices.Contains(tags, want) {
			t.Errorf("tags = %v, missing %q", tags, want)
		}
	}
}

// A name nobody has pushed yet is the ordinary state of a name, not a failure.
func TestListTagsOnAnUnpublishedRepository(t *testing.T) {
	f := newFakeRegistry(t)
	tags, err := ListTags(context.Background(), f.base(), "cellranger")
	if err != nil {
		t.Fatalf("ListTags: %v", err)
	}
	if len(tags) != 0 {
		t.Errorf("tags = %v, want none", tags)
	}
}

func TestListTagsReportsAnUnreadableRegistry(t *testing.T) {
	f := newFakeRegistry(t)
	f.status = http.StatusUnauthorized

	if _, err := ListTags(context.Background(), f.base(), "cellranger"); !errors.Is(err, ErrUnauthorized) {
		t.Fatalf("err = %v, want ErrUnauthorized", err)
	}
}
