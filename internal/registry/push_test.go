package registry

import (
	"context"
	"encoding/json"
	"errors"
	"path/filepath"
	"slices"
	"strings"
	"testing"

	ocispec "github.com/opencontainers/image-spec/specs-go/v1"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
)

func TestVisibilityAccepts(t *testing.T) {
	manifest := func(typ catalog.Type, buildType string) meta.Manifest {
		return meta.Manifest{Name: "x/1.0", Type: typ, BuildType: buildType}
	}
	tests := []struct {
		why        string
		visibility Visibility
		m          meta.Manifest
		wantErr    bool
	}{
		{"a base is ours to publish", Public, manifest(catalog.TypeBase, "def"), false},
		{"apt packages from a public distribution", Public, manifest(catalog.TypeOS, "script"), false},
		{"public reference data", Public, manifest(catalog.TypeData, "script"), false},
		{"an app redistributes someone else's binaries", Public, manifest(catalog.TypeApp, "script"), true},
		{"a Conda build embeds its own solve inputs", Public, manifest(catalog.TypeApp, "conda"), true},
		// Independent reasons: a Conda build is refused even where the type is fine.
		{"a Conda build of a permitted type", Public, manifest(catalog.TypeData, "conda"), true},
		{"an internal endpoint takes an app", Internal, manifest(catalog.TypeApp, "script"), false},
		{"an internal endpoint takes a Conda build", Internal, manifest(catalog.TypeApp, "conda"), false},
	}
	for _, tt := range tests {
		if err := tt.visibility.Accepts(tt.m); (err != nil) != tt.wantErr {
			t.Errorf("%s: %s.Accepts(%s/%s) = %v, wantErr %v",
				tt.why, tt.visibility, tt.m.Type, tt.m.BuildType, err, tt.wantErr)
		}
	}
}

// The default is the restrictive one: an endpoint that never declared its
// visibility must not be treated as private.
func TestPublishDefaultsToPublic(t *testing.T) {
	requireSquashfsTools(t)
	source, _ := packImage(t, imageSpec{name: "hello/1.0", typ: catalog.TypeApp, recipe: "#!/bin/bash\n"})
	f := newFakeRegistry(t)

	err := Publish(context.Background(), PublishRequest{Path: source, Base: f.base()})
	if err == nil {
		t.Fatal("an app was published to an endpoint with no declared visibility")
	}
	if !strings.Contains(err.Error(), string(Public)) {
		t.Errorf("the refusal does not say the endpoint was treated as public: %v", err)
	}
}

// The end-to-end contract: what push writes is what pull installs, byte for byte,
// through the real distribution protocol on both sides.
func TestPublishAndPullRoundTrip(t *testing.T) {
	requireSquashfsTools(t)
	source, m := packImage(t, imageSpec{
		name: "grch38/genome/gencode49", typ: catalog.TypeData, recipe: "#!/bin/bash\nbuild index\n",
	})
	f := newFakeRegistry(t)
	ctx := context.Background()

	if err := Publish(ctx, PublishRequest{Path: source, Base: f.base()}); err != nil {
		t.Fatalf("Publish: %v", err)
	}

	repo, tag, err := PullReference(m.Type, m.Name)
	if err != nil {
		t.Fatal(err)
	}
	desc, ann, err := ResolveArtifact(ctx, f.base(), repo, tag)
	if err != nil {
		t.Fatalf("ResolveArtifact: %v", err)
	}
	if err := Check(ann, Want{Name: m.Name, Identity: m.Keys.Identity}); err != nil {
		t.Fatalf("Check: %v", err)
	}

	dest := filepath.Join(t.TempDir(), "grch38-genome--gencode49.sqf")
	if err := Pull(ctx, f.base(), repo, desc, ann, dest); err != nil {
		t.Fatalf("Pull: %v", err)
	}
	if !sameBytes(t, source, dest) {
		t.Error("the pulled artifact is not the one that was published")
	}
}

// Nothing can be pinned to an artifact with no scheme-backed keys, so publishing
// one would put something in a registry that no lock could ever name.
func TestPublishRefusesAnArtifactWithNoKeys(t *testing.T) {
	requireSquashfsTools(t)
	source, _ := packImage(t, imageSpec{
		name: "ubuntu24/build-essential", typ: catalog.TypeOS, recipe: "#!/bin/bash\n", noKeys: true,
	})
	f := newFakeRegistry(t)

	err := Publish(context.Background(), PublishRequest{Path: source, Base: f.base()})
	if err == nil {
		t.Fatal("an artifact with no keys was published")
	}
	if !strings.Contains(err.Error(), "identity") {
		t.Errorf("the refusal does not name the missing identity: %v", err)
	}
}

// An artifact is held to its own claim before anyone downstream has to trust it.
func TestPublishRefusesKeysThatDoNotRegenerate(t *testing.T) {
	requireSquashfsTools(t)
	source, _ := packImage(t, imageSpec{
		name: "ubuntu24/build-essential", typ: catalog.TypeOS, recipe: "#!/bin/bash\n",
		tamper: func(m *meta.Manifest) { m.Keys.Identity.SHA256 = strings.Repeat("f", 64) },
	})
	f := newFakeRegistry(t)

	err := Publish(context.Background(), PublishRequest{Path: source, Base: f.base()})
	if !errors.Is(err, ErrInvalidArtifact) {
		t.Fatalf("err = %v, want ErrInvalidArtifact", err)
	}
	if len(f.tags) != 0 {
		t.Errorf("a refused artifact was published anyway: %v", f.tags)
	}
}

// A versioned tag is immutable, because the one thing a version promises is that
// it does not change under someone who already pulled it.
func TestPublishRefusesAnExistingVersionedTag(t *testing.T) {
	requireSquashfsTools(t)
	source, _ := packImage(t, imageSpec{
		name: "grch38/genome/gencode49", typ: catalog.TypeData, recipe: "#!/bin/bash\n",
	})
	f := newFakeRegistry(t)
	ctx := context.Background()

	if err := Publish(ctx, PublishRequest{Path: source, Base: f.base()}); err != nil {
		t.Fatalf("first publish: %v", err)
	}
	err := Publish(ctx, PublishRequest{Path: source, Base: f.base()})
	if err == nil {
		t.Fatal("a versioned tag was silently replaced")
	}
	if !strings.Contains(err.Error(), "--force") {
		t.Errorf("the refusal does not name the way through: %v", err)
	}
	if err := Publish(ctx, PublishRequest{Path: source, Base: f.base(), Force: true}); err != nil {
		t.Errorf("--force did not permit the replacement: %v", err)
	}
}

// A version-less artifact is addressed by build date and a rolling tag;
// republishing it is the point of that scheme, so nothing is protected.
func TestPublishVersionLessRepublishesUnderBothTags(t *testing.T) {
	requireSquashfsTools(t)
	source, m := packImage(t, imageSpec{
		name: "ubuntu24/build-essential", typ: catalog.TypeOS, recipe: "#!/bin/bash\n",
	})
	f := newFakeRegistry(t)
	ctx := context.Background()

	if err := Publish(ctx, PublishRequest{Path: source, Base: f.base()}); err != nil {
		t.Fatalf("Publish: %v", err)
	}
	if err := Publish(ctx, PublishRequest{Path: source, Base: f.base()}); err != nil {
		t.Errorf("a version-less artifact could not be republished: %v", err)
	}

	tags, err := ListTags(ctx, f.base(), "ubuntu24/build-essential")
	if err != nil {
		t.Fatal(err)
	}
	dateTag := m.Build.Created.UTC().Format(dateTagLayout)
	for _, want := range []string{dateTag, RollingTag} {
		if !slices.Contains(tags, want) {
			t.Errorf("tags = %v, missing %q", tags, want)
		}
	}
}

// A native artifact is published under an index, so a second architecture can add
// itself to the same tag.
func TestPublishWrapsANativeArtifactInAnIndex(t *testing.T) {
	requireSquashfsTools(t)
	source, m := packImage(t, imageSpec{
		name: "grch38/genome/gencode49", typ: catalog.TypeData, recipe: "#!/bin/bash\n",
	})
	f := newFakeRegistry(t)
	ctx := context.Background()

	if err := Publish(ctx, PublishRequest{Path: source, Base: f.base()}); err != nil {
		t.Fatal(err)
	}
	repo, tag, _ := PullReference(m.Type, m.Name)
	if got := f.mediaTypeAt(t, repo, tag); got != ocispec.MediaTypeImageIndex {
		t.Errorf("tag resolves to %q, want an image index", got)
	}
}

// #ARCH:noarch has no platform spelling, and an index over a single child would
// imply the payload varies by architecture when it does not.
func TestPublishTagsANoarchArtifactDirectly(t *testing.T) {
	requireSquashfsTools(t)
	source, m := packImage(t, imageSpec{
		name: "grch38/genome/gencode49", typ: catalog.TypeData,
		arch: meta.ArchNone, recipe: "#!/bin/bash\n",
	})
	f := newFakeRegistry(t)
	ctx := context.Background()

	if err := Publish(ctx, PublishRequest{Path: source, Base: f.base()}); err != nil {
		t.Fatal(err)
	}
	repo, tag, _ := PullReference(m.Type, m.Name)
	if got := f.mediaTypeAt(t, repo, tag); got != ocispec.MediaTypeImageManifest {
		t.Errorf("tag resolves to %q, want a bare image manifest", got)
	}

	_, ann, err := ResolveArtifact(ctx, f.base(), repo, tag)
	if err != nil {
		t.Fatal(err)
	}
	if !IsNoarch(ann) {
		t.Errorf("annotations = %v, want the noarch mark", ann)
	}
}

// Pushing this architecture must not unpublish another one that was already
// there: the whole point of the index is that one tag serves both.
func TestPublishPreservesAnotherArchitecture(t *testing.T) {
	requireSquashfsTools(t)
	source, m := packImage(t, imageSpec{
		name: "grch38/genome/gencode49", typ: catalog.TypeData, recipe: "#!/bin/bash\n",
	})
	f := newFakeRegistry(t)
	ctx := context.Background()

	native := nativePlatform(t)
	foreign := ocispec.Platform{OS: "linux", Architecture: "s390x"}
	if native.Architecture == foreign.Architecture {
		t.Skip("running on the architecture used as the foreign one")
	}
	repo, tag, _ := PullReference(m.Type, m.Name)
	publishArtifact(t, f, repo, foreign, Annotations(m, ""), tag)

	if err := Publish(ctx, PublishRequest{Path: source, Base: f.base()}); err != nil {
		t.Fatalf("Publish: %v", err)
	}

	got := f.indexPlatformsAt(t, repo, tag)
	for _, want := range []string{platformKey(&native), platformKey(&foreign)} {
		if !slices.Contains(got, want) {
			t.Errorf("index carries %v, missing %q", got, want)
		}
	}
}

func TestPushRefusesAWritableOverlay(t *testing.T) {
	f := newFakeRegistry(t)
	err := Push(context.Background(), "/images/dev.img", f.base(), "dev", []string{"1.0"}, nil, false)
	if err == nil {
		t.Fatal("a writable overlay was pushed")
	}
	if !strings.Contains(err.Error(), "never distributed") {
		t.Errorf("the refusal does not say why: %v", err)
	}
}

func TestPushRefusesNoTags(t *testing.T) {
	f := newFakeRegistry(t)
	if err := Push(context.Background(), "/images/x.sqf", f.base(), "x", nil, nil, false); err == nil {
		t.Fatal("a push with no tags was accepted")
	}
}

// mediaTypeAt reports what a tag resolves to in the fake registry.
func (f *fakeRegistry) mediaTypeAt(t *testing.T, repo, tag string) string {
	t.Helper()
	f.mu.Lock()
	defer f.mu.Unlock()
	desc, ok := f.content[repo+"/"+tag]
	if !ok {
		t.Fatalf("nothing published at %s:%s", repo, tag)
	}
	return desc.MediaType
}

// indexPlatformsAt reports the platforms carried by the index at a tag.
func (f *fakeRegistry) indexPlatformsAt(t *testing.T, repo, tag string) []string {
	t.Helper()
	f.mu.Lock()
	defer f.mu.Unlock()
	desc, ok := f.content[repo+"/"+tag]
	if !ok {
		t.Fatalf("nothing published at %s:%s", repo, tag)
	}
	var idx ocispec.Index
	if err := json.Unmarshal(f.blobs[desc.Digest], &idx); err != nil {
		t.Fatalf("%s:%s is not an index: %v", repo, tag, err)
	}
	return indexPlatforms(idx)
}
