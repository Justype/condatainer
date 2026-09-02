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

func TestAudienceAccepts(t *testing.T) {
	yes, no := true, false
	manifest := func(typ catalog.Type, buildType meta.BuildType) meta.Manifest {
		return meta.Manifest{Name: "x/1.0", Type: typ, BuildType: buildType}
	}
	declaring := func(typ catalog.Type, buildType meta.BuildType, answer *bool) meta.Manifest {
		m := manifest(typ, buildType)
		m.Redistribute = answer
		return m
	}
	tests := []struct {
		why      string
		audience Audience
		m        meta.Manifest
		wantErr  bool
	}{
		// Undeclared: the type default stands in for an unanswered question.
		{"a base is ours to publish", Public, manifest(catalog.TypeBase, meta.BuildTypeDef), false},
		{"apt packages from a public distribution", Public, manifest(catalog.TypeOS, meta.BuildTypeScript), false},
		{"public reference data", Public, manifest(catalog.TypeData, meta.BuildTypeScript), false},
		{"an app redistributes someone else's binaries", Public, manifest(catalog.TypeApp, meta.BuildTypeScript), true},

		// A Conda build embeds no recipe, so it can never carry the declaration
		// the app default would demand. Gating it would refuse it forever.
		{"a Conda app cannot declare, so it is not asked", Public, manifest(catalog.TypeApp, meta.BuildTypeConda), false},
		{"a Conda build of a permitted type", Public, manifest(catalog.TypeData, meta.BuildTypeConda), false},

		// Declared: the recipe answers and the type default does not apply.
		{"a declared app publishes", Public, declaring(catalog.TypeApp, meta.BuildTypeScript, &yes), false},
		{"a refusal beats a permissive type", Public, declaring(catalog.TypeData, meta.BuildTypeScript, &no), true},
		{"a refusal beats the Conda exemption", Public, declaring(catalog.TypeApp, meta.BuildTypeConda, &no), true},

		// A snapshot embeds no recipe, so it can carry no declaration — the same
		// mechanical reason a Conda build is exempt.
		{"a snapshot publishes with nothing to declare in", Public, manifest(catalog.TypeEnv, meta.BuildTypeSnapshot), false},

		{"a restricted endpoint takes an app", Restricted, manifest(catalog.TypeApp, meta.BuildTypeScript), false},
		{"a restricted endpoint takes a Conda build", Restricted, manifest(catalog.TypeApp, meta.BuildTypeConda), false},
		{"a restricted endpoint takes a refused artifact", Restricted, declaring(catalog.TypeApp, meta.BuildTypeScript, &no), false},
	}
	for _, tt := range tests {
		if err := tt.audience.Accepts(tt.m); (err != nil) != tt.wantErr {
			t.Errorf("%s: %s.Accepts(%s/%s) = %v, wantErr %v",
				tt.why, tt.audience, tt.m.Type, tt.m.BuildType, err, tt.wantErr)
		}
	}
}

// The refusal must name both ways out, or a user reads it as "CondaTainer will
// not publish apps" and stops.
func TestAudienceRefusalNamesBothRemedies(t *testing.T) {
	err := Public.Accepts(meta.Manifest{Name: "star/2.7.11b", Type: catalog.TypeApp, BuildType: meta.BuildTypeScript})
	if err == nil {
		t.Fatal("an undeclared app must be refused at a public endpoint")
	}
	for _, want := range []string{"#REDISTRIBUTE: yes", "restricted"} {
		if !strings.Contains(err.Error(), want) {
			t.Errorf("refusal does not mention %q: %v", want, err)
		}
	}
}

// The default is the restrictive one: an endpoint that never declared its
// audience must not be treated as private.
func TestPublishDefaultsToPublic(t *testing.T) {
	requireSquashfsTools(t)
	source, _ := packImage(t, imageSpec{name: "hello/1.0", typ: catalog.TypeApp, recipe: "#!/bin/bash\n"})
	f := newFakeRegistry(t)

	_, err := Publish(context.Background(), PublishRequest{Path: source, Base: f.base()})
	if err == nil {
		t.Fatal("an app was published to an endpoint with no declared audience")
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

	if _, err := Publish(ctx, PublishRequest{Path: source, Base: f.base()}); err != nil {
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

	_, err := Publish(context.Background(), PublishRequest{Path: source, Base: f.base()})
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

	_, err := Publish(context.Background(), PublishRequest{Path: source, Base: f.base()})
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

	if _, err := Publish(ctx, PublishRequest{Path: source, Base: f.base()}); err != nil {
		t.Fatalf("first publish: %v", err)
	}
	_, err := Publish(ctx, PublishRequest{Path: source, Base: f.base()})
	if err == nil {
		t.Fatal("a versioned tag was silently replaced")
	}
	if !strings.Contains(err.Error(), "--force") {
		t.Errorf("the refusal does not name the way through: %v", err)
	}
	if _, err := Publish(ctx, PublishRequest{Path: source, Base: f.base(), Force: true}); err != nil {
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

	if _, err := Publish(ctx, PublishRequest{Path: source, Base: f.base()}); err != nil {
		t.Fatalf("Publish: %v", err)
	}
	if _, err := Publish(ctx, PublishRequest{Path: source, Base: f.base()}); err != nil {
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

	if _, err := Publish(ctx, PublishRequest{Path: source, Base: f.base()}); err != nil {
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

	if _, err := Publish(ctx, PublishRequest{Path: source, Base: f.base()}); err != nil {
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

	if _, err := Publish(ctx, PublishRequest{Path: source, Base: f.base()}); err != nil {
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
	_, err := Push(context.Background(), "/images/dev.img", f.base(), "dev", []string{"1.0"}, nil, false)
	if err == nil {
		t.Fatal("a writable overlay was pushed")
	}
	if !strings.Contains(err.Error(), "never distributed") {
		t.Errorf("the refusal does not say why: %v", err)
	}
}

func TestPushRefusesNoTags(t *testing.T) {
	f := newFakeRegistry(t)
	if _, err := Push(context.Background(), "/images/x.sqf", f.base(), "x", nil, nil, false); err == nil {
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

// A project keeps artifacts from several collections in one repository, so that
// package's source is the project rather than whichever collection built each
// artifact. The publisher supplies it; without one the artifact's own answer
// stands.
func TestPublishSourceOverridesTheRecordedCollection(t *testing.T) {
	requireSquashfsTools(t)
	source, m := packImage(t, imageSpec{
		name: "star/2.7.11b", typ: catalog.TypeApp, recipe: "#!/bin/bash\nbuild\n",
	})
	ctx := context.Background()

	f := newFakeRegistry(t)
	project := "https://github.com/my-lab/rnaseq-2026"
	if _, err := Publish(ctx, PublishRequest{
		Path: source, Base: f.base(), Audience: Restricted,
		Placement: &Placement{Repo: "cnt", Tags: []string{"star--2.7.11b"}},
		Source:    project,
	}); err != nil {
		t.Fatalf("Publish: %v", err)
	}
	_, ann, err := ResolveArtifact(ctx, f.base(), "cnt", "star--2.7.11b")
	if err != nil {
		t.Fatalf("ResolveArtifact: %v", err)
	}
	if got := ann[AnnSource]; got != project {
		t.Errorf("%s = %q, want %q", AnnSource, got, project)
	}

	// No override: whatever the artifact recorded, which for this fixture is
	// nothing, so the annotation is absent rather than empty.
	g := newFakeRegistry(t)
	if _, err := Publish(ctx, PublishRequest{
		Path: source, Base: g.base(), Audience: Restricted,
	}); err != nil {
		t.Fatalf("Publish: %v", err)
	}
	_, ann, err = ResolveArtifact(ctx, g.base(), "star", "2.7.11b")
	if err != nil {
		t.Fatalf("ResolveArtifact: %v", err)
	}
	if got, ok := ann[AnnSource]; ok && got != m.Build.Source {
		t.Errorf("%s = %q, want the artifact's own %q", AnnSource, got, m.Build.Source)
	}
}

// The annotation exists so a registry can link a package to a repository, which
// it does by exact URL match. Anything that is not an http(s) URL links nothing,
// so it never reaches a manifest.
func TestPublishKeepsTheSourceAnnotationUsable(t *testing.T) {
	requireSquashfsTools(t)
	ctx := context.Background()

	// A caller that names one and gets it wrong is told.
	clean, _ := packImage(t, imageSpec{
		name: "star/2.7.11b", typ: catalog.TypeApp, recipe: "#!/bin/bash\nbuild\n",
	})
	f := newFakeRegistry(t)
	_, err := Publish(ctx, PublishRequest{
		Path: clean, Base: f.base(), Audience: Restricted, Source: "ftp://example.invalid/p",
	})
	if err == nil || !strings.Contains(err.Error(), "http(s)") {
		t.Fatalf("error = %v, want a refusal naming the scheme", err)
	}

	// An artifact whose own recorded source is unusable still publishes, because
	// a descriptor that predates the check is not a reason to block distribution.
	// The annotation is simply absent.
	dirty, _ := packImage(t, imageSpec{
		name: "star/2.7.11b", typ: catalog.TypeApp, recipe: "#!/bin/bash\nbuild\n",
		tamper: func(m *meta.Manifest) { m.Build.Source = "not a url" },
	})
	g := newFakeRegistry(t)
	if _, err := Publish(ctx, PublishRequest{
		Path: dirty, Base: g.base(), Audience: Restricted,
	}); err != nil {
		t.Fatalf("Publish: %v", err)
	}
	_, ann, err := ResolveArtifact(ctx, g.base(), "star", "2.7.11b")
	if err != nil {
		t.Fatalf("ResolveArtifact: %v", err)
	}
	if got, ok := ann[AnnSource]; ok {
		t.Errorf("%s = %q, want it dropped", AnnSource, got)
	}
}
