package publish

import (
	"context"
	"encoding/json"
	"net/http"
	"net/http/httptest"
	"path/filepath"
	"strings"
	"testing"

	"github.com/opencontainers/go-digest"
	ocispec "github.com/opencontainers/image-spec/specs-go/v1"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/capsule"
	"github.com/Justype/condatainer/internal/artifact/key"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/project/lock"
	"github.com/Justype/condatainer/internal/registry"
	"github.com/Justype/condatainer/internal/utils"
)

const collectionRepo = "https://github.com/lab/recipes"

// vendorData writes a data artifact into cnt-lock/provenance/ and returns its
// relative path with the manifest it recorded. Data rather than app, so the
// public-endpoint skip does not fire on every fixture.
func vendorData(t *testing.T, root, name, recipe string, deps ...meta.Dependency) (string, meta.Manifest) {
	t.Helper()
	manifest := meta.Manifest{
		SchemaVersion: meta.SchemaVersion,
		Name:          name,
		Type:          catalog.TypeData,
		BuildType:     "script",
		Platform:      meta.Platform{OS: "linux", Arch: "amd64"},
		Source:        meta.Source{Files: []string{meta.RecipeFileName}},
		Dependencies:  deps,
		Build:         meta.Build{Source: collectionRepo},
	}
	derived, err := key.Generate(manifest, key.Sources{meta.RecipeFileName: []byte(recipe)})
	if err != nil {
		t.Fatalf("generate keys: %v", err)
	}
	manifest.Keys = derived.Keys()

	body, err := json.MarshalIndent(manifest, "", "  ")
	if err != nil {
		t.Fatal(err)
	}
	relative, err := lock.StageEntry(root, capsule.EntryName(name, manifest.Keys.Identity.Digest()),
		map[string][]byte{meta.FileName: body, meta.RecipeFileName: []byte(recipe)})
	if err != nil {
		t.Fatalf("stage entry: %v", err)
	}
	return relative, manifest
}

// vendorApp writes an app artifact, which a public endpoint refuses unless its
// recipe declares itself redistributable.
func vendorApp(t *testing.T, root, name, recipe string) string {
	t.Helper()
	manifest := meta.Manifest{
		SchemaVersion: meta.SchemaVersion,
		Name:          name,
		Type:          catalog.TypeApp,
		BuildType:     "script",
		Platform:      meta.Platform{OS: "linux", Arch: "amd64"},
		Source:        meta.Source{Files: []string{meta.RecipeFileName}},
		Build:         meta.Build{Source: collectionRepo},
	}
	derived, err := key.Generate(manifest, key.Sources{meta.RecipeFileName: []byte(recipe)})
	if err != nil {
		t.Fatalf("generate keys: %v", err)
	}
	manifest.Keys = derived.Keys()
	body, err := json.MarshalIndent(manifest, "", "  ")
	if err != nil {
		t.Fatal(err)
	}
	relative, err := lock.StageEntry(root, capsule.EntryName(name, manifest.Keys.Identity.Digest()),
		map[string][]byte{meta.FileName: body, meta.RecipeFileName: []byte(recipe)})
	if err != nil {
		t.Fatalf("stage entry: %v", err)
	}
	return relative
}

func projectRoot(t *testing.T) string {
	t.Helper()
	root := t.TempDir()
	if err := utils.MkdirAllShared(filepath.Join(root, lock.DirName)); err != nil {
		t.Fatal(err)
	}
	return root
}

// serveManifest stands up a registry that publishes one plain manifest under
// every tag, carrying the identity annotations given. Plain rather than an
// index, so the fixture does not depend on the running architecture.
func serveManifest(t *testing.T, annotations map[string]string) string {
	t.Helper()
	manifest := ocispec.Manifest{
		MediaType:    ocispec.MediaTypeImageManifest,
		ArtifactType: registry.ArtifactTypeOverlay,
		Config:       ocispec.DescriptorEmptyJSON,
		Layers:       []ocispec.Descriptor{{MediaType: registry.MediaTypeOverlayBlob, Digest: digest.FromString("payload"), Size: 7}},
		Annotations:  annotations,
	}
	body, err := json.Marshal(manifest)
	if err != nil {
		t.Fatal(err)
	}
	dgst := digest.FromBytes(body)

	server := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		if !strings.Contains(r.URL.Path, "/manifests/") {
			w.WriteHeader(http.StatusNotFound)
			return
		}
		w.Header().Set("Content-Type", ocispec.MediaTypeImageManifest)
		w.Header().Set("Docker-Content-Digest", dgst.String())
		if r.Method == http.MethodHead {
			w.Header().Set("Content-Length", itoa(len(body)))
			return
		}
		w.Write(body) //nolint:errcheck
	}))
	t.Cleanup(server.Close)
	return strings.TrimPrefix(server.URL, "http://")
}

func itoa(n int) string {
	if n == 0 {
		return "0"
	}
	var b []byte
	for n > 0 {
		b = append([]byte{byte('0' + n%10)}, b...)
		n /= 10
	}
	return string(b)
}

func collection(t *testing.T, endpoint string) catalog.Catalog {
	t.Helper()
	return catalog.Catalog{{
		Name: "lab",
		Desc: catalog.Descriptor{
			Schema: 1,
			Source: collectionRepo,
			OCI:    catalog.OCI{Pull: []string{endpoint}, Audience: string(registry.Public)},
		},
	}}
}

// The identity the endpoint advertises is the whole test: a lock remote is an
// address for one exact build, so anything else must record nothing.
func TestUpstreamRecordsOnlyAnExactIdentityMatch(t *testing.T) {
	root := projectRoot(t)
	artifact, manifest := vendorData(t, root, "grch38/genome/gencode49", "echo index\n")

	exact := map[string]string{
		registry.AnnTitle:          manifest.Name,
		registry.AnnSchema:         "1",
		registry.AnnIdentityScheme: manifest.Keys.Identity.Scheme,
		registry.AnnIdentitySHA:    manifest.Keys.Identity.SHA256,
		registry.AnnEquivScheme:    manifest.Keys.Equiv.Scheme,
		registry.AnnEquivSHA:       manifest.Keys.Equiv.SHA256,
	}
	got := Upstream(context.Background(), root, []string{artifact}, collection(t, serveManifest(t, exact)))
	if len(got[artifact]) != 1 {
		t.Fatalf("an exact identity match recorded %#v", got)
	}
	if !strings.HasSuffix(got[artifact][0].Repository, "/grch38/genome") {
		t.Errorf("repository = %q, want the catalog scheme's repo path", got[artifact][0].Repository)
	}
	if !strings.HasPrefix(got[artifact][0].ManifestDigest, "sha256:") {
		t.Errorf("manifest digest = %q", got[artifact][0].ManifestDigest)
	}

	// Same equivalence, different identity: a different build wearing the right
	// label. Restore would fetch it and then reject it against the lock, so it
	// must never be recorded as an address for this one.
	substitute := map[string]string{}
	for k, v := range exact {
		substitute[k] = v
	}
	substitute[registry.AnnIdentitySHA] = strings.Repeat("c", 64)
	if got := Upstream(context.Background(), root, []string{artifact}, collection(t, serveManifest(t, substitute))); len(got) != 0 {
		t.Errorf("an equivalent-but-different build was recorded: %#v", got)
	}
}

// Locking has to work offline. Every one of these is silence, not an error.
func TestUpstreamIsBestEffort(t *testing.T) {
	root := projectRoot(t)
	artifact, _ := vendorData(t, root, "grch38/genome/gencode49", "echo index\n")
	ctx := context.Background()

	cases := map[string]catalog.Catalog{
		"no configured source":    nil,
		"an unreachable registry": collection(t, "127.0.0.1:1"),
		"a collection that declares no pull endpoint": {{
			Name: "lab",
			Desc: catalog.Descriptor{Schema: 1, Source: collectionRepo},
		}},
		"a collection that published nothing this artifact came from": {{
			Name: "other",
			Desc: catalog.Descriptor{Schema: 1, Source: "https://github.com/somebody/else",
				OCI: catalog.OCI{Pull: []string{"127.0.0.1:1"}}},
		}},
	}
	for why, cat := range cases {
		if got := Upstream(ctx, root, []string{artifact}, cat); len(got) != 0 {
			t.Errorf("%s recorded %#v, want nothing", why, got)
		}
	}
	if got := Upstream(ctx, root, []string{"provenance/does--not@000000000000"}, collection(t, "127.0.0.1:1")); len(got) != 0 {
		t.Errorf("an unreadable artifact recorded %#v", got)
	}
}

// Two collections claiming one repository: nothing can say which published the
// artifact, and picking either would be a guess written into a tracked file.
func TestUpstreamRefusesAnAmbiguousCollection(t *testing.T) {
	root := projectRoot(t)
	artifact, manifest := vendorData(t, root, "grch38/genome/gencode49", "echo index\n")
	endpoint := serveManifest(t, map[string]string{
		registry.AnnTitle:          manifest.Name,
		registry.AnnSchema:         "1",
		registry.AnnIdentityScheme: manifest.Keys.Identity.Scheme,
		registry.AnnIdentitySHA:    manifest.Keys.Identity.SHA256,
	})
	cat := append(collection(t, endpoint), collection(t, endpoint)...)
	if got := Upstream(context.Background(), root, []string{artifact}, cat); len(got) != 0 {
		t.Errorf("an ambiguous collection recorded %#v", got)
	}
}

// vendorEnv writes a frozen environment, whose entry is a manifest alone: a
// snapshot names no sources, and its keys are the tree hash of the payload the
// pack produced rather than anything regenerable from a checkout.
func vendorEnv(t *testing.T, root, identity string) string {
	t.Helper()
	ref := meta.KeyRef{Scheme: string(key.PayloadTreeV1), SHA256: identity}
	manifest := meta.Manifest{
		SchemaVersion: meta.SchemaVersion,
		Name:          meta.EnvName,
		Type:          catalog.TypeEnv,
		BuildType:     meta.BuildTypeSnapshot,
		Platform:      meta.Platform{OS: "linux", Arch: "amd64"},
		Keys:          meta.Keys{Identity: ref, Equiv: ref},
	}
	body, err := json.MarshalIndent(manifest, "", "  ")
	if err != nil {
		t.Fatal(err)
	}
	relative, err := lock.StageEntry(root, capsule.EntryName(meta.EnvName, ref.Digest()),
		map[string][]byte{meta.FileName: body})
	if err != nil {
		t.Fatalf("stage entry: %v", err)
	}
	return relative
}
