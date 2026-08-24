package compare

import (
	"strings"
	"testing"

	"github.com/Justype/condatainer/internal/artifact/meta"
)

// An Artifact digest is "sha256:<hex>" and a meta.KeyRef holds the bare hex.
// Crossing that by hand yields a key equal to nothing, which fails silently:
// a lookup simply never matches. It has been got wrong at four call sites, so
// the crossing is pinned here rather than at each of them.
func TestRefsRoundTripBackToTheDigest(t *testing.T) {
	const hex = "41ab1c2d3e4f5061728394a5b6c7d8e9f00112233445566778899aabbccddeeff"
	a := Artifact{
		Identity: "sha256:" + hex, IdentityScheme: "script-identity-v1",
		Equiv: "sha256:" + hex, EquivScheme: "script-equiv-v1",
	}

	if got := a.IdentityRef().Digest(); got != a.Identity {
		t.Errorf("IdentityRef().Digest() = %q, want %q", got, a.Identity)
	}
	if got := a.EquivRef().Digest(); got != a.Equiv {
		t.Errorf("EquivRef().Digest() = %q, want %q", got, a.Equiv)
	}
	if strings.HasPrefix(a.IdentityRef().SHA256, "sha256:") {
		t.Errorf("IdentityRef kept the prefix: %q", a.IdentityRef().SHA256)
	}
	if strings.HasPrefix(a.EquivRef().SHA256, "sha256:") {
		t.Errorf("EquivRef kept the prefix: %q", a.EquivRef().SHA256)
	}
}

// The scheme travels with the digest. A KeyRef missing one is Empty, so a ref
// that dropped it would read as "no key recorded" rather than as a mismatch.
func TestRefsCarryTheirSchemes(t *testing.T) {
	a := Artifact{
		Identity: "sha256:aa", IdentityScheme: "conda-explicit-v1",
		Equiv: "sha256:bb", EquivScheme: "conda-environment-v1",
	}
	if got := a.IdentityRef(); got.Scheme != "conda-explicit-v1" || got.SHA256 != "aa" {
		t.Errorf("IdentityRef() = %+v", got)
	}
	if got := a.EquivRef(); got.Scheme != "conda-environment-v1" || got.SHA256 != "bb" {
		t.Errorf("EquivRef() = %+v", got)
	}
}

// An unkeyed artifact must produce an empty ref, not a half-populated one that
// compares unequal to every real key while looking like a claim.
func TestRefsOfAnUnkeyedArtifactAreEmpty(t *testing.T) {
	var a Artifact
	if !a.IdentityRef().Empty() || !a.EquivRef().Empty() {
		t.Errorf("refs of an unkeyed artifact should be empty: %+v %+v", a.IdentityRef(), a.EquivRef())
	}
	if got := (meta.KeyRef{}).Digest(); got != "" {
		t.Errorf("an empty ref should render no digest, got %q", got)
	}
}
