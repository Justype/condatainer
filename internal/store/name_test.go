package store

import (
	"errors"
	"strings"
	"testing"

	"github.com/Justype/condatainer/internal/artifact/meta"
)

func TestKeyRefAndIdentityQuery(t *testing.T) {
	hex := strings.Repeat("a", 64)
	ref, err := ParseKeyRef("script-identity-v1@sha256:" + hex)
	if err != nil {
		t.Fatal(err)
	}
	if ref.Scheme != "script-identity-v1" || ref.SHA256 != hex {
		t.Fatalf("unexpected key: %#v", ref)
	}
	if got := FormatKeyRef(ref); got != "script-identity-v1@sha256:"+hex {
		t.Fatalf("FormatKeyRef() = %q", got)
	}
	q, err := ParseIdentityQuery("A31F902C12AB")
	if err != nil {
		t.Fatal(err)
	}
	if q.Scheme != "" || q.SHA256 != "a31f902c12ab" {
		t.Fatalf("unexpected query: %#v", q)
	}
	if _, err := ParseKeyRef("sha256:" + hex); !errors.Is(err, ErrInvalidIdentity) {
		t.Fatalf("scheme-free complete key error = %v", err)
	}
	for _, invalid := range []string{"", "scheme@sha512:abc", "scheme@sha256:xyz", strings.Repeat("a", 65)} {
		if _, err := ParseIdentityQuery(invalid); !errors.Is(err, ErrInvalidIdentity) {
			t.Errorf("ParseIdentityQuery(%q) error = %v", invalid, err)
		}
	}
}

func TestStoreFilenameRoundTrip(t *testing.T) {
	identity := meta.KeyRef{Scheme: "script-identity-v1", SHA256: strings.Repeat("b", 64)}
	name, err := Filename("grch38/genome/gencode49", identity, 16)
	if err != nil {
		t.Fatal(err)
	}
	if name != "grch38--genome--gencode49@bbbbbbbbbbbbbbbb.sqf" {
		t.Fatalf("Filename() = %q", name)
	}
	parsed, err := ParseFilename(name)
	if err != nil {
		t.Fatal(err)
	}
	if parsed.Name != "grch38/genome/gencode49" || parsed.Prefix != strings.Repeat("b", 16) {
		t.Fatalf("unexpected parsed filename: %#v", parsed)
	}
	for _, invalid := range []string{"x@abc.sqf", "x@" + strings.Repeat("g", 12) + ".sqf", "x.sqf", "../x@" + strings.Repeat("a", 12) + ".sqf"} {
		if _, err := ParseFilename(invalid); err == nil {
			t.Errorf("ParseFilename(%q) unexpectedly succeeded", invalid)
		}
	}
}
