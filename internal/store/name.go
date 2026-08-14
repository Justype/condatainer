// Package store resolves immutable, identity-addressed overlays alongside the
// traditional flat image namespace.
package store

import (
	"errors"
	"fmt"
	"path/filepath"
	"regexp"
	"strings"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/image"
)

const (
	// DirName is the immutable-overflow directory below an image root.
	DirName = "store"
	// DefaultPrefixChars is the initial identity prefix length used in names.
	DefaultPrefixChars = 12
)

var (
	ErrInvalidIdentity = errors.New("invalid artifact identity")
	hexPattern         = regexp.MustCompile(`^[0-9a-f]+$`)
)

// IdentityQuery is either a complete identity or a user-supplied digest
// prefix. An empty Scheme means that all verified schemes are candidates.
type IdentityQuery struct {
	Scheme string
	SHA256 string
}

// Complete reports whether q is a complete scheme-backed identity.
func (q IdentityQuery) Complete() bool { return q.Scheme != "" && len(q.SHA256) == 64 }

// Matches reports whether ref satisfies the scheme and digest prefix in q.
func (q IdentityQuery) Matches(ref meta.KeyRef) bool {
	return (q.Scheme == "" || q.Scheme == ref.Scheme) && strings.HasPrefix(ref.SHA256, q.SHA256)
}

// ParseKeyRef parses the canonical complete-key spelling
// "scheme@sha256:<64 lowercase hex characters>".
func ParseKeyRef(raw string) (meta.KeyRef, error) {
	q, err := ParseIdentityQuery(raw)
	if err != nil {
		return meta.KeyRef{}, err
	}
	if !q.Complete() {
		return meta.KeyRef{}, fmt.Errorf("%w: complete keys require a scheme and 64 hexadecimal characters", ErrInvalidIdentity)
	}
	return meta.KeyRef{Scheme: q.Scheme, SHA256: q.SHA256}, nil
}

// ParseIdentityQuery accepts a canonical key, sha256 digest, or bare digest
// prefix. Scheme-free forms are conveniences and must be resolved against
// verified candidates before use.
func ParseIdentityQuery(raw string) (IdentityQuery, error) {
	raw = strings.TrimSpace(raw)
	if raw == "" {
		return IdentityQuery{}, fmt.Errorf("%w: empty value", ErrInvalidIdentity)
	}
	var q IdentityQuery
	if before, after, ok := strings.Cut(raw, "@"); ok {
		if before == "" || strings.Contains(after, "@") {
			return IdentityQuery{}, fmt.Errorf("%w: %q", ErrInvalidIdentity, raw)
		}
		q.Scheme, raw = before, after
	}
	if strings.HasPrefix(raw, "sha256:") {
		raw = strings.TrimPrefix(raw, "sha256:")
	} else if strings.Contains(raw, ":") {
		return IdentityQuery{}, fmt.Errorf("%w: only sha256 digests are supported", ErrInvalidIdentity)
	}
	raw = strings.ToLower(raw)
	if len(raw) == 0 || len(raw) > 64 || !hexPattern.MatchString(raw) {
		return IdentityQuery{}, fmt.Errorf("%w: digest must contain 1 to 64 hexadecimal characters", ErrInvalidIdentity)
	}
	q.SHA256 = raw
	return q, nil
}

// FormatKeyRef renders the complete, scheme-backed key used by the CLI.
func FormatKeyRef(ref meta.KeyRef) string {
	if ref.Scheme == "" || ref.SHA256 == "" {
		return ""
	}
	return ref.Scheme + "@sha256:" + ref.SHA256
}

// Filename returns a store filename using prefixChars identity characters.
func Filename(name string, identity meta.KeyRef, prefixChars int) (string, error) {
	if _, err := ParseKeyRef(FormatKeyRef(identity)); err != nil {
		return "", err
	}
	if prefixChars < DefaultPrefixChars || prefixChars > len(identity.SHA256) {
		return "", fmt.Errorf("identity prefix length must be between %d and %d", DefaultPrefixChars, len(identity.SHA256))
	}
	name = catalog.Normalize(name)
	if name == "" || strings.ContainsAny(name, `\\`) {
		return "", fmt.Errorf("invalid artifact name %q", name)
	}
	return image.EncodeArtifactName(name) + "@" + identity.SHA256[:prefixChars] + ".sqf", nil
}

// ParsedFilename is the untrusted addressing material in a store filename.
// Its full identity and scheme must still be regenerated from the image.
type ParsedFilename struct {
	Name   string
	Prefix string
}

// ParseFilename parses one store filename without trusting it as identity.
func ParseFilename(filename string) (ParsedFilename, error) {
	if filepath.Base(filename) != filename || filepath.Ext(filename) != ".sqf" {
		return ParsedFilename{}, fmt.Errorf("invalid store filename %q", filename)
	}
	stem := strings.TrimSuffix(filename, ".sqf")
	encoded, prefix, ok := strings.Cut(stem, "@")
	if !ok || encoded == "" || strings.Contains(prefix, "@") || len(prefix) < DefaultPrefixChars || len(prefix) > 64 || !hexPattern.MatchString(prefix) {
		return ParsedFilename{}, fmt.Errorf("invalid store filename %q", filename)
	}
	name := image.DecodeArtifactName(encoded)
	if image.EncodeArtifactName(name) != encoded {
		return ParsedFilename{}, fmt.Errorf("invalid encoded artifact name in %q", filename)
	}
	return ParsedFilename{Name: name, Prefix: prefix}, nil
}
