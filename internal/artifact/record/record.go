// Package record is the grammar of the text records CondaTainer embeds beside an
// image's payload: identity.record and equiv.record. It builds them, parses them
// back, and hashes them, and it knows nothing about images, recipes, or what
// belongs in a key — that is the key package's job.
//
// A record is UTF-8 text: a format line, then key=value lines, LF-terminated, no
// trailing whitespace, no blank lines, no comments. Text rather than JSON,
// deliberately: no canonicalization question, diffable in a terminal, hashable by
// hand with sha256sum.
package record

import (
	"crypto/sha256"
	"encoding/hex"
	"errors"
	"fmt"
	"sort"
	"strings"

	"github.com/Justype/condatainer/catalog"
)

// Kind is which of the two keys a record computes.
type Kind string

const (
	// KindIdentity answers "which recorded build is this?" and moves when
	// anything about the recorded build changes.
	KindIdentity Kind = "identity"
	// KindEquiv answers "may this substitute for what I asked for?" and moves
	// only when the artifact would behave differently for the asker.
	KindEquiv Kind = "equiv"
)

// Tag is the format line a record of this kind begins with. It versions the
// whole derivation: the field set, the ordering, and every rule that decides what
// reaches a line. Nothing may change under a tag without moving every key
// computed under it.
func (k Kind) Tag() string {
	switch k {
	case KindIdentity:
		return "cnt-identity-v1"
	case KindEquiv:
		return "cnt-equiv-v1"
	}
	return ""
}

// ErrInvalid reports a record that does not satisfy the grammar.
var ErrInvalid = errors.New("invalid record")

// DigestPrefix is how a record spells a hash inline.
const DigestPrefix = "sha256:"

// Env is one #ENV: contribution, with {prefix} left unsubstituted.
type Env struct {
	Key   string
	Value string
}

// Placeholder is one selected #PH: value.
type Placeholder struct {
	Name  string
	Value string
}

// Dep is one dependency line. What the fields mean is per-kind and decided
// elsewhere; here they are opaque space-free tokens, which is what lets identity
// and equivalence use the same line shape for different content.
type Dep struct {
	Type   catalog.Type
	Fields []string
}

// text renders the dependency as it appears after "dep=".
func (d Dep) text() string {
	return string(d.Type) + " " + strings.Join(d.Fields, " ")
}

// Record is one identity or equivalence record. Fields are written in this
// order — type, env, source, ph, dep — and repeated keys sort within their own
// group. No record ever carries a name or a prefix: a key describes what was
// built, and the name is what it was called.
type Record struct {
	Kind Kind
	Type catalog.Type
	Env  []Env
	// Recipe is the digest of the recipe body, empty for a build with no recipe.
	Recipe string
	// From is the digest of the upstream image a definition bootstrapped from,
	// empty when there is none. Identity records only: the reference as written
	// is already inside the recipe, so this adds which bytes it meant that day.
	From         string
	Placeholders []Placeholder
	Deps         []Dep
}

// Marshal renders the record in canonical form, sorting each group. The bytes it
// returns are the key's preimage, so two records with the same content marshal
// identically however their fields were ordered on the way in.
func Marshal(r Record) ([]byte, error) {
	if err := validate(r); err != nil {
		return nil, err
	}

	env := append([]Env(nil), r.Env...)
	sort.Slice(env, func(i, j int) bool { return env[i].Key < env[j].Key })
	ph := append([]Placeholder(nil), r.Placeholders...)
	sort.Slice(ph, func(i, j int) bool { return ph[i].Name < ph[j].Name })
	deps := append([]Dep(nil), r.Deps...)
	sort.SliceStable(deps, func(i, j int) bool { return deps[i].text() < deps[j].text() })

	var sb strings.Builder
	sb.WriteString(r.Kind.Tag())
	sb.WriteByte('\n')
	fmt.Fprintf(&sb, "type=%s\n", r.Type)
	for _, e := range env {
		fmt.Fprintf(&sb, "env=%s=%s\n", e.Key, e.Value)
	}
	if r.Recipe != "" {
		fmt.Fprintf(&sb, "recipe=%s\n", r.Recipe)
	}
	if r.From != "" {
		fmt.Fprintf(&sb, "from=%s\n", r.From)
	}
	for _, p := range ph {
		fmt.Fprintf(&sb, "ph=%s=%s\n", p.Name, p.Value)
	}
	for _, d := range deps {
		fmt.Fprintf(&sb, "dep=%s\n", d.text())
	}
	return []byte(sb.String()), nil
}

// Key renders the record and returns the digest of those bytes — the value the
// image is addressed and compared by.
func Key(r Record) (string, error) {
	data, err := Marshal(r)
	if err != nil {
		return "", err
	}
	return Digest(data), nil
}

// Digest returns the SHA-256 of data in the form a record writes it inline.
func Digest(data []byte) string { return DigestPrefix + Sum(data) }

// Sum returns the bare lowercase hex SHA-256, which is what sha256sum prints and
// what manifest.keys stores.
func Sum(data []byte) string {
	sum := sha256.Sum256(data)
	return hex.EncodeToString(sum[:])
}

// ValidDigest reports whether s is a well-formed inline digest, sha256: followed
// by 64 lowercase hex characters.
func ValidDigest(s string) bool {
	rest, ok := strings.CutPrefix(s, DigestPrefix)
	if !ok || len(rest) != 2*sha256.Size {
		return false
	}
	for _, c := range rest {
		if (c < '0' || c > '9') && (c < 'a' || c > 'f') {
			return false
		}
	}
	return true
}

// validate reports whether r can be rendered: a known kind, a type that can hold
// records, and field text that survives the round trip.
func validate(r Record) error {
	if r.Kind.Tag() == "" {
		return fmt.Errorf("%w: unknown kind %q", ErrInvalid, r.Kind)
	}
	if err := validSubjectType(r.Type); err != nil {
		return err
	}
	if r.From != "" && r.Kind != KindIdentity {
		return fmt.Errorf("%w: from belongs to identity only", ErrInvalid)
	}

	seenEnv := make(map[string]bool, len(r.Env))
	for _, e := range r.Env {
		if !validToken(e.Key) || strings.Contains(e.Key, "=") {
			return fmt.Errorf("%w: %q is not a usable env key", ErrInvalid, e.Key)
		}
		if !validValue(e.Value) {
			return fmt.Errorf("%w: env %s has an unusable value %q", ErrInvalid, e.Key, e.Value)
		}
		if seenEnv[e.Key] {
			return fmt.Errorf("%w: env %s appears more than once", ErrInvalid, e.Key)
		}
		seenEnv[e.Key] = true
	}

	if r.Recipe != "" && !ValidDigest(r.Recipe) {
		return fmt.Errorf("%w: recipe %q is not a sha256 digest", ErrInvalid, r.Recipe)
	}

	// A digest when the upstream resolved, a bare marker when it did not — the
	// same shape dependency fields use.
	if r.From != "" {
		if !validToken(r.From) {
			return fmt.Errorf("%w: %q is not a usable from value", ErrInvalid, r.From)
		}
		if strings.HasPrefix(r.From, DigestPrefix) && !ValidDigest(r.From) {
			return fmt.Errorf("%w: from %q is not a sha256 digest", ErrInvalid, r.From)
		}
	}

	seenPH := make(map[string]bool, len(r.Placeholders))
	for _, p := range r.Placeholders {
		if !validToken(p.Name) || strings.Contains(p.Name, "=") {
			return fmt.Errorf("%w: %q is not a usable placeholder name", ErrInvalid, p.Name)
		}
		if !validValue(p.Value) {
			return fmt.Errorf("%w: placeholder %s has an unusable value %q", ErrInvalid, p.Name, p.Value)
		}
		if seenPH[p.Name] {
			return fmt.Errorf("%w: placeholder %s appears more than once", ErrInvalid, p.Name)
		}
		seenPH[p.Name] = true
	}

	// Repeated dependency lines are kept, not collapsed: a recipe that mounted
	// two inputs which happen to be equivalent says so twice.
	for _, d := range r.Deps {
		if err := validDepType(d.Type); err != nil {
			return fmt.Errorf("%w (dependency)", err)
		}
		if len(d.Fields) == 0 {
			return fmt.Errorf("%w: dependency %s carries no fields", ErrInvalid, d.Type)
		}
		for _, f := range d.Fields {
			if !validToken(f) {
				return fmt.Errorf("%w: %q is not a usable dependency field", ErrInvalid, f)
			}
			if strings.HasPrefix(f, DigestPrefix) && !ValidDigest(f) {
				return fmt.Errorf("%w: %q is not a sha256 digest", ErrInvalid, f)
			}
		}
	}
	return nil
}

// validSubjectType reports whether an artifact type can be what a record
// describes. Every type can, a base included.
func validSubjectType(t catalog.Type) error {
	switch t {
	case catalog.TypeApp, catalog.TypeData, catalog.TypeOS, catalog.TypeBase:
		return nil
	default:
		return fmt.Errorf("%w: unknown type %q", ErrInvalid, t)
	}
}

// validDepType reports whether an artifact type can appear on a dep line. A base
// cannot: it is the environment a build runs inside rather than an input it
// consumed, so it never enters another artifact's history.
func validDepType(t catalog.Type) error {
	switch t {
	case catalog.TypeApp, catalog.TypeData, catalog.TypeOS:
		return nil
	case catalog.TypeBase:
		return fmt.Errorf("%w: a base is never a dependency", ErrInvalid)
	default:
		return fmt.Errorf("%w: unknown type %q", ErrInvalid, t)
	}
}

// validToken reports whether s can stand as a space-free field.
func validToken(s string) bool {
	if s == "" {
		return false
	}
	return !strings.ContainsFunc(s, func(r rune) bool { return r <= ' ' || r == 0x7f })
}

// validValue reports whether s can stand as the tail of a key=value line: it may
// contain spaces and '=', but not a line break, a control character, or edge
// whitespace that a reader would have to guess about.
func validValue(s string) bool {
	if s != strings.TrimSpace(s) {
		return false
	}
	return !strings.ContainsFunc(s, func(r rune) bool { return r < ' ' || r == 0x7f })
}
