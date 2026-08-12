// Package compare decides how one artifact relates to another: is this the
// build that was asked for, may it substitute for it, or is it something else.
//
// Everything downstream — store selection, `build --update`, dependency
// resolution, an explicit verify — calls this rather than restating the rules.
// It never runs on the execution path: mounting reads runtime.json and stops.
package compare

import (
	"fmt"
	"sort"
	"strings"

	"github.com/Justype/condatainer/internal/artifact/key"
	"github.com/Justype/condatainer/internal/artifact/meta"
)

// Verdict is how a candidate relates to what was asked for.
type Verdict string

const (
	// Exact: the same recorded build.
	Exact Verdict = "exact"
	// Equivalent: a different build that may substitute.
	Equivalent Verdict = "equivalent"
	// Different: comparable, and not the same.
	Different Verdict = "different"
	// Unverifiable: something is missing or malformed, so no claim is made.
	// It is not a synonym for different — it means the question was not answered.
	Unverifiable Verdict = "unverifiable"
)

// Diff is one input that moved between two artifacts.
type Diff struct {
	Field string // recipe, from, ph:<name>, env:<KEY>, dep:<name>, packages, channels, arch, type, name
	Want  string
	Got   string
}

func (d Diff) String() string { return fmt.Sprintf("%s: %s -> %s", d.Field, d.Want, d.Got) }

// Result is a verdict and why.
type Result struct {
	Verdict Verdict
	// Reason is set for a gate failure and for Unverifiable, where there is a
	// cause rather than a difference to list.
	Reason string
	// Diffs are the inputs that moved, for Different.
	Diffs []Diff
}

// Compare reports how a candidate relates to what was asked for, on this host.
func Compare(want, got Artifact) Result { return compareOn(want, got, meta.NativeArch()) }

// compareOn is Compare with the host architecture supplied, so a test can ask
// what a foreign host would decide.
func compareOn(want, got Artifact, hostArch string) Result {
	if r, failed := gates(want, got, hostArch); failed {
		return r
	}
	// Before the runtime gate, so an image whose records cannot be read is
	// reported for that rather than for the disagreement it causes downstream.
	if reason := got.usable(); reason != "" {
		return Result{Verdict: Unverifiable, Reason: reason}
	}
	if reason := want.usable(); reason != "" {
		return Result{Verdict: Unverifiable, Reason: "requested artifact: " + reason}
	}
	// A free integrity check: runtime.json was never part of the hashed preimage,
	// so an edited one cannot agree with the record beside it.
	if reason := got.runtimeAgrees(); reason != "" {
		return Result{Verdict: Unverifiable, Reason: reason}
	}
	if want.Format != got.Format {
		return Result{
			Verdict: Different,
			Diffs:   []Diff{{Field: "build_type", Want: want.Format, Got: got.Format}},
		}
	}
	if want.IdentityScheme != got.IdentityScheme || want.EquivScheme != got.EquivScheme {
		return Result{
			Verdict: Different,
			Reason:  "keys derived by different schemes are not directly comparable",
			Diffs: []Diff{
				{Field: "identity_scheme", Want: want.IdentityScheme, Got: got.IdentityScheme},
				{Field: "equiv_scheme", Want: want.EquivScheme, Got: got.EquivScheme},
			},
		}
	}

	switch {
	case want.Identity == got.Identity:
		return Result{Verdict: Exact}
	case want.Equiv == got.Equiv:
		return Result{Verdict: Equivalent, Diffs: diff(want, got)}
	default:
		return Result{Verdict: Different, Diffs: diff(want, got)}
	}
}

// gates are the structural checks that run before any key is looked at. A gate
// failure is a Different verdict with a specific diff, never a fallthrough to
// equivalence: two artifacts that fail a gate are not comparable, so a key match
// between them would be meaningless rather than reassuring.
func gates(want, got Artifact, hostArch string) (Result, bool) {
	// Architecture first, and against the host rather than against each other:
	// an artifact that cannot run here is not a candidate whatever it contains.
	if err := runsOn(got.Arch, hostArch); err != nil {
		return Result{
			Verdict: Different,
			Reason:  err.Error(),
			Diffs:   []Diff{{Field: "arch", Want: hostArch, Got: got.Arch}},
		}, true
	}
	if want.Type != got.Type {
		return Result{
			Verdict: Different,
			Reason:  "an app and a dataset are never interchangeable",
			Diffs:   []Diff{{Field: "type", Want: string(want.Type), Got: string(got.Type)}},
		}, true
	}
	// A substitution question is always about a requested name.
	if want.Name != got.Name {
		return Result{
			Verdict: Different,
			Reason:  "a substitution question is always about one name",
			Diffs:   []Diff{{Field: "name", Want: want.Name, Got: got.Name}},
		}, true
	}
	return Result{}, false
}

// runsOn reports whether an artifact of this architecture may run on this host.
// noarch runs anywhere; anything else runs only where it was built. Nothing
// infers portability — an artifact is portable only where its recipe said so.
func runsOn(arch, hostArch string) error {
	switch arch {
	case "":
		return fmt.Errorf("the image records no architecture")
	case meta.ArchNone, hostArch:
		return nil
	default:
		return fmt.Errorf("built for %s, this host is %s", arch, hostArch)
	}
}

// MountAllowed is the one comparison rule that runs on the execution path, and
// only because runtime.json is already in hand: an image built for another
// architecture mounts cleanly and fails somewhere further downstream, so the
// string comparison is worth doing at the point of mount.
func MountAllowed(rt meta.Runtime) error { return runsOn(rt.Platform.Arch, meta.NativeArch()) }

// diff reports which inputs moved. It reads whatever the two artifacts can
// support: canonical models name the exact field, and a Conda app — which has no
// canonical model — falls back to naming the export that differs.
func diff(want, got Artifact) []Diff {
	if want.identity == nil || got.identity == nil {
		return condaDiff(want, got)
	}
	out := modelDiff(*want.identity, *got.identity)

	// A dependency's own name is what a human needs, and equivalence models
	// deliberately drop it — so the manifest's list is what names the mover.
	if d := dependencyDiff(want.Dependencies, got.Dependencies); len(d) > 0 {
		out = append(out, d...)
	}
	return out
}

// modelDiff compares two identity models field by field.
func modelDiff(want, got key.Model) []Diff {
	var out []Diff
	if want.Recipe != got.Recipe {
		out = append(out, Diff{Field: "recipe", Want: short(want.Recipe), Got: short(got.Recipe)})
	}
	// The common case for an OS or a base: one tag, and an upstream that moved
	// underneath it.
	if want.From != got.From {
		out = append(out, Diff{Field: "from", Want: short(want.From), Got: short(got.From)})
	}
	for _, name := range union(placeholderNames(want), placeholderNames(got)) {
		a, b := placeholder(want, name), placeholder(got, name)
		if a != b {
			out = append(out, Diff{Field: "ph:" + name, Want: a, Got: b})
		}
	}
	for _, key := range union(envKeys(want), envKeys(got)) {
		a, b := env(want, key), env(got, key)
		if a != b {
			out = append(out, Diff{Field: "env:" + key, Want: a, Got: b})
		}
	}
	return out
}

// dependencyDiff names dependencies that moved, appeared, or vanished.
func dependencyDiff(want, got []meta.Dependency) []Diff {
	index := func(deps []meta.Dependency) map[string]meta.Dependency {
		out := make(map[string]meta.Dependency, len(deps))
		for _, d := range deps {
			out[d.Name] = d
		}
		return out
	}
	a, b := index(want), index(got)

	var names []string
	for name := range a {
		names = append(names, name)
	}
	for name := range b {
		if _, seen := a[name]; !seen {
			names = append(names, name)
		}
	}
	sort.Strings(names)

	var out []Diff
	for _, name := range names {
		x, inWant := a[name]
		y, inGot := b[name]
		switch {
		case !inWant:
			out = append(out, Diff{Field: "dep:" + name, Want: "absent", Got: describe(y)})
		case !inGot:
			out = append(out, Diff{Field: "dep:" + name, Want: describe(x), Got: "absent"})
		case x.Identity != y.Identity:
			out = append(out, Diff{Field: "dep:" + name, Want: describe(x), Got: describe(y)})
		}
	}
	return out
}

// condaDiff names which export moved, since a Conda app has no records to walk.
func condaDiff(want, got Artifact) []Diff {
	var out []Diff
	if want.Identity != got.Identity {
		out = append(out, Diff{Field: "packages", Want: short(want.Identity), Got: short(got.Identity)})
	}
	if want.Equiv != got.Equiv {
		out = append(out, Diff{Field: "environment", Want: short(want.Equiv), Got: short(got.Equiv)})
	}
	return out
}

func describe(d meta.Dependency) string {
	if d.Records == meta.Unrecorded {
		return meta.Unrecorded
	}
	return short(d.Identity.Digest())
}

// short trims a digest to the length a human reads, matching capsule addressing.
func short(digest string) string {
	trimmed := strings.TrimPrefix(digest, key.DigestPrefix)
	if len(trimmed) > 12 {
		return trimmed[:12]
	}
	return trimmed
}

func placeholderNames(r key.Model) []string {
	out := make([]string, 0, len(r.Placeholders))
	for _, p := range r.Placeholders {
		out = append(out, p.Name)
	}
	return out
}

func placeholder(r key.Model, name string) string {
	for _, p := range r.Placeholders {
		if p.Name == name {
			return p.Value
		}
	}
	return "absent"
}

func envKeys(r key.Model) []string {
	out := make([]string, 0, len(r.Env))
	for _, e := range r.Env {
		out = append(out, e.Key)
	}
	return out
}

func env(r key.Model, key string) string {
	for _, e := range r.Env {
		if e.Key == key {
			return e.Value
		}
	}
	return "absent"
}

// union merges two name lists, sorted and deduplicated.
func union(a, b []string) []string {
	seen := make(map[string]bool, len(a)+len(b))
	var out []string
	for _, list := range [][]string{a, b} {
		for _, name := range list {
			if !seen[name] {
				seen[name] = true
				out = append(out, name)
			}
		}
	}
	sort.Strings(out)
	return out
}
