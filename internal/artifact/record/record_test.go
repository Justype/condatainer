package record

import (
	"crypto/sha256"
	"encoding/hex"
	"errors"
	"reflect"
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
)

// dataIdentity is the worked example from the design: a STAR index built from a
// template recipe against three dependencies, one of which predates records.
func dataIdentity() Record {
	return Record{
		Kind:   KindIdentity,
		Type:   catalog.TypeData,
		Env:    []Env{{Key: "STAR_INDEX_DIR", Value: "{prefix}"}},
		Recipe: digestOf("recipe body"),
		Placeholders: []Placeholder{
			{Name: "gencode_version", Value: "49"},
			{Name: "read_length", Value: "101"},
			{Name: "star_version", Value: "2.7.11b"},
		},
		Deps: []Dep{
			{Type: catalog.TypeApp, Fields: []string{"samtools/1.23.1", "unrecorded"}},
			{Type: catalog.TypeApp, Fields: []string{"star/2.7.11b", digestOf("star")}},
			{Type: catalog.TypeData, Fields: []string{"grch38/gtf-gencode/49", digestOf("gtf")}},
		},
	}
}

func digestOf(s string) string {
	sum := sha256.Sum256([]byte(s))
	return DigestPrefix + hex.EncodeToString(sum[:])
}

func mustMarshal(t *testing.T, r Record) []byte {
	t.Helper()
	data, err := Marshal(r)
	if err != nil {
		t.Fatalf("Marshal: %v", err)
	}
	return data
}

// The rendered shape is the format, so it is asserted literally rather than
// through a round trip: a reader with sha256sum has to be able to reproduce a key
// by hand.
func TestMarshalShape(t *testing.T) {
	got := string(mustMarshal(t, dataIdentity()))
	want := strings.Join([]string{
		"cnt-identity-v1",
		"type=data",
		"env=STAR_INDEX_DIR={prefix}",
		"recipe=" + digestOf("recipe body"),
		"ph=gencode_version=49",
		"ph=read_length=101",
		"ph=star_version=2.7.11b",
		"dep=app samtools/1.23.1 unrecorded",
		"dep=app star/2.7.11b " + digestOf("star"),
		"dep=data grch38/gtf-gencode/49 " + digestOf("gtf"),
		"",
	}, "\n")
	if got != want {
		t.Errorf("marshalled record:\n%s\nwant:\n%s", got, want)
	}
}

// An equiv record for the same artifact drops the unnamed app entirely, carries
// no name for its data dependency, and no digest for its app one.
func TestMarshalEquivShape(t *testing.T) {
	r := Record{
		Kind:   KindEquiv,
		Type:   catalog.TypeData,
		Env:    []Env{{Key: "STAR_INDEX_DIR", Value: "{prefix}"}},
		Recipe: digestOf("recipe body"),
		Placeholders: []Placeholder{
			{Name: "star_version", Value: "2.7.11b"},
		},
		Deps: []Dep{
			{Type: catalog.TypeData, Fields: []string{digestOf("gtf equiv")}},
			{Type: catalog.TypeApp, Fields: []string{"star/2.7.11b"}},
		},
	}
	got := string(mustMarshal(t, r))
	want := strings.Join([]string{
		"cnt-equiv-v1",
		"type=data",
		"env=STAR_INDEX_DIR={prefix}",
		"recipe=" + digestOf("recipe body"),
		"ph=star_version=2.7.11b",
		"dep=app star/2.7.11b",
		"dep=data " + digestOf("gtf equiv"),
		"",
	}, "\n")
	if got != want {
		t.Errorf("marshalled record:\n%s\nwant:\n%s", got, want)
	}
}

// Marshal sorts, so a caller may build a record in whatever order its inputs
// arrive and still get the one preimage the key is computed from.
func TestMarshalSortsIndependentOfInputOrder(t *testing.T) {
	shuffled := dataIdentity()
	shuffled.Env = append([]Env{{Key: "AAA_FIRST", Value: "x"}}, shuffled.Env...)
	slices := shuffled.Placeholders
	shuffled.Placeholders = []Placeholder{slices[2], slices[0], slices[1]}
	shuffled.Deps = []Dep{shuffled.Deps[2], shuffled.Deps[1], shuffled.Deps[0]}

	ordered := dataIdentity()
	ordered.Env = append(ordered.Env, Env{Key: "AAA_FIRST", Value: "x"})

	if a, b := string(mustMarshal(t, shuffled)), string(mustMarshal(t, ordered)); a != b {
		t.Errorf("input order changed the preimage:\n%s\nvs\n%s", a, b)
	}
}

// Two distinct dependencies that happen to be equivalent produce two identical
// lines: the recipe mounted two inputs, and the record says so.
func TestMarshalKeepsRepeatedDependencies(t *testing.T) {
	r := Record{
		Kind: KindEquiv,
		Type: catalog.TypeData,
		Deps: []Dep{
			{Type: catalog.TypeData, Fields: []string{digestOf("same")}},
			{Type: catalog.TypeData, Fields: []string{digestOf("same")}},
		},
	}
	data := mustMarshal(t, r)
	if n := strings.Count(string(data), "dep=data "+digestOf("same")); n != 2 {
		t.Errorf("repeated dependency was collapsed (%d lines):\n%s", n, data)
	}
	back, err := Parse(data)
	if err != nil {
		t.Fatalf("Parse: %v", err)
	}
	if len(back.Deps) != 2 {
		t.Errorf("round trip dropped a repeated dependency: %+v", back.Deps)
	}
}

func TestRoundTrip(t *testing.T) {
	for _, r := range []Record{
		dataIdentity(),
		// A minimal record: an OS artifact with nothing but a recipe.
		{Kind: KindEquiv, Type: catalog.TypeOS, Recipe: digestOf("def")},
		// A script app: no dependencies, one env contribution.
		{Kind: KindIdentity, Type: catalog.TypeApp, Recipe: digestOf("sh"),
			Env: []Env{{Key: "TOOL_HOME", Value: "{prefix}"}}},
	} {
		data := mustMarshal(t, r)
		back, err := Parse(data)
		if err != nil {
			t.Fatalf("Parse(%s): %v", data, err)
		}
		if !reflect.DeepEqual(back, r) {
			t.Errorf("round trip changed the record:\ngot  %+v\nwant %+v", back, r)
		}
		again, err := Marshal(back)
		if err != nil {
			t.Fatalf("re-Marshal: %v", err)
		}
		if string(again) != string(data) {
			t.Errorf("re-marshalling changed the bytes:\n%s\nvs\n%s", again, data)
		}
	}
}

// A record that parses but does not re-marshal identically would hash to
// something other than the key it carries, so non-canonical input is rejected
// rather than repaired.
func TestParseRejects(t *testing.T) {
	canonical := string(mustMarshal(t, dataIdentity()))
	tests := []struct {
		name string
		text string
	}{
		{"empty", ""},
		{"no trailing newline", strings.TrimSuffix(canonical, "\n")},
		{"unknown format line", "cnt-identity-v2\ntype=data\n"},
		{"no format line", "type=data\n"},
		{"blank line", "cnt-equiv-v1\ntype=data\n\nrecipe=" + digestOf("x") + "\n"},
		{"trailing whitespace", "cnt-equiv-v1\ntype=data \n"},
		{"comment", "cnt-equiv-v1\n# a comment\ntype=data\n"},
		{"not key=value", "cnt-equiv-v1\ntype=data\nlonely\n"},
		{"unknown key", "cnt-equiv-v1\ntype=data\nsubdir=linux-64\n"},
		{"no type", "cnt-equiv-v1\nrecipe=" + digestOf("x") + "\n"},
		{"type twice", "cnt-equiv-v1\ntype=data\ntype=app\n"},
		{"recipe twice", "cnt-equiv-v1\ntype=data\nrecipe=" + digestOf("a") + "\nrecipe=" + digestOf("b") + "\n"},
		{"groups out of order", "cnt-equiv-v1\ntype=data\nph=a=1\nenv=K=v\n"},
		{"type after env", "cnt-equiv-v1\nenv=K=v\ntype=data\n"},
		{"env unsorted", "cnt-equiv-v1\ntype=data\nenv=ZZ=1\nenv=AA=2\n"},
		{"env duplicated", "cnt-equiv-v1\ntype=data\nenv=AA=1\nenv=AA=2\n"},
		{"ph unsorted", "cnt-equiv-v1\ntype=data\nph=zz=1\nph=aa=2\n"},
		{"deps unsorted", "cnt-equiv-v1\ntype=data\ndep=data " + digestOf("b") + "\ndep=app star/1\n"},
		{"unknown type", "cnt-equiv-v1\ntype=bundle\n"},
		{"unknown dependency type", "cnt-equiv-v1\ntype=data\ndep=bundle x\n"},
		{"a base is never a dependency", "cnt-equiv-v1\ntype=data\ndep=base ubuntu24/base " + digestOf("b") + "\n"},
		{"from in an equivalence record", "cnt-equiv-v1\ntype=os\nfrom=" + digestOf("u") + "\n"},
		{"from before recipe", "cnt-identity-v1\ntype=os\nfrom=" + digestOf("u") + "\nrecipe=" + digestOf("r") + "\n"},
		{"from twice", "cnt-identity-v1\ntype=os\nfrom=" + digestOf("u") + "\nfrom=" + digestOf("v") + "\n"},
		{"truncated from digest", "cnt-identity-v1\ntype=os\nfrom=sha256:019e\n"},
		{"dependency with no fields", "cnt-equiv-v1\ntype=data\ndep=app\n"},
		{"truncated digest", "cnt-equiv-v1\ntype=data\nrecipe=sha256:1a4f\n"},
		{"truncated digest in a dependency", "cnt-equiv-v1\ntype=data\ndep=data sha256:41ab\n"},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if _, err := Parse([]byte(tt.text)); err == nil {
				t.Fatalf("accepted %s:\n%s", tt.name, tt.text)
			} else if !errors.Is(err, ErrInvalid) {
				t.Errorf("err = %v, want ErrInvalid", err)
			}
		})
	}
}

func TestMarshalRejects(t *testing.T) {
	tests := []struct {
		name   string
		mutate func(*Record)
	}{
		{"unknown kind", func(r *Record) { r.Kind = "provenance" }},
		{"unknown type", func(r *Record) { r.Type = "bundle" }},
		{"a base as a dependency", func(r *Record) { r.Deps[0].Type = catalog.TypeBase }},
		{"from that is not a digest", func(r *Record) { r.From = "sha256:019e" }},
		{"from with a space", func(r *Record) { r.From = "two words" }},
		{"env key with a space", func(r *Record) { r.Env[0].Key = "TWO WORDS" }},
		{"env key with an equals", func(r *Record) { r.Env[0].Key = "A=B" }},
		{"env value with a newline", func(r *Record) { r.Env[0].Value = "one\ntwo" }},
		{"duplicate env key", func(r *Record) { r.Env = append(r.Env, r.Env[0]) }},
		{"duplicate placeholder", func(r *Record) { r.Placeholders = append(r.Placeholders, r.Placeholders[0]) }},
		{"placeholder value with a newline", func(r *Record) { r.Placeholders[0].Value = "a\nb" }},
		{"recipe that is not a digest", func(r *Record) { r.Recipe = "1a4f" }},
		{"dependency field with a space", func(r *Record) { r.Deps[0].Fields = []string{"two words"} }},
		{"dependency with no fields", func(r *Record) { r.Deps[0].Fields = nil }},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			r := dataIdentity()
			tt.mutate(&r)
			if _, err := Marshal(r); err == nil {
				t.Fatal("rendered a record that cannot round trip")
			} else if !errors.Is(err, ErrInvalid) {
				t.Errorf("err = %v, want ErrInvalid", err)
			}
		})
	}
}

// The key is the digest of exactly the bytes on disk, so sha256sum of the file
// reproduces it.
func TestKeyIsTheDigestOfTheBytes(t *testing.T) {
	r := dataIdentity()
	data := mustMarshal(t, r)
	key, err := Key(r)
	if err != nil {
		t.Fatalf("Key: %v", err)
	}
	if want := Digest(data); key != want {
		t.Errorf("Key = %s, want %s", key, want)
	}
	if !strings.HasSuffix(key, Sum(data)) {
		t.Errorf("Key %s does not end in the bare sum %s", key, Sum(data))
	}
	if !ValidDigest(key) {
		t.Errorf("Key %s is not a well-formed digest", key)
	}
}

// Content decides the key, not the order fields arrived in, and any change to
// content moves it.
func TestKeyMovesOnlyOnContent(t *testing.T) {
	base, err := Key(dataIdentity())
	if err != nil {
		t.Fatal(err)
	}

	same := dataIdentity()
	same.Deps = []Dep{same.Deps[1], same.Deps[0], same.Deps[2]}
	if got, _ := Key(same); got != base {
		t.Error("reordering the input moved the key")
	}

	moves := map[string]func(*Record){
		"a placeholder value": func(r *Record) { r.Placeholders[0].Value = "50" },
		"the recipe digest":   func(r *Record) { r.Recipe = digestOf("edited body") },
		"an env value":        func(r *Record) { r.Env[0].Value = "{prefix}/index" },
		"an env name":         func(r *Record) { r.Env[0].Key = "STAR_INDEX" },
		"a dependency":        func(r *Record) { r.Deps[1].Fields[1] = digestOf("rebuilt star") },
		"the kind":            func(r *Record) { r.Kind = KindEquiv },
	}
	for what, mutate := range moves {
		r := dataIdentity()
		mutate(&r)
		got, err := Key(r)
		if err != nil {
			t.Fatalf("%s: %v", what, err)
		}
		if got == base {
			t.Errorf("changing %s did not move the key", what)
		}
	}
}

func TestValidDigest(t *testing.T) {
	good := digestOf("x")
	if !ValidDigest(good) {
		t.Errorf("%s rejected", good)
	}
	for _, bad := range []string{
		"", "sha256:", "sha256:1a4f", strings.TrimPrefix(good, DigestPrefix),
		"md5:" + strings.TrimPrefix(good, DigestPrefix),
		DigestPrefix + strings.ToUpper(strings.TrimPrefix(good, DigestPrefix)),
		good + "0",
	} {
		if ValidDigest(bad) {
			t.Errorf("%q accepted as a digest", bad)
		}
	}
}
