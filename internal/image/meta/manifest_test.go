package meta

import (
	"errors"
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
)

// valid returns a manifest that passes Validate, for a test to then break in
// exactly one way.
func valid() Manifest {
	return Manifest{
		SchemaVersion: SchemaVersion,
		Name:          "samtools/1.23.1",
		Type:          catalog.TypeApp,
		BuildType:     "conda",
		Runtime: Runtime{
			Prefix: "/cnt/samtools/1.23.1",
			Env:    []EnvVar{{Key: "SAMTOOLS_HOME", Value: "{prefix}"}},
		},
	}
}

func TestValidateAcceptsEachType(t *testing.T) {
	for _, typ := range []catalog.Type{catalog.TypeApp, catalog.TypeData} {
		m := valid()
		m.Type = typ
		if err := Validate(m); err != nil {
			t.Errorf("%s with a prefix rejected: %v", typ, err)
		}
	}
	// base and os apply at the container root, so they carry no prefix.
	for _, typ := range []catalog.Type{catalog.TypeBase, catalog.TypeOS} {
		m := valid()
		m.Type = typ
		m.Runtime.Prefix = ""
		if err := Validate(m); err != nil {
			t.Errorf("%s without a prefix rejected: %v", typ, err)
		}
	}
}

func TestValidateRejects(t *testing.T) {
	tests := []struct {
		name   string
		mutate func(*Manifest)
		want   error
	}{
		{"unsupported schema", func(m *Manifest) { m.SchemaVersion = SchemaVersion + 1 }, ErrUnsupportedSchema},
		{"zero schema", func(m *Manifest) { m.SchemaVersion = 0 }, ErrUnsupportedSchema},
		{"empty name", func(m *Manifest) { m.Name = "  " }, ErrInvalid},
		{"app without prefix", func(m *Manifest) { m.Runtime.Prefix = "" }, ErrInvalid},
		{"relative prefix", func(m *Manifest) { m.Runtime.Prefix = "cnt/samtools" }, ErrInvalid},
		{"empty env key", func(m *Manifest) { m.Runtime.Env[0].Key = "" }, ErrInvalid},
		{"env key with a dash", func(m *Manifest) { m.Runtime.Env[0].Key = "NOT-VALID" }, ErrInvalid},
		{"env key starting with a digit", func(m *Manifest) { m.Runtime.Env[0].Key = "1BAD" }, ErrInvalid},
		{"duplicate env key", func(m *Manifest) {
			m.Runtime.Env = append(m.Runtime.Env, EnvVar{Key: "SAMTOOLS_HOME", Value: "/elsewhere"})
		}, ErrInvalid},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			m := valid()
			tt.mutate(&m)
			err := Validate(m)
			if err == nil {
				t.Fatalf("accepted %s", tt.name)
			}
			if !errors.Is(err, tt.want) {
				t.Fatalf("error = %v, want %v", err, tt.want)
			}
		})
	}
}

// An absent or unrecognized type means app, which is what the old "Bundle
// Overlay" classification became.
func TestNormalizeDefaultsToApp(t *testing.T) {
	for _, raw := range []catalog.Type{"", "bundle", "module", "APP"} {
		m := valid()
		m.Type = raw
		m.Normalize()
		if m.Type != catalog.TypeApp {
			t.Errorf("type %q normalized to %q, want app", raw, m.Type)
		}
	}
	// A recognized type is left alone.
	for _, raw := range []catalog.Type{catalog.TypeBase, catalog.TypeOS, catalog.TypeApp, catalog.TypeData} {
		m := valid()
		m.Type = raw
		m.Normalize()
		if m.Type != raw {
			t.Errorf("type %q normalized to %q", raw, m.Type)
		}
	}
}

func TestPrefix(t *testing.T) {
	tests := []struct {
		name string
		typ  catalog.Type
		want string
	}{
		{"samtools/1.23.1", catalog.TypeApp, "/cnt/samtools/1.23.1"},
		{"grch38/star/2.7.11b/gencode47", catalog.TypeData, "/cnt/grch38/star/2.7.11b/gencode47"},
		{"ubuntu24/base", catalog.TypeBase, ""},
		{"ubuntu24/r-essential", catalog.TypeOS, ""},
	}
	for _, tt := range tests {
		if got := Prefix(tt.name, tt.typ); got != tt.want {
			t.Errorf("Prefix(%q, %s) = %q, want %q", tt.name, tt.typ, got, tt.want)
		}
	}
}

func TestEnvVarResolved(t *testing.T) {
	v := EnvVar{Key: "PATH_EXTRA", Value: "{prefix}/bin:{prefix}/libexec"}
	if got, want := v.Resolved("/cnt/x/1"), "/cnt/x/1/bin:/cnt/x/1/libexec"; got != want {
		t.Errorf("Resolved() = %q, want %q", got, want)
	}
	// A value with no token survives untouched.
	plain := EnvVar{Key: "K", Value: "/absolute"}
	if got := plain.Resolved("/cnt/x/1"); got != "/absolute" {
		t.Errorf("Resolved() = %q, want /absolute", got)
	}
}

// {prefix} has to reach the manifest intact: the install prefix is not known
// when the image is built, only when it is loaded.
func TestMarshalKeepsPrefixToken(t *testing.T) {
	data, err := Marshal(valid())
	if err != nil {
		t.Fatalf("Marshal: %v", err)
	}
	if !strings.Contains(string(data), `"value": "{prefix}"`) {
		t.Errorf("{prefix} did not survive marshalling:\n%s", data)
	}
	if !strings.HasSuffix(string(data), "\n") {
		t.Error("marshalled manifest has no trailing newline")
	}
}

func TestMarshalIsStable(t *testing.T) {
	first, err := Marshal(valid())
	if err != nil {
		t.Fatalf("Marshal: %v", err)
	}
	second, err := Marshal(valid())
	if err != nil {
		t.Fatalf("Marshal: %v", err)
	}
	if string(first) != string(second) {
		t.Error("marshalling the same manifest twice produced different bytes")
	}
}
