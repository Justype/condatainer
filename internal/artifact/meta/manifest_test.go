package meta

import (
	"encoding/json"
	"errors"
	"path/filepath"
	"reflect"
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
)

func TestValidateManifestAcceptsEachType(t *testing.T) {
	for _, typ := range []catalog.Type{catalog.TypeApp, catalog.TypeData, catalog.TypeOS, catalog.TypeBase} {
		m := validManifest()
		m.Type = typ
		if err := ValidateManifest(m); err != nil {
			t.Errorf("%s rejected: %v", typ, err)
		}
	}
}

func TestValidateManifestRejects(t *testing.T) {
	tests := []struct {
		name   string
		mutate func(*Manifest)
		want   error
	}{
		{"unsupported schema", func(m *Manifest) { m.SchemaVersion = SchemaVersion + 1 }, ErrUnsupportedSchema},
		{"zero schema", func(m *Manifest) { m.SchemaVersion = 0 }, ErrUnsupportedSchema},
		{"empty name", func(m *Manifest) { m.Name = "  " }, ErrInvalid},
		{"no architecture", func(m *Manifest) { m.Platform.Arch = "" }, ErrInvalid},
		{"unknown type", func(m *Manifest) { m.Type = "bundle" }, ErrInvalid},
		{"only identity key", func(m *Manifest) {
			m.Keys.Identity = KeyRef{Scheme: "conda-explicit-v1", SHA256: strings.Repeat("a", 64)}
		}, ErrInvalid},
		{"old keys without schemes", func(m *Manifest) {
			m.Keys.Identity = KeyRef{SHA256: strings.Repeat("a", 64)}
			m.Keys.Equiv = KeyRef{SHA256: strings.Repeat("b", 64)}
		}, ErrInvalid},
		{"invalid key digest", func(m *Manifest) {
			m.Keys.Identity = KeyRef{Scheme: "conda-explicit-v1", SHA256: "ABC"}
			m.Keys.Equiv = KeyRef{Scheme: "conda-environment-v1", SHA256: strings.Repeat("b", 64)}
		}, ErrInvalid},
		{"tools without Condatainer", func(m *Manifest) {
			m.Build.Tools = validBuildTools()
			m.Build.Tools.Condatainer = Tool{}
		}, ErrInvalid},
		{"tools without Apptainer implementation", func(m *Manifest) {
			m.Build.Tools = validBuildTools()
			m.Build.Tools.Apptainer.Name = ""
		}, ErrInvalid},
		{"Conda tools without Micromamba", func(m *Manifest) {
			m.Build.Tools = validBuildTools()
			m.Build.Tools.Micromamba = Tool{}
		}, ErrInvalid},
		{"script tools with Micromamba", func(m *Manifest) {
			m.BuildType = "script"
			m.Build.Tools = validBuildTools()
		}, ErrInvalid},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			m := validManifest()
			tt.mutate(&m)
			err := ValidateManifest(m)
			if err == nil {
				t.Fatalf("accepted %s", tt.name)
			}
			if !errors.Is(err, tt.want) {
				t.Fatalf("error = %v, want %v", err, tt.want)
			}
		})
	}
}

func TestManifestBuildToolsRoundTrip(t *testing.T) {
	m := validManifest()
	m.Build.Tools = validBuildTools()

	data, err := MarshalManifest(m)
	if err != nil {
		t.Fatal(err)
	}
	for _, want := range []string{`"tools"`, `"condatainer"`, `"apptainer"`, `"micromamba"`} {
		if !strings.Contains(string(data), want) {
			t.Errorf("manifest omits %s:\n%s", want, data)
		}
	}
	if strings.Contains(string(data), `"base"`) {
		t.Errorf("manifest unexpectedly records build base:\n%s", data)
	}

	var got Manifest
	if err := json.Unmarshal(data, &got); err != nil {
		t.Fatal(err)
	}
	if !reflect.DeepEqual(got.Build.Tools, m.Build.Tools) {
		t.Errorf("tools = %+v, want %+v", got.Build.Tools, m.Build.Tools)
	}
}

func validBuildTools() BuildTools {
	return BuildTools{
		Condatainer: Tool{Version: "1.4.2"},
		Apptainer:   Tool{Name: "apptainer", Version: "1.4.2"},
		Micromamba:  Tool{Version: "2.3.0"},
	}
}

// The manifest is off the mount path, so it carries no runtime block at all —
// nothing about a prefix, an environment, or what loading the image does.
func TestManifestCarriesNoRuntime(t *testing.T) {
	m := validManifest()
	m.Description = "SAMtools alignment toolkit"
	data, err := MarshalManifest(m)
	if err != nil {
		t.Fatalf("MarshalManifest: %v", err)
	}
	for _, field := range []string{"runtime", "prefix", "\"env\""} {
		if strings.Contains(string(data), field) {
			t.Errorf("manifest carries %s:\n%s", field, data)
		}
	}
	if !strings.HasSuffix(string(data), "\n") {
		t.Error("marshalled manifest has no trailing newline")
	}
}

func TestReadManifestRoundTrip(t *testing.T) {
	requireSquashfsTools(t)
	withTempCache(t)

	want := validManifest()
	want.Description = "SAMtools alignment toolkit"
	want.URL = "https://www.htslib.org/"

	got, err := ReadManifest(stagedImage(t, validRuntime(), want))
	if err != nil {
		t.Fatalf("ReadManifest: %v", err)
	}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("manifest = %+v, want %+v", got, want)
	}
}

// An image with runtime metadata but no manifest reports exactly that, distinct
// from a missing runtime document.
func TestReadManifestWithoutManifest(t *testing.T) {
	requireSquashfsTools(t)
	withTempCache(t)

	root := t.TempDir()
	if err := StageRuntime(filepath.Join(root, DirName), validRuntime()); err != nil {
		t.Fatalf("StageRuntime: %v", err)
	}
	if _, err := ReadManifest(packSqf(t, root)); !errors.Is(err, ErrNoManifest) {
		t.Fatalf("err = %v, want ErrNoManifest", err)
	}
}
