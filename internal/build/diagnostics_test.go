package build

import (
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
)

func TestNormalizedToolVersion(t *testing.T) {
	tests := []struct {
		name string
		raw  string
		want string
		ok   bool
	}{
		{name: "trimmed", raw: "  2.3.0\n", want: "2.3.0", ok: true},
		{name: "empty", raw: " \n", ok: false},
		{name: "multiline", raw: "2.3.0\nextra", ok: false},
		{name: "too long", raw: strings.Repeat("x", maxRecordedToolVersion+1), ok: false},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, ok := normalizedToolVersion(tt.raw)
			if got != tt.want || ok != tt.ok {
				t.Errorf("normalizedToolVersion() = %q, %v; want %q, %v", got, ok, tt.want, tt.ok)
			}
		})
	}
}

func TestBuildObjectManifestAddsCapturedTools(t *testing.T) {
	b := &BuildObject{
		spec: Spec{
			Image:  ImageSpec{Name: "samtools/1.23.1", Type: catalog.TypeApp},
			Source: SourceSpec{Script: &ScriptSource{}},
		},
		buildTools: meta.BuildTools{
			Condatainer: meta.Tool{Version: "1.4.2"},
			Apptainer:   meta.Tool{Name: "apptainer", Version: "1.4.2"},
		},
	}

	if !b.spec.Manifest().Build.Tools.Empty() {
		t.Fatal("Spec manifest contains execution diagnostics")
	}
	if got := b.Manifest().Build.Tools; got != b.buildTools {
		t.Errorf("manifest tools = %+v, want %+v", got, b.buildTools)
	}
}
