package cmd

import (
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
)

func TestPresentListOverlayTrustsRecordedType(t *testing.T) {
	tests := []struct {
		name     string
		encoded  string
		runtime  meta.Runtime
		recorded bool
		want     listPresentation
	}{
		{
			name: "deep explicitly app", encoded: "team--tool--1.0", recorded: true,
			runtime: meta.Runtime{Name: "team/tool/1.0", Type: catalog.TypeApp},
			want:    listPresentation{name: "team/tool", version: "1.0"},
		},
		{
			name: "shallow explicitly data", encoded: "reference--1", recorded: true,
			runtime: meta.Runtime{Name: "reference/1", Type: catalog.TypeData},
			want:    listPresentation{name: "reference/1", data: true},
		},
		{
			name: "os", encoded: "ubuntu24--tool--1", recorded: true,
			runtime: meta.Runtime{Name: "ubuntu24/tool/1", Type: catalog.TypeOS},
			want:    listPresentation{name: "ubuntu24/tool/1", version: "(system app)"},
		},
		{
			name: "legacy depth fallback", encoded: "reference--genome--1", recorded: false,
			want: listPresentation{name: "reference/genome/1", data: true},
		},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			if got := presentListOverlay(tc.encoded, tc.runtime, tc.recorded); got != tc.want {
				t.Fatalf("presentListOverlay() = %#v, want %#v", got, tc.want)
			}
		})
	}
}
