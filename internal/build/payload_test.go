package build

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
)

// Only a recipe's build carries a payload key: a Conda environment and a
// definition are pinned by their own keys, and keying an OS rootfs is the
// slowest walk there is.
func TestOnlyAScriptBuildIsPayloadKeyed(t *testing.T) {
	for _, tc := range []struct {
		buildType BuildType
		want      bool
	}{
		{BuildTypeScript, true},
		{BuildTypeConda, false},
		{BuildTypeDef, false},
	} {
		b := newPackObject(t, catalog.TypeApp)
		b.buildType = tc.buildType
		if err := os.MkdirAll(b.ws.CntDir, 0o755); err != nil {
			t.Fatal(err)
		}
		if err := os.WriteFile(filepath.Join(b.ws.CntDir, "tool"), []byte("x"), 0o755); err != nil {
			t.Fatal(err)
		}
		if err := meta.StageManifest(b.ws.MetaDir, b.Manifest()); err != nil {
			t.Fatal(err)
		}

		if err := b.stagePayloadKey(t.Context(), b.ws.CntDir, b.ws.MetaDir); err != nil {
			t.Fatalf("%s: %v", tc.buildType, err)
		}

		data, err := os.ReadFile(filepath.Join(b.ws.MetaDir, meta.FileName))
		if err != nil {
			t.Fatal(err)
		}
		m, err := meta.DecodeManifest(data, "staged")
		if err != nil {
			t.Fatal(err)
		}
		if got := !m.Keys.Payload.Empty(); got != tc.want {
			t.Errorf("%s: payload key recorded = %v, want %v", tc.buildType, got, tc.want)
		}
	}
}
