package catalog

import (
	"errors"
	"strings"
	"testing"
)

// A recipe may restate the type it already has, but a declaration that disagrees
// is a silent lie today — it reads as the path-based default — so it is rejected.
func TestValidateDeclaredType(t *testing.T) {
	cases := []struct {
		name     string
		declared string
		typ      Type
		wantErr  string
	}{
		{"agrees/app", "app", TypeApp, ""},
		{"agrees/data", "data", TypeData, ""},
		{"a def restating its own type", "os", TypeOS, ""},
		{"nothing declared", "", TypeApp, ""},
		{"env is never declarable", "env", TypeApp, "overlay freeze"},
		{"disagrees", "data", TypeApp, "declares #TYPE:data but is app"},
		{"nonsense", "container", TypeApp, "declares #TYPE:container"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			r := &Recipe{DeclaredType: tc.declared}
			r.Name = "x/1.0"
			r.Type = tc.typ
			err := r.Validate()
			switch {
			case tc.wantErr == "" && err != nil:
				t.Fatalf("unexpected error: %v", err)
			case tc.wantErr == "":
				return
			case err == nil:
				t.Fatalf("expected an error mentioning %q", tc.wantErr)
			case !errors.Is(err, ErrInvalidRecipe):
				t.Fatalf("error is not ErrInvalidRecipe: %v", err)
			case !strings.Contains(err.Error(), tc.wantErr):
				t.Fatalf("error %q does not mention %q", err, tc.wantErr)
			}
		})
	}
}
