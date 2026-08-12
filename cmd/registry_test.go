package cmd

import (
	"bytes"
	"errors"
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/registry"
	"github.com/spf13/cobra"
)

func TestInferPushSource(t *testing.T) {
	good := &catalog.Source{Name: "lab", Desc: catalog.Descriptor{
		Repository: "https://example.invalid/recipes/",
		OCI:        catalog.OCI{Push: "registry.invalid/lab", Visibility: "internal"},
	}}
	for _, tc := range []struct {
		name    string
		cat     catalog.Catalog
		want    *catalog.Source
		wantErr string
	}{
		{"unique", catalog.Catalog{good}, good, ""},
		{"absent", nil, nil, "matches 0"},
		{"ambiguous", catalog.Catalog{good, &catalog.Source{Name: "other", Desc: good.Desc}}, nil, "matches 2"},
		{"invalid descriptor", catalog.Catalog{&catalog.Source{Name: "bad", Desc: good.Desc, DescriptorErr: errors.New("bad json")}}, nil, "invalid source descriptor"},
		{"stale", catalog.Catalog{&catalog.Source{Name: "old", Desc: good.Desc, Stale: true}}, nil, "unavailable or stale"},
		{"no push", catalog.Catalog{&catalog.Source{Name: "nopush", Desc: catalog.Descriptor{Repository: good.Desc.Repository}}}, nil, "no OCI push endpoint"},
	} {
		t.Run(tc.name, func(t *testing.T) {
			got, err := inferPushSource("https://example.invalid/recipes", tc.cat)
			if tc.wantErr != "" {
				if err == nil || !strings.Contains(err.Error(), tc.wantErr) {
					t.Fatalf("error = %v, want %q", err, tc.wantErr)
				}
				return
			}
			if err != nil || got != tc.want {
				t.Fatalf("source = %p, error = %v", got, err)
			}
		})
	}
}

func TestRegistryCommandHasThePlannedSurface(t *testing.T) {
	cmd := newRegistryCommand()
	want := map[string]bool{
		"push": false, "pull": false, "tags": false,
		"resolve": false, "login": false, "logout": false,
	}
	for _, child := range cmd.Commands() {
		if _, ok := want[child.Name()]; ok {
			want[child.Name()] = true
		}
	}
	for name, found := range want {
		if !found {
			t.Errorf("registry command is missing %q", name)
		}
	}
}

func TestRegistryPushCompletionStopsAfterArtifact(t *testing.T) {
	got, directive := registryPushCompletion(nil, []string{"hello/1.0"}, "")
	if len(got) != 0 || directive != cobra.ShellCompDirectiveNoFileComp {
		t.Fatalf("completion = %v, %v", got, directive)
	}
}

func TestCompletionScriptsRecognizeRegistryPushFZF(t *testing.T) {
	tests := []struct {
		name, input, marker string
		process             func(string) string
	}{
		{
			name: "bash", input: `__condatainer_debug "The completions are: ${out}"`,
			marker: `"$sub" == "registry" && "${words[2]}" == "push"`, process: postProcessBashCompletion,
		},
		{
			name: "zsh", input: `__condatainer_debug "completions: ${out}"`,
			marker: `"$sub" == "registry" && "${words[3]}" == "push"`, process: postProcessZshCompletion,
		},
		{
			name: "fish", input: `set -l results (eval $requestComp 2> /dev/null)`,
			marker: `test "$args[2]" = registry; and test "$args[3]" = push`, process: postProcessFishCompletion,
		},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			got := tc.process(tc.input)
			if !strings.Contains(got, tc.marker) || !strings.Contains(got, "Select artifact") {
				t.Fatalf("processed completion lacks registry push fzf hook:\n%s", got)
			}
		})
	}
}

func TestRegistryHelpersValidateInputs(t *testing.T) {
	if _, err := requireRegistryBase(""); err == nil || !strings.Contains(err.Error(), "--registry") {
		t.Errorf("missing registry error = %v", err)
	}
	if got, err := requireRegistryBase("oci://ghcr.io/lab/cnt/"); err != nil || got != "ghcr.io/lab/cnt" {
		t.Errorf("requireRegistryBase = (%q, %v)", got, err)
	}
	if _, err := parseVisibility("private"); err == nil {
		t.Error("parseVisibility accepted an unknown value")
	}
	if got, err := parseVisibility("INTERNAL"); err != nil || got != registry.Internal {
		t.Errorf("parseVisibility = (%q, %v)", got, err)
	}
}

func TestAddressPlacementName(t *testing.T) {
	digest := "sha256:" + strings.Repeat("a", 64)
	for _, tc := range []struct {
		name, selector, title, want string
	}{
		{"hello", "1.0", "hello/1.0", "hello/1.0"},
		{"hello", "1.0", "different/1.0", "hello/1.0"},
		{"ubuntu24/base", digest, "ubuntu24/base", "ubuntu24/base"},
		{"hello", digest, "hello/1.0", "hello/1.0"},
		{"hello", digest, "", ""},
	} {
		if got := addressPlacementName(tc.name, tc.selector, tc.title); got != tc.want {
			t.Errorf("addressPlacementName(%q, %q, %q) = %q, want %q",
				tc.name, tc.selector, tc.title, got, tc.want)
		}
	}
}

func TestRegistryPullDestinationPrecedenceAndValidation(t *testing.T) {
	t.Setenv("XDG_DATA_HOME", t.TempDir())
	images := filepath.Join(config.GetUserDataDir(), "images")

	dest, err := registryPullDestination("renamed/2.0", "", "hello", "1.0", "hello/1.0", ".sqf")
	if err != nil || dest != filepath.Join(images, "renamed--2.0.sqf") {
		t.Fatalf("--name destination = (%q, %v)", dest, err)
	}
	dest, err = registryPullDestination("", filepath.Join(t.TempDir(), "custom"), "hello", "1.0", "hello/1.0", ".sqf")
	if err != nil || filepath.Base(dest) != "custom.sqf" {
		t.Fatalf("--prefix destination = (%q, %v)", dest, err)
	}
	if _, err := registryPullDestination("", "bad--name", "hello", "1.0", "hello/1.0", ".sqf"); err == nil {
		t.Error("--prefix accepted the reserved name separator")
	}
	if _, err := registryPullDestination("", "base.sqf", "ubuntu24/base", "latest", "ubuntu24/base", ".sif"); err == nil {
		t.Error("--prefix accepted the wrong image extension")
	}
	dest, err = registryPullDestination("", filepath.Join(t.TempDir(), "custom.v1"), "hello", "1.0", "hello/1.0", ".sqf")
	if err != nil || filepath.Base(dest) != "custom.v1.sqf" {
		t.Errorf("extensionless dotted prefix = (%q, %v)", dest, err)
	}
	if _, err := registryPullDestination("bad--name", "", "hello", "1.0", "hello/1.0", ".sqf"); err == nil {
		t.Error("--name accepted the reserved name separator")
	}
	if _, err := registryPullDestination("", "", "hello", "sha256:"+strings.Repeat("b", 64), "", ".sqf"); err == nil {
		t.Error("versioned digest without a title guessed an install name")
	}
}

func TestRegistryPasswordFromStdin(t *testing.T) {
	cmd := newRegistryCommand()
	cmd.SetIn(bytes.NewBufferString("secret\n"))
	got, err := registryPassword(cmd, &registryOptions{passwordStdin: true})
	if err != nil || got != "secret" {
		t.Fatalf("registryPassword = (%q, %v)", got, err)
	}
	if _, err := registryPassword(cmd, &registryOptions{passwordStdin: true, password: "also"}); err == nil {
		t.Error("registryPassword accepted two password sources")
	}
}

func TestRegistryPullRefusesABareNameBeforeNetwork(t *testing.T) {
	cmd := newRegistryCommand()
	cmd.SilenceErrors = true
	cmd.SetArgs([]string{"pull", "hello", "--registry", "registry.invalid/lab"})
	err := cmd.Execute()
	if err == nil || !strings.Contains(err.Error(), "use create") {
		t.Fatalf("bare-name pull error = %v", err)
	}
}
