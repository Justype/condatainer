package meta

import (
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
)

func envRuntime(prefix string) Runtime {
	return Runtime{
		SchemaVersion: SchemaVersion,
		Name:          EnvName,
		Type:          catalog.TypeEnv,
		Platform:      Platform{Arch: "amd64"},
		Prefix:        prefix,
	}
}

// An environment records EnvPrefix or nothing. Recording it is what makes PATH
// and {prefix} resolve as they did for the writable .img, and what makes two
// environments — or an environment and the .img it came from — collide through
// the ordinary prefix check rather than needing a rule of their own.
func TestEnvRuntimePrefix(t *testing.T) {
	if err := ValidateRuntime(envRuntime(EnvPrefix)); err != nil {
		t.Errorf("%s was rejected: %v", EnvPrefix, err)
	}
	if err := ValidateRuntime(envRuntime("")); err != nil {
		t.Errorf("an environment with no conda prefix was rejected: %v", err)
	}
	err := ValidateRuntime(envRuntime("/cnt/rnaseq/1.0"))
	if err == nil {
		t.Fatal("an environment claimed a /cnt/<name> subtree it does not own")
	}
	if !strings.Contains(err.Error(), EnvPrefix) {
		t.Errorf("refusal does not name the prefix it allows: %v", err)
	}
}
