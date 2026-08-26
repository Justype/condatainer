package publish

import (
	"context"
	"strings"
	"testing"

	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/project/lock"
)

// buildPlan vendors the artifacts, pins the named ones, and plans a push.
func buildPlan(t *testing.T, root string, l *lock.Lock, opts Options) *Plan {
	t.Helper()
	verified, problems := lock.Verify(root, l)
	if len(problems) > 0 {
		t.Fatalf("fixture does not verify: %v", problems)
	}
	plan, err := Build(context.Background(), root, l, verified, opts)
	if err != nil {
		t.Fatalf("Build: %v", err)
	}
	return plan
}

func lockWith(t *testing.T, root string, push string) *lock.Lock {
	t.Helper()
	l := lock.New()
	l.OCI = lock.OCI{Push: push}
	return l
}

func TestBuildDefaultsToSelectionsAndAddsClosureOnDemand(t *testing.T) {
	root := projectRoot(t)
	dep, depManifest := vendorData(t, root, "cutadapt/5.0", "echo cutadapt\n")
	pin, _ := vendorData(t, root, "grch38/genome/gencode49", "echo index\n", meta.Dependency{
		Name:     depManifest.Name,
		Type:     depManifest.Type,
		Identity: depManifest.Keys.Identity,
		Equiv:    depManifest.Keys.Equiv,
		Role:     "data",
	})

	l := lockWith(t, root, "ghcr.io/lab/p/cnt")
	l.Pins["grch38/genome/gencode49"] = lock.PinEntry{Artifact: pin}

	// Restore prunes a closure node whenever its dependent is fetched, so the
	// build dependency is reachable only on the path a push exists to avoid.
	plan := buildPlan(t, root, l, Options{})
	if len(plan.Steps) != 1 || plan.Steps[0].Artifact != pin {
		t.Fatalf("the default set is not the pins alone: %#v", plan.Steps)
	}
	if !plan.Steps[0].Pinned {
		t.Error("a pin is not marked selected, so it would lose its plain tag")
	}

	plan = buildPlan(t, root, l, Options{Closure: true})
	if len(plan.Steps) != 2 {
		t.Fatalf("--closure did not add the build dependency: %#v", plan.Steps)
	}
	for _, step := range plan.Steps {
		if step.Artifact == dep && step.Pinned {
			t.Error("a build dependency is marked selected, so it would take the plain name tag")
		}
		if step.Artifact == dep && len(step.Tags) != 1 {
			t.Errorf("a build dependency should carry only its qualified tag, got %v", step.Tags)
		}
	}
}

// A push with no destination cannot be planned, and the message has to name the
// command that fixes it.
func TestBuildRefusesWithoutADestination(t *testing.T) {
	root := projectRoot(t)
	l := lock.New()
	verified, _ := lock.Verify(root, l)
	_, err := Build(context.Background(), root, l, verified, Options{})
	if err == nil || !strings.Contains(err.Error(), "project registry set") {
		t.Errorf("err = %v, want one naming `project registry set`", err)
	}
}

// --registry publishes elsewhere without the lock having recorded it.
func TestBuildAcceptsAMirrorOverride(t *testing.T) {
	root := projectRoot(t)
	l := lock.New()
	verified, _ := lock.Verify(root, l)
	plan, err := Build(context.Background(), root, l, verified, Options{Repository: "registry.invalid/mirror"})
	if err != nil {
		t.Fatalf("Build: %v", err)
	}
	if plan.Repository != "registry.invalid/mirror" {
		t.Errorf("repository = %q", plan.Repository)
	}
}

// A public endpoint refuses an undeclared app, and the refusal has to reach the
// plan's problems rather than being silently dropped from the set.
func TestBuildRefusesWhatTheEndpointWillNotTake(t *testing.T) {
	root := projectRoot(t)
	app := vendorApp(t, root, "star/2.7.11b", "echo star\n")

	l := lockWith(t, root, "ghcr.io/lab/p/cnt")
	l.Pins["star/2.7.11b"] = lock.PinEntry{Artifact: app}

	plan := buildPlan(t, root, l, Options{})
	if len(plan.Steps) != 1 || plan.Steps[0].Disposition != Refused {
		t.Fatalf("an undeclared app was not refused: %#v", plan.Steps)
	}
	if plan.Complete() {
		t.Error("a refusal did not reach the plan's problems")
	}

	// Declared, and it publishes. The type rule is a default for an unanswered
	// question, not an absolute.
	l.OCI.Audience = "restricted"
	if plan := buildPlan(t, root, l, Options{}); plan.Steps[0].Disposition != Upload {
		t.Errorf("a restricted endpoint refused an app: %#v", plan.Steps[0])
	}
}

// Two pins sharing a manifest name: neither can hold the plain tag, since
// it is a claim about which one the project uses.
func TestAmbiguousNamesLoseThePlainTag(t *testing.T) {
	root := projectRoot(t)
	one, _ := vendorData(t, root, "grch38/genome/gencode49", "echo one\n")
	two, _ := vendorData(t, root, "grch38/genome/gencode49", "echo two\n")
	if one == two {
		t.Fatal("the fixture produced one artifact, not two identities of one name")
	}

	l := lockWith(t, root, "ghcr.io/lab/p/cnt")
	l.Pins["grch38/genome/gencode49"] = lock.PinEntry{Artifact: one}
	l.Pins["path:overlays/gencode49.sqf"] = lock.PinEntry{Artifact: two}

	plan := buildPlan(t, root, l, Options{})
	dropped := dropAmbiguousPlainTags(plan)
	if len(dropped) != 2 {
		t.Fatalf("dropped = %v, want both pins", dropped)
	}
	for _, step := range plan.Steps {
		for _, tag := range step.Tags {
			if !strings.Contains(tag, "__") {
				t.Errorf("%s kept the ambiguous plain tag %q", step.Name, tag)
			}
		}
		if len(step.Tags) != 1 {
			t.Errorf("%s should keep exactly its qualified tag, got %v", step.Name, step.Tags)
		}
	}
}

// A tag whose spelling depended on the rest of the set would be unfindable from
// the artifact alone, which is what lets a re-push skip an upload. So a
// collision refuses rather than lengthening the prefix.
func TestTagCollisionRefuses(t *testing.T) {
	plan := &Plan{Steps: []Step{
		{Artifact: "provenance/a", Tags: []string{"x--1.0__aaaaaaaaaaaa"}},
		{Artifact: "provenance/b", Tags: []string{"x--1.0__aaaaaaaaaaaa"}},
	}}
	if err := checkTagCollisions(plan); err == nil {
		t.Error("two artifacts composing one tag were not refused")
	}
}
