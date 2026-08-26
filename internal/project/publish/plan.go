package publish

import (
	"context"
	"fmt"
	"sort"
	"strings"

	"github.com/Justype/condatainer/internal/project/lock"
	"github.com/Justype/condatainer/internal/registry"
)

// Disposition is what a push will do with one artifact.
type Disposition string

const (
	// Upload sends the payload.
	Upload Disposition = "upload"
	// Present skips it: the destination already holds this exact identity, so
	// the tags are re-applied and the digest recorded without a transfer.
	Present Disposition = "present"
	// Served skips it: a collection already publishes this exact identity, so
	// the location is recorded and nothing is uploaded.
	Served Disposition = "served"
	// Refused will not be published, and Reason says why.
	Refused Disposition = "refused"
)

// Step is one artifact's place in a push.
type Step struct {
	// Artifact is the vendored directory, relative to the lock directory.
	Artifact string `json:"artifact"`
	Name     string `json:"name"`
	Identity string `json:"identity"`
	// Pinned reports that a pin points at this artifact rather than only
	// a dependency edge. It decides the plain name tag.
	Pinned      bool        `json:"pinned"`
	Disposition Disposition `json:"disposition"`
	// Path is the local artifact to upload, empty for anything not uploaded.
	Path string `json:"path,omitempty"`
	// Destination is the project-relative path a `path:` pin lives at.
	// Empty means the artifact is addressed by identity in an images root.
	Destination string `json:"destination,omitempty"`
	// Tags are what this artifact publishes under, canonical first.
	Tags []string `json:"tags,omitempty"`
	// Remote is the location to record: what an upload will create, or where an
	// existing copy was found.
	Remote *lock.Remote `json:"remote,omitempty"`
	// Reason explains a refusal.
	Reason string `json:"reason,omitempty"`
	// Channels are a Conda build's recorded channels, reported and never judged:
	// no allowlist could be maintained honestly and no licence string could be
	// interpreted safely, so the operator is shown what went into the solve.
	Channels []string `json:"channels,omitempty"`
}

// Plan is what a push will publish, and everything that would prevent it.
type Plan struct {
	Repository string   `json:"repository"`
	Audience   string   `json:"audience"`
	Source     string   `json:"source,omitempty"`
	Steps      []Step   `json:"steps"`
	Problems   []string `json:"problems,omitempty"`
}

// Complete reports whether the plan can run as computed.
func (p *Plan) Complete() bool { return len(p.Problems) == 0 }

// Uploads reports how many artifacts would actually be transferred.
func (p *Plan) Uploads() int {
	var n int
	for _, step := range p.Steps {
		if step.Disposition == Upload {
			n++
		}
	}
	return n
}

// Options tunes what a push publishes.
type Options struct {
	// All includes pins an upstream collection already serves. The default
	// skips them: restore tries the recorded upstream location first and finds
	// it there, so a second copy buys nothing but the bytes.
	//
	// It is for the project that must not depend on a collection's retention —
	// a paper's artifacts outliving the lab's recipe repository — where holding
	// a second copy is exactly the point.
	All bool
	// Closure adds every vendored build dependency. Restore prunes a closure
	// node whenever its dependent is adopted or fetched, so these are reachable
	// only on the rebuild path, which is the path a push exists to avoid.
	Closure bool
	// Repository overrides the recorded destination, for publishing a mirror.
	// It does not change what the lock records as this project's destination.
	Repository string
}

// Build computes what a push would do, without contacting anything.
//
// Everything network-shaped is left to the caller: this decides the set, the
// tags, and the refusals from the checkout alone, so `--dry-run` can answer the
// expensive questions cheaply and a real push can refuse before its first byte.
func Build(ctx context.Context, root string, l *lock.Lock, verified *lock.Verified, opts Options) (*Plan, error) {
	destination := strings.TrimSpace(opts.Repository)
	if destination == "" {
		destination = strings.TrimSpace(l.OCI.Push)
	}
	if destination == "" {
		return nil, fmt.Errorf("this project records no registry; set one with `condatainer project registry set`")
	}
	audience := registry.Audience(strings.ToLower(strings.TrimSpace(l.OCI.Audience)))
	if audience == "" {
		audience = registry.Public
	}

	plan := &Plan{Repository: destination, Audience: string(audience), Source: l.Source}
	pinned := l.PinnedArtifacts()
	destinations := projectPaths(l)

	for _, artifact := range sortedArtifacts(verified) {
		entry := verified.Entries[artifact]
		isPinned := pinned[artifact]
		if !isPinned && !opts.Closure {
			continue
		}
		step := Step{
			Artifact:    artifact,
			Name:        entry.Manifest.Name,
			Identity:    entry.Identity.Digest(),
			Pinned:      isPinned,
			Destination: destinations[artifact],
			Channels:    entry.Manifest.Build.Channels,
		}
		plan.Steps = append(plan.Steps, planStep(step, entry, audience))
	}
	if len(plan.Steps) == 0 {
		plan.Problems = append(plan.Problems, "nothing to publish")
	}
	for _, step := range plan.Steps {
		if step.Disposition == Refused {
			plan.Problems = append(plan.Problems, fmt.Sprintf("%s: %s", step.Name, step.Reason))
		}
	}
	if err := checkTagCollisions(plan); err != nil {
		plan.Problems = append(plan.Problems, err.Error())
	}
	return plan, nil
}

// planStep decides one artifact's tags and whether the endpoint will take it.
func planStep(step Step, entry *lock.Entry, audience registry.Audience) Step {
	if err := audience.Accepts(entry.Manifest); err != nil {
		step.Disposition, step.Reason = Refused, err.Error()
		return step
	}
	tags, err := registry.ProjectTags(entry.Manifest, step.Pinned)
	if err != nil {
		step.Disposition, step.Reason = Refused, err.Error()
		return step
	}
	step.Disposition, step.Tags = Upload, tags
	return step
}

// checkTagCollisions refuses a set in which two artifacts compose one tag with
// different identities.
//
// It refuses rather than lengthening the identity prefix, which is what the
// store does for a filename. A 12-hex collision does not happen; a tag whose
// spelling depended on the rest of the push set would, and it would make the
// tag unfindable from the artifact alone — which is what makes a re-push able
// to skip an upload.
func checkTagCollisions(plan *Plan) error {
	seen := make(map[string]string, len(plan.Steps)*2)
	for _, step := range plan.Steps {
		for _, tag := range step.Tags {
			if other, taken := seen[tag]; taken && other != step.Artifact {
				return fmt.Errorf("%s and %s both publish as %q", other, step.Artifact, tag)
			}
			seen[tag] = step.Artifact
		}
	}
	return nil
}

// dropAmbiguousPlainTags removes the unqualified name tag where two artifacts in
// one set share a manifest name, and reports what it dropped.
//
// The plain tag is a human handle that says "this is what the project uses", and
// two artifacts cannot both be that. The qualified tags still identify each, and
// nothing addresses a published artifact by tag anyway.
func dropAmbiguousPlainTags(plan *Plan) []string {
	byName := map[string]int{}
	for _, step := range plan.Steps {
		if step.Pinned && step.Disposition != Refused {
			byName[step.Name]++
		}
	}
	var dropped []string
	for i, step := range plan.Steps {
		if !step.Pinned || byName[step.Name] < 2 {
			continue
		}
		var kept []string
		for _, tag := range step.Tags {
			if strings.Contains(tag, "__") {
				kept = append(kept, tag)
			}
		}
		plan.Steps[i].Tags = kept
		dropped = append(dropped, step.Name)
	}
	return dropped
}

// projectPaths maps each artifact a `path:` pin points at to where the
// project keeps it. Such an artifact is never in an images root — the project
// owns the file — so a push has to read it from where the project mounts it.
func projectPaths(l *lock.Lock) map[string]string {
	out := map[string]string{}
	for request, pin := range l.Pins {
		if destination, ok := strings.CutPrefix(request, lock.PathPrefix); ok {
			out[pin.Artifact] = destination
		}
	}
	return out
}

func sortedArtifacts(verified *lock.Verified) []string {
	out := make([]string, 0, len(verified.Entries))
	for artifact := range verified.Entries {
		if verified.Reachable[artifact] {
			out = append(out, artifact)
		}
	}
	sort.Strings(out)
	return out
}
