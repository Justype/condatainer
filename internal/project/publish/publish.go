package publish

import (
	"context"
	"errors"
	"fmt"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/project"
	"github.com/Justype/condatainer/internal/project/lock"
	"github.com/Justype/condatainer/internal/registry"
)

// ErrIncomplete reports that a push did not publish everything it planned.
var ErrIncomplete = errors.New("push is incomplete")

// Report is what a push did.
type Report struct {
	Repository string   `json:"repository"`
	Published  []Step   `json:"published,omitempty"`
	Failures   []string `json:"failures,omitempty"`
	// AmbiguousNames are pins that share a manifest name, so no plain name
	// tag could say which one the project uses.
	AmbiguousNames []string `json:"ambiguous_names,omitempty"`
}

// Refine asks the network the two questions the checkout cannot answer: whether
// a collection already serves an artifact, and whether the destination already
// holds it.
//
// Both only ever *remove* work. Neither can add an artifact to the set or change
// what it publishes as, so a refusal to answer — an unreachable registry, an
// unauthorized one — leaves the plan as the checkout computed it and the push
// proceeds. That is deliberate: the cost of guessing wrong is one redundant
// upload, and the cost of failing here is a push nobody can complete offline.
func Refine(ctx context.Context, root string, plan *Plan, cat catalog.Catalog, opts Options) {
	if !opts.All {
		markUpstreamServed(ctx, root, plan, cat)
	}
	markAlreadyPublished(ctx, plan)
}

// markUpstreamServed drops artifacts a collection already publishes at this
// exact identity, recording where instead of uploading a second copy.
//
// Evaluated now rather than trusted from the lock: a recorded location that no
// longer resolves does not count, so the default set republishes it rather than
// leaving the project pointing at something gone.
func markUpstreamServed(ctx context.Context, root string, plan *Plan, cat catalog.Catalog) {
	var candidates []string
	for _, step := range plan.Steps {
		if step.Disposition == Upload {
			candidates = append(candidates, step.Artifact)
		}
	}
	served := Upstream(ctx, root, candidates, cat)
	for i, step := range plan.Steps {
		remotes, ok := served[step.Artifact]
		if !ok || len(remotes) == 0 {
			continue
		}
		plan.Steps[i].Disposition = Served
		plan.Steps[i].Remote = &remotes[0]
		plan.Steps[i].Tags = nil
	}
}

// markAlreadyPublished skips an upload the destination already holds.
//
// The recorded digest is asked first: it resolving is proof the content is
// there, in one request, and it is the address restore will actually use. The
// qualified tag is the fallback, and it is the repair path — a lock whose
// remotes were lost finds its own artifacts again because the tag is derived
// from the identity.
//
// Tags are re-applied either way. A manifest that is present but has lost its
// retention anchor is exactly the state this exists to repair.
func markAlreadyPublished(ctx context.Context, plan *Plan) {
	log := logging.FromContext(ctx)
	base, repo, err := registry.SplitCoordinate(plan.Repository)
	if err != nil {
		return
	}
	for i, step := range plan.Steps {
		if step.Disposition != Upload || len(step.Tags) == 0 {
			continue
		}
		desc, annotations, err := registry.ResolveArtifact(ctx, base, repo, step.Tags[0])
		if err != nil {
			log.Debug("destination does not hold this artifact yet", "name", step.Name, "err", err)
			continue
		}
		if got := registry.Identity(annotations); got.Digest() != step.Identity {
			// A content-keyed tag holding different content is a corrupted
			// package, not something to overwrite.
			plan.Steps[i].Disposition = Refused
			plan.Steps[i].Reason = fmt.Sprintf("%s already publishes %s, which is not this artifact",
				step.Tags[0], describe(got))
			plan.Problems = append(plan.Problems, fmt.Sprintf("%s: %s", step.Name, plan.Steps[i].Reason))
			continue
		}
		plan.Steps[i].Disposition = Present
		plan.Steps[i].Remote = &lock.Remote{
			Repository:     strings.TrimRight(registry.TrimBaseScheme(plan.Repository), "/"),
			ManifestDigest: desc.Digest.String(),
		}
	}
}

// Run publishes the plan and records where each artifact landed.
//
// One lock transaction per artifact, not one at the end: a project push is tens
// of gigabytes over a link that drops, and losing the record of nine completed
// uploads because the tenth timed out is worse than nine small rewrites. Each
// transaction is atomic on its own, and a location recorded for an artifact that
// did publish stays true whatever happens afterwards.
func Run(ctx context.Context, root string, plan *Plan, opts Options) (*Report, error) {
	log := logging.FromContext(ctx)
	report := &Report{Repository: plan.Repository, AmbiguousNames: plan.Ambiguous}
	base, repo, err := registry.SplitCoordinate(plan.Repository)
	if err != nil {
		return report, err
	}
	audience := registry.Audience(plan.Audience)

	for _, step := range plan.Steps {
		switch step.Disposition {
		case Refused:
			report.Failures = append(report.Failures, fmt.Sprintf("%s: %s", step.Name, step.Reason))
			continue
		case Served, Present:
			if err := record(root, step); err != nil {
				report.Failures = append(report.Failures, fmt.Sprintf("%s: %v", step.Name, err))
				continue
			}
			report.Published = append(report.Published, step)
			continue
		}

		path, err := locate(root, step)
		if err != nil {
			report.Failures = append(report.Failures, fmt.Sprintf("%s: %v", step.Name, err))
			continue
		}
		log.Info("publishing", "kind", "note", "name", step.Name, "tags", strings.Join(step.Tags, ", "))
		desc, err := registry.Publish(ctx, registry.PublishRequest{
			Path:      path,
			Base:      base,
			Audience:  audience,
			Placement: &registry.Placement{Repo: repo, Tags: step.Tags},
			Source:    plan.Source,
		})
		if err != nil {
			report.Failures = append(report.Failures, fmt.Sprintf("%s: %v", step.Name, err))
			if ctx.Err() != nil {
				break
			}
			continue
		}
		step.Path = path
		step.Remote = &lock.Remote{
			Repository:     strings.TrimRight(registry.TrimBaseScheme(plan.Repository), "/"),
			ManifestDigest: desc.Digest.String(),
		}
		if err := record(root, step); err != nil {
			report.Failures = append(report.Failures, fmt.Sprintf("%s: %v", step.Name, err))
			continue
		}
		report.Published = append(report.Published, step)
	}

	if len(report.Failures) > 0 {
		return report, ErrIncomplete
	}
	return report, nil
}

// record appends one artifact's location to the lock, atomically, reloading
// first so a concurrent writer's pins are not lost.
func record(root string, step Step) error {
	if step.Remote == nil {
		return nil
	}
	current, err := lock.Load(root)
	if err != nil {
		return err
	}
	if err := current.AddRemote(step.Artifact, *step.Remote); err != nil {
		return err
	}
	return lock.Publish(root, current)
}

// locate finds the artifact to upload on this machine, at the exact identity the
// lock pins.
//
// Identity and nothing looser: an equivalent artifact is a legitimate answer for
// what this machine runs and a lie about what the lock names, so publishing one
// would put a substitution behind an address the lock claims is exact.
func locate(root string, step Step) (string, error) {
	entry, err := lock.ReadEntry(root, step.Artifact)
	if err != nil {
		return "", err
	}
	// A project path is where the project *keeps* this artifact, so it is the
	// only place to read it from. A copy of the same identity in an images root
	// is a different file the project does not mount.
	if step.Destination != "" {
		at := filepath.Join(root, filepath.FromSlash(step.Destination))
		candidate, ok := project.LookupAt(at, entry.Manifest.Name, entry.Manifest.Keys, project.MatchIdentity)
		if !ok {
			return "", fmt.Errorf("%s is not at %s at %s; run `condatainer project select-match identity`, then `project restore`",
				entry.Manifest.Name, step.Destination, entry.Identity.Digest())
		}
		return candidate.Path, nil
	}
	candidate, ok := project.LookupLocal(entry.Manifest.Name, entry.Manifest.Keys, project.MatchIdentity, nil)
	if !ok {
		return "", fmt.Errorf("%s is not installed at %s; run `condatainer project select-match identity`, then `project restore`",
			entry.Manifest.Name, entry.Identity.Digest())
	}
	return candidate.Path, nil
}
