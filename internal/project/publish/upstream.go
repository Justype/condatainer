// Package publish resolves where a project's artifacts already live in an OCI
// registry, and where it puts the ones that do not.
//
// It knows about locks, collections, and the registry transport; none of them
// knows about it. In particular internal/project/lock gains no registry import,
// which is what keeps `project validate` checkout-local by construction rather
// than by discipline.
package publish

import (
	"context"
	"strings"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/project/lock"
	"github.com/Justype/condatainer/internal/registry"
)

// Upstream finds the artifacts among these vendored directories that a recipe
// collection already publishes at the exact identity the lock pins, and returns
// one remote per artifact that matched.
//
// It is a resolve-and-compare, never a search: the name and the identity are
// both already known, so this asks one coordinate per declared endpoint and
// verifies the answer. Nothing enumerates tags.
//
// Everything about it is best effort. No network, no configured source, no
// declared endpoint, no match, an unreachable registry, an unauthorized one —
// each records nothing and none is an error. Locking has to work offline, and an
// absent remote costs a rebuild rather than a failure.
func Upstream(ctx context.Context, root string, artifacts []string, cat catalog.Catalog) map[string][]lock.Remote {
	if len(artifacts) == 0 || len(cat) == 0 {
		return nil
	}
	log := logging.FromContext(ctx)
	out := make(map[string][]lock.Remote)
	for _, artifact := range artifacts {
		entry, err := lock.ReadEntry(root, artifact)
		if err != nil {
			log.Debug("cannot read a vendored artifact while looking for upstream copies", "artifact", artifact, "err", err)
			continue
		}
		remote, ok := upstreamOf(ctx, entry, cat)
		if !ok {
			continue
		}
		out[artifact] = []lock.Remote{remote}
		log.Debug("upstream publishes this exact identity", "name", entry.Manifest.Name,
			"repository", remote.Repository, "digest", remote.ManifestDigest)
	}
	if len(out) == 0 {
		return nil
	}
	return out
}

// upstreamOf walks one artifact's collection endpoints in declared order and
// returns the first that advertises this exact identity.
func upstreamOf(ctx context.Context, entry *lock.Entry, cat catalog.Catalog) (lock.Remote, bool) {
	source := collectionOf(entry.Manifest, cat)
	if source == nil {
		return lock.Remote{}, false
	}
	repo, tag, err := registry.PullReference(entry.Manifest.Type, entry.Manifest.Name)
	if err != nil {
		return lock.Remote{}, false
	}
	log := logging.FromContext(ctx)
	for _, endpoint := range source.Desc.OCI.Pull {
		// A public endpoint never carries an app, so asking is a round trip that
		// can only answer no. The same skip build.tryPrebuilt makes.
		if source.Desc.OCI.Audience == string(registry.Public) && entry.Manifest.Type == catalog.TypeApp {
			continue
		}
		desc, annotations, err := registry.ResolveArtifact(ctx, endpoint, repo, tag)
		if err != nil {
			log.Debug("endpoint did not answer for this artifact", "endpoint", endpoint, "name", entry.Manifest.Name, "err", err)
			continue
		}
		// Identity, not equivalence. A lock remote is an address for one exact
		// build: an equivalent artifact at that digest is a different build
		// wearing the right label, and restore would spend the bytes fetching it
		// and then reject it against the lock.
		if got := registry.Identity(annotations); got != entry.Identity {
			log.Debug("endpoint publishes a different build", "endpoint", endpoint,
				"name", entry.Manifest.Name, "published", describe(got), "locked", entry.Identity.Digest())
			continue
		}
		return lock.Remote{
			Repository:     strings.TrimRight(registry.TrimBaseScheme(endpoint), "/") + "/" + repo,
			ManifestDigest: desc.Digest.String(),
		}, true
	}
	return lock.Remote{}, false
}

// collectionOf finds the configured source an artifact was built from, by the
// collection repository its manifest recorded.
//
// A collection that is unreadable, stale, or declares no pull endpoint answers
// nothing: its registry defaults cannot be trusted, and a coordinate guessed
// from a broken descriptor would be written into a tracked file.
func collectionOf(m meta.Manifest, cat catalog.Catalog) *catalog.Source {
	want := strings.TrimRight(strings.TrimSpace(m.Build.Source), "/")
	if want == "" {
		return nil
	}
	var found *catalog.Source
	for _, source := range cat {
		if strings.TrimRight(strings.TrimSpace(source.Desc.Source), "/") != want {
			continue
		}
		if found != nil {
			// Two collections claiming one repository: nothing here can say which
			// published the artifact, and picking either would be a guess.
			return nil
		}
		found = source
	}
	if found == nil || found.DescriptorErr != nil || found.Err != nil || found.Stale {
		return nil
	}
	if len(found.Desc.OCI.Pull) == 0 {
		return nil
	}
	return found
}

func describe(k meta.KeyRef) string {
	if k.Empty() {
		return "no key"
	}
	return k.Digest()
}
