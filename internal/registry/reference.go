package registry

import (
	"fmt"
	"regexp"
	"strings"

	"github.com/opencontainers/go-digest"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/capsule"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/image"
)

// RollingTag moves with every push of a version-less artifact, so a puller that
// does not know a build date still resolves the current one.
const RollingTag = "latest"

// dateTagLayout formats a version-less artifact's build-date tag.
const dateTagLayout = "20060102"

var (
	ociRepoSegmentPattern = regexp.MustCompile(`^[a-z0-9]+(?:[._-][a-z0-9]+)*$`)
	ociTagPattern         = regexp.MustCompile(`^[A-Za-z0-9_][A-Za-z0-9_.-]{0,127}$`)
)

// TrimBaseScheme strips an optional oci:// or oras:// scheme and any trailing
// slash from a registry *base*. A base already includes the owner and prefix;
// repository paths are appended to it.
//
// This is not [Normalize], which expands a short *image reference* the way a
// container runtime does. Two unrelated jobs: this one removes a scheme a human
// typed, that one supplies a registry host that was left out.
func TrimBaseScheme(base string) string {
	base = strings.TrimSpace(base)
	for _, scheme := range []string{"oci://", "oras://"} {
		base = strings.TrimPrefix(base, scheme)
	}
	return strings.TrimRight(base, "/")
}

// FullRef joins a registry base, a repository path, and a tag or digest into one
// pullable reference. A digest joins with "@", a tag with ":".
func FullRef(base, repo, tag string) string {
	if strings.HasPrefix(tag, digestPrefix) {
		return TrimBaseScheme(base) + "/" + repo + "@" + tag
	}
	return TrimBaseScheme(base) + "/" + repo + ":" + tag
}

// digestPrefix is the only digest algorithm this package writes or accepts.
const digestPrefix = "sha256:"

// isVersionLess reports whether an artifact is addressed by build date and a
// rolling tag rather than by a version segment of its own name.
//
// One rule: **version-less means the name carries no version segment**, so there
// is nothing to tag with and a date stands in for it. What differs by type is how
// many segments a name needs before one of them is a version:
//
//   - a base never has one. Its trailing segment is a role — ubuntu24/base — and
//     it is rebuilt in place under that one name.
//   - an OS needs three. ubuntu24/build-essential is a name whose payload is
//     whatever the apt mirror served that day; ubuntu24/r/4.4.4 pins a version.
//   - an app or data recipe needs two. cellranger/9.0.1 pins its own download,
//     while a bare `myenv` pins nothing.
//
// The catalog parses every name as <name>/<version>, but that is an addressing
// convention rather than a claim about meaning: `build-essential` is a name, not
// a version. So the type is what separates ubuntu24/build-essential from
// hello/1.0, which are otherwise the same shape.
func isVersionLess(typ catalog.Type, name string) bool {
	segments := strings.Count(catalog.Normalize(name), "/") + 1
	switch typ {
	case catalog.TypeBase:
		return true
	case catalog.TypeOS:
		return segments < 3
	default:
		return segments < 2
	}
}

// ValidateIdentity rejects a name that cannot map reversibly onto an OCI
// repository and tag.
//
// Repository segments must already be lowercase. Nothing is down-cased here: two
// distinct names that differ only in case would collide into one repository, and
// silently publishing one over the other is worse than refusing both.
func ValidateIdentity(nameVersion string, versionLess bool) error {
	nv := catalog.Normalize(nameVersion)
	if nv == "" || strings.ContainsAny(nv, ":@") {
		return fmt.Errorf("registry identity %q is empty or carries a selector", nameVersion)
	}
	parts := strings.Split(nv, "/")
	repoParts := parts
	if !versionLess {
		if len(parts) < 2 {
			return fmt.Errorf("versioned artifact %q must be name/version", nameVersion)
		}
		tag := parts[len(parts)-1]
		if !ociTagPattern.MatchString(tag) || tag == "." || tag == ".." {
			return fmt.Errorf("identity tag segment %q is not OCI-safe", tag)
		}
		repoParts = parts[:len(parts)-1]
	}
	for _, segment := range repoParts {
		if segment == "." || segment == ".." || !ociRepoSegmentPattern.MatchString(segment) {
			return fmt.Errorf("identity repository segment %q must already be lowercase and OCI-safe", segment)
		}
	}
	return nil
}

// SplitPullSpec separates an artifact name from an optional selector. The forms
// are "name", "name:tag", and "name@sha256:<hex>".
//
// A colon counts as a tag separator only after the final slash, so a registry
// host with a port stays the separate base argument's problem.
func SplitPullSpec(spec string) (name, selector string, err error) {
	spec = strings.TrimSpace(spec)
	if at := strings.LastIndex(spec, "@"); at >= 0 {
		name, selector = catalog.Normalize(spec[:at]), spec[at+1:]
		if name == "" || digest.Digest(selector).Validate() != nil {
			return "", "", fmt.Errorf("invalid digest pull reference %q", spec)
		}
		return name, selector, nil
	}
	if slash, colon := strings.LastIndex(spec, "/"), strings.LastIndex(spec, ":"); colon > slash {
		name, selector = catalog.Normalize(spec[:colon]), spec[colon+1:]
		if name == "" || selector == "" {
			return "", "", fmt.Errorf("invalid tagged pull reference %q", spec)
		}
		return name, selector, nil
	}
	return catalog.Normalize(spec), "", nil
}

// PushReference derives the repository and ordered tags for publishing m. The
// first tag is canonical; a version-less artifact also gets [RollingTag].
//
// Tags are architecture-independent. A native artifact's tag resolves to an image
// index whose children carry per-arch platform descriptors, so pushing on each
// architecture builds one multi-arch tag.
func PushReference(m meta.Manifest) (repo string, tags []string, err error) {
	nv := catalog.Normalize(m.Name)
	versionLess := isVersionLess(m.Type, nv)
	if err := ValidateIdentity(nv, versionLess); err != nil {
		return "", nil, err
	}

	if versionLess {
		// No fallback to the push time: a date tag is an address, and one that
		// records when it was uploaded rather than when it was built is a lie
		// that only shows up as a wrong artifact months later.
		if m.Build.Created.IsZero() {
			return "", nil, fmt.Errorf("%s records no build time, so it has no date tag", m.Name)
		}
		return nv, []string{m.Build.Created.UTC().Format(dateTagLayout), RollingTag}, nil
	}

	idx := strings.LastIndex(nv, "/")
	return nv[:idx], []string{nv[idx+1:]}, nil
}

// PullReference derives the repository and single tag to fetch nameVersion,
// mirroring PushReference's canonical tag.
//
// A puller does not know a version-less artifact's build date, so it resolves
// RollingTag; an explicit date or digest arrives through SplitPullSpec instead.
// typ comes from the recipe rather than from a manifest, because the whole point
// of pulling is that the artifact is not here to read.
func PullReference(typ catalog.Type, nameVersion string) (repo, tag string, err error) {
	nv := catalog.Normalize(nameVersion)
	versionLess := isVersionLess(typ, nv)
	if err := ValidateIdentity(nv, versionLess); err != nil {
		return "", "", err
	}
	if versionLess {
		return nv, RollingTag, nil
	}
	idx := strings.LastIndex(nv, "/")
	return nv[:idx], nv[idx+1:], nil
}

// projectTagSeparator joins an encoded artifact name to its identity prefix in a
// project's flat tag namespace.
//
// Doubled for the same reason `--` works as the name separator: catalog's
// segment grammar allows a single `.`, `_` or `-` between alphanumerics and
// never two, so a doubled separator cannot occur inside a name. `@`, which the
// store filename uses, is not legal in an OCI tag at all.
const projectTagSeparator = "__"

// ProjectTags renders the tags a project publishes one artifact under, canonical
// first. A project keeps every artifact in one repository, so the name lives in
// the tag rather than in the repository path.
//
// The qualified tag is the retention anchor and every artifact gets one. Nothing
// fetches by tag — a lock records the platform manifest digest — so what a tag
// has to do is keep a manifest referenced: what no tag reaches is unreferenced,
// and unreferenced content is the registry's to reclaim. A plain name tag can
// only keep one identity alive, so re-selecting would orphan the build an older
// lock still points at.
//
// selected additionally writes the plain name tag. It is a human handle and a
// moving pointer, like `latest`: it says which identity this project uses now,
// and dropping it would break no restore.
func ProjectTags(m meta.Manifest, selected bool) ([]string, error) {
	encoded := image.EncodeArtifactName(catalog.Normalize(m.Name))
	if encoded == "" {
		return nil, fmt.Errorf("artifact has no name to publish under")
	}
	sha := strings.TrimPrefix(m.Keys.Identity.SHA256, "sha256:")
	if len(sha) < capsule.IdentityChars {
		return nil, fmt.Errorf("%s records no identity to publish under", m.Name)
	}
	qualified := encoded + projectTagSeparator + sha[:capsule.IdentityChars]

	tags := []string{qualified}
	if selected {
		tags = append(tags, encoded)
	}
	for _, tag := range tags {
		// Refused, never truncated: a truncated tag is a different artifact's
		// address, and the name is what a puller reads back out of it.
		if len(tag) > maxTagLength || !ociTagPattern.MatchString(tag) {
			return nil, fmt.Errorf("%s does not fit an OCI tag as %q", m.Name, tag)
		}
	}
	return tags, nil
}

// ParseProjectTag decodes a project tag back into the artifact name it carries
// and the identity prefix qualifying it, if any.
//
// A catalog tag is a single version segment and never contains `--`, so the two
// namespaces cannot be confused: this reports false for one rather than
// inventing a name from it.
func ParseProjectTag(tag string) (name, sha string, ok bool) {
	encoded := tag
	if base, prefix, found := strings.Cut(tag, projectTagSeparator); found {
		if base == "" || !isHex(prefix) {
			return "", "", false
		}
		encoded, sha = base, prefix
	}
	if !strings.Contains(encoded, "--") {
		return "", "", false
	}
	name = image.DecodeArtifactName(encoded)
	if catalog.Normalize(name) != name || name == "" {
		return "", "", false
	}
	return name, sha, true
}

func isHex(s string) bool {
	if s == "" {
		return false
	}
	return strings.TrimLeft(s, "0123456789abcdef") == ""
}
