package publish

import (
	"context"
	"os/exec"
	"strings"

	"github.com/Justype/condatainer/internal/logging"
)

// OriginURL reports the project's code repository as an https URL, derived from
// the checkout's `origin` remote, or "" when it cannot say.
//
// Only origin: which of several remotes is "the" project is not something a tool
// should decide silently.
//
// Empty is an ordinary answer, never an error. No git, no checkout, no origin,
// or a remote naming no browsable repository all mean the same thing to the
// caller — the annotation is omitted, and --source supplies one.
func OriginURL(ctx context.Context, root string) string {
	cmd := exec.CommandContext(ctx, "git", "-C", root, "remote", "get-url", "origin")
	out, err := cmd.Output()
	if err != nil {
		logging.FromContext(ctx).Debug("no git origin to derive the project source from", "root", root, "err", err)
		return ""
	}
	return NormalizeRemoteURL(string(out))
}

// NormalizeRemoteURL folds a git remote into the one https URL an annotation
// should carry, and returns "" for a remote that names no browsable repository.
//
// It accepts the https, ssh:// and scp-style forms, since `host:owner/repo` is
// the layout web git hosts share, and trims `.git` and a trailing slash.
//
// The spellings are not cosmetic variants: `git@github.com:o/r.git` and
// `https://github.com/o/r` are what SSH and HTTPS clones leave behind, so two
// collaborators on one project routinely have different ones. Publishing
// whichever the pusher happened to have would put two annotations on one
// project's packages, and GHCR links a package to a repository only on an exact
// match with the https form.
//
// A path that is absolute (`git@server:/srv/git/x.git`) or has one component
// names no owner and no web UI, so it returns "" rather than a URL nothing
// serves.
func NormalizeRemoteURL(raw string) string {
	url := strings.TrimSpace(raw)
	switch {
	case strings.HasPrefix(url, "https://"):
		url = strings.TrimPrefix(url, "https://")
	case strings.HasPrefix(url, "http://"):
		url = strings.TrimPrefix(url, "http://")
	case strings.HasPrefix(url, "ssh://git@"):
		url = strings.TrimPrefix(url, "ssh://git@")
	case strings.HasPrefix(url, "git@"):
		// scp-style: the colon separates host from path and is not a port.
		url = strings.Replace(strings.TrimPrefix(url, "git@"), ":", "/", 1)
	default:
		return ""
	}
	url = strings.TrimSuffix(strings.Trim(url, "/"), ".git")

	host, repoPath, found := strings.Cut(url, "/")
	if !found || host == "" || repoPath == "" || strings.HasPrefix(repoPath, "/") {
		return ""
	}
	// At least owner/repo. Nested groups make a deeper path a repository on most
	// hosts; on github.com it is a page inside one, which GHCR links to nothing.
	segments := strings.Split(repoPath, "/")
	if len(segments) < 2 || (host == "github.com" && len(segments) != 2) {
		return ""
	}
	for _, segment := range segments {
		if segment == "" || segment == "." || segment == ".." {
			return ""
		}
	}
	return "https://" + host + "/" + repoPath
}
