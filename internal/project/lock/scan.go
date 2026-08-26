package lock

import (
	"fmt"
	"io/fs"
	"os"
	"path"
	"path/filepath"
	"sort"
	"strings"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/utils"
)

// OverlaysDirName is the conventional directory a project keeps its own
// overlays in, and the one directory of project source that scanning skips.
//
// The scripts there are the recipes that built those overlays — `create -f
// overlays/tool.sh` writes overlays/tool.sqf beside it — and a recipe's #DEP:
// are its overlay's build dependencies, already recorded in that overlay's
// provenance. Reading them would make the project pin what it never mounts.
//
// A directory is the rule rather than anything read out of the tree, because
// only a stated convention answers the same every time: a sibling .sqf appears
// when the overlay is built, and a script's own annotations say nothing about
// whether it is a recipe. Declaring `#DEP: overlays/tool.sqf` from a project
// script is unaffected — that declaration lives in the script, not here.
const OverlaysDirName = "overlays"

// UnpinnedMarker is the note that declares a dependency cannot be pinned:
// `#DEP: env.img  ## unpinned — reason`. The reason is free text this package
// never interprets; only its presence is checked, and only the marker is
// required.
const UnpinnedMarker = "unpinned"

// Kind is what a declaration addresses.
type Kind string

const (
	// KindName addresses an artifact by name, optionally constrained.
	KindName Kind = "name"
	// KindPath addresses an immutable overlay by project-relative path.
	KindPath Kind = "path"
	// KindWritable addresses a writable .img, which has no identity to pin.
	KindWritable Kind = "writable"
	// KindExternal addresses a .sqf outside the project — absolute, or reached
	// by climbing out of the root.
	KindExternal Kind = "external"
)

// Pinnable reports whether an artifact can be pinned for this kind.
//
// A name is pinnable because the store resolves it by identity, and a
// project-relative .sqf because restore owns that path and writes the artifact
// there. The other two are not, for the same underlying reason: there is
// nowhere for restore to put the answer. A writable .img has no identity to
// pin, and an external .sqf is someone else's file — restore does not own that
// path and must never write to it, so pinning one would record a promise it
// could not keep.
func (k Kind) Pinnable() bool { return k == KindName || k == KindPath }

// Request is one declaration a project makes, merged across every script that
// makes it.
type Request struct {
	// Key is the canonical pin key: catalog.Dep.String() for a name, and
	// PathPrefix plus the cleaned relative path for a path.
	Key  string
	Kind Kind
	// Dep is the parsed declaration, set for KindName only.
	Dep catalog.Dep
	// Path is the cleaned project-relative path, set for the path kinds.
	Path string
	// Unpinned reports the `## unpinned` marker, and Reason is whatever
	// followed it, unparsed and possibly empty.
	Unpinned bool
	Reason   string
	// Scripts are the project-relative scripts that declared this, sorted.
	// Displayed, never serialized: rescanning finds them again.
	Scripts []string
	// first locates the declaration this request was created from, so a problem
	// only decidable after every script has been read still points at a line.
	first Finding
}

// ConstraintReason reports why a version constraint cannot appear in a project
// declaration, or "".
//
// A range is a build-recipe feature: it lets a build reuse a satisfying version
// that is already installed instead of producing another, and for reference
// data the tool version frequently does not matter — `samtools faidx` writes the
// same .fai whichever recent samtools ran. An analysis reports results from what
// it mounts, where "any version in this range" is not a claim worth attaching to
// a result, and a lock exists to say which one.
//
// It is also what keeps one module from having two keys. The key is
// NameVersion()+Op+Min, so `star/2.7.11b` and `star/2.7.11b>=2.7.0` would be
// separate pins that may point at different artifacts, with nothing in
// the lock saying which a given script meant.
func ConstraintReason(dep catalog.Dep, declaration string) string {
	if dep.Op == "" {
		return ""
	}
	return fmt.Sprintf("%s carries a version constraint, which only a build recipe may use; declare %s exactly",
		declaration, dep.NameVersion())
}

// undeclaredReason explains why an unpinnable declaration needs the marker.
// Kind-specific, because the two are unpinnable for different reasons and the
// reader can only act on the one that applies.
func undeclaredReason(kind Kind, target string) string {
	switch kind {
	case KindWritable:
		return fmt.Sprintf("%s is writable, so it has no identity to pin; declare it `## %s`", target, UnpinnedMarker)
	default:
		// A `../` declaration is the case worth spelling out: it looks
		// script-relative, and under one anchor for the whole project it is not.
		anchor := ""
		if !path.IsAbs(target) {
			anchor = " (a path in a project is relative to the project root, not to the script declaring it)"
		}
		return fmt.Sprintf("%s is outside the project%s, so restore cannot own that path; declare it `## %s`",
			target, anchor, UnpinnedMarker)
	}
}

// Finding is a declaration problem the scanner can see without leaving the
// checkout. It is reported rather than raised so the caller decides severity:
// `project lock` lists them, `project validate` fails on them.
type Finding struct {
	Script string
	Line   int
	Text   string
	Reason string
}

// ScanResult is what a project declares.
type ScanResult struct {
	// Requests are unique declarations, sorted by key.
	Requests []Request
	// Scripts are every scanned file, project-relative and sorted.
	Scripts  []string
	Findings []Finding
}

// ScanScript reads one script's declarations, merged and finalized exactly as
// Scan does for a whole project.
//
// Execution reads a script through this rather than through its own parse, so a
// run and a lock cannot disagree about what the script declares — which is the
// only reason a lock means anything at run time. The script must be inside root.
func ScanScript(root, script string) (*ScanResult, error) {
	root, err := filepath.Abs(root)
	if err != nil {
		return nil, err
	}
	script, err = filepath.Abs(script)
	if err != nil {
		return nil, err
	}
	rel, err := filepath.Rel(root, script)
	if err != nil || rel == ".." || strings.HasPrefix(rel, ".."+string(filepath.Separator)) {
		return nil, fmt.Errorf("%s is not inside %s", script, root)
	}
	result := &ScanResult{Scripts: []string{filepath.ToSlash(rel)}}
	merged := map[string]*Request{}
	if err := scanScript(script, filepath.ToSlash(rel), merged, result); err != nil {
		return nil, err
	}
	finalize(merged, result)
	return result, nil
}

// ScanOptions tunes discovery.
type ScanOptions struct {
	// ExcludeDirs are additional directory names to skip anywhere in the tree.
	// cnt-lock, overlays and every dot-directory are always skipped.
	ExcludeDirs []string
}

// Scan walks a project root and returns what its scripts declare.
//
// A declaration counts wherever it is written; position carries no meaning, so
// the scanner and the runtime read a script the same way.
func Scan(root string, opts ScanOptions) (*ScanResult, error) {
	skip := map[string]bool{DirName: true, OverlaysDirName: true}
	for _, dir := range opts.ExcludeDirs {
		if dir = strings.TrimSpace(dir); dir != "" {
			skip[dir] = true
		}
	}

	merged := map[string]*Request{}
	result := &ScanResult{}

	walk := func(path string, entry fs.DirEntry, err error) error {
		if err != nil {
			return err
		}
		rel, relErr := filepath.Rel(root, path)
		if relErr != nil {
			return relErr
		}
		if entry.IsDir() {
			// A dot-directory is tool state, not project source: .git, .venv,
			// .snakemake, .tox and a local conda env all carry shell scripts of
			// their own, and a #DEP: in one of those is not this project's
			// declaration. The root itself is read even when it is hidden.
			// overlays/ is skipped for a different reason — see OverlaysDirName.
			if path != root && (skip[entry.Name()] || strings.HasPrefix(entry.Name(), ".")) {
				return filepath.SkipDir
			}
			return nil
		}
		// A directory symlink is never followed, and a symlinked script is not
		// read: either can point outside the checkout, and a lock describes the
		// checkout.
		if entry.Type()&fs.ModeSymlink != 0 || !entry.Type().IsRegular() {
			return nil
		}
		if !isShellScript(entry.Name()) {
			return nil
		}
		relSlash := filepath.ToSlash(rel)
		result.Scripts = append(result.Scripts, relSlash)
		return scanScript(path, relSlash, merged, result)
	}

	if err := filepath.WalkDir(root, walk); err != nil {
		return nil, err
	}

	finalize(merged, result)
	return result, nil
}

// finalize turns merged declarations into a sorted result and reports any
// unpinnable one left undeclared.
//
// The marker merges across scripts, so whether one is missing is decidable only
// once every script in scope has been read — which is why this runs here rather
// than as each declaration is parsed.
func finalize(merged map[string]*Request, result *ScanResult) {
	for _, request := range merged {
		sort.Strings(request.Scripts)
		if !request.Kind.Pinnable() && !request.Unpinned {
			finding := request.first
			finding.Reason = undeclaredReason(request.Kind, request.Path)
			result.Findings = append(result.Findings, finding)
		}
		result.Requests = append(result.Requests, *request)
	}
	sort.Slice(result.Requests, func(i, j int) bool { return result.Requests[i].Key < result.Requests[j].Key })
	sort.Strings(result.Scripts)
	sort.Slice(result.Findings, func(i, j int) bool {
		if result.Findings[i].Script != result.Findings[j].Script {
			return result.Findings[i].Script < result.Findings[j].Script
		}
		return result.Findings[i].Line < result.Findings[j].Line
	})
}

// isShellScript reports whether a file should be read for declarations.
//
// Extension only. `run` executes every project script with /bin/bash, so no
// other shell's script could run here anyway, and sniffing a shebang to catch
// extensionless files cost more than it bought: the test was a substring match
// that fired on any interpreter path containing "sh", which on an HPC filesystem
// means /home/shared/... and every user called josh.
func isShellScript(name string) bool {
	switch strings.ToLower(filepath.Ext(name)) {
	case ".sh", ".bash":
		return true
	}
	return false
}

// scanScript reads one script's declarations.
func scanScript(path, rel string, merged map[string]*Request, result *ScanResult) error {
	text, err := os.ReadFile(path)
	if err != nil {
		return err
	}

	for _, annotation := range catalog.Select(catalog.ScanAnnotations(text), "#DEP") {
		if annotation.Value == "" {
			continue
		}
		request, reason := ParseDeclaration(annotation.Value, annotation.Note)
		if reason != "" {
			result.Findings = append(result.Findings, Finding{
				Script: rel, Line: annotation.Line, Text: annotation.Value, Reason: reason})
			continue
		}
		have, ok := merged[request.Key]
		if !ok {
			request.Scripts = []string{rel}
			request.first = Finding{Script: rel, Line: annotation.Line, Text: annotation.Value}
			merged[request.Key] = &request
			continue
		}
		if have.Scripts[len(have.Scripts)-1] != rel {
			have.Scripts = append(have.Scripts, rel)
		}
		// One script marking a dependency unpinned marks it for the project:
		// the claim is about the artifact, not about the file.
		if request.Unpinned && !have.Unpinned {
			have.Unpinned, have.Reason = true, request.Reason
		}
	}
	return nil
}

// ParseDeclaration turns one declaration into a Request, or returns why it
// cannot be one. The second return is empty exactly when the first is usable.
//
// This is the grammar a `#DEP:` uses and the grammar `exec -o` accepts, which is
// why it is exported: a name typed on the command line inside a project has to
// classify the same way the scanner classifies the same text in a script, or the
// two would disagree about what the lock covers. note is the `##` comment, empty
// for a command-line argument.
func ParseDeclaration(value, note string) (Request, string) {
	unpinned, reason := parseNote(note)

	if utils.IsOverlay(value) {
		clean := filepath.ToSlash(filepath.Clean(value))
		kind := KindPath
		switch {
		case utils.IsImg(value):
			kind = KindWritable
		case path.IsAbs(clean), clean == "..", strings.HasPrefix(clean, "../"):
			// Outside the project, so restore has no path it may write.
			kind = KindExternal
		}
		return Request{
			Key: PathPrefix + clean, Kind: kind, Path: clean,
			Unpinned: unpinned, Reason: reason,
		}, ""
	}

	dep, err := catalog.ParseDep(value)
	if err != nil {
		return Request{}, "not a usable dependency: " + err.Error()
	}
	if why := ConstraintReason(dep, value); why != "" {
		return Request{}, why
	}
	return Request{
		Key: dep.String(), Kind: KindName, Dep: dep,
		Unpinned: unpinned, Reason: reason,
	}, ""
}

// parseRequest turns a canonical pin key back into the dependency it
// renders. It is Request.Key's inverse for a named request, and the one place
// that knows a key is not simply a name.
func parseRequest(key string) (catalog.Dep, error) {
	if path, ok := strings.CutPrefix(key, PathPrefix); ok {
		return catalog.Dep{}, fmt.Errorf("%q addresses a path, not a name", path)
	}
	return catalog.ParseDep(key)
}

// parseNote reads the `## unpinned` marker. The marker must come first; the
// rest is a reason nothing interprets. A note that is not the marker is just a
// comment and says nothing about pinning.
func parseNote(note string) (unpinned bool, reason string) {
	if note == "" {
		return false, ""
	}
	first, rest, _ := strings.Cut(note, " ")
	if !strings.EqualFold(strings.TrimSpace(first), UnpinnedMarker) {
		return false, ""
	}
	return true, strings.TrimSpace(strings.TrimLeft(strings.TrimSpace(rest), "-—:"))
}
