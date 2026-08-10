package conda

import (
	"fmt"
	"slices"
	"sort"
	"strings"

	"go.yaml.in/yaml/v3"
)

// The two exports an image embeds, captured from the environment that was
// actually installed rather than from a second solve.
const (
	ExplicitFileName    = "explicit.txt"
	EnvironmentFileName = "environment.yml"
)

// CanonicalExplicit re-emits `micromamba env export --explicit --no-md5` output:
// the @EXPLICIT marker, then one bare package URL per line, sorted, LF-terminated.
//
// Re-emitting rather than storing the tool's bytes is the point. A Micromamba
// upgrade that reorders or re-spaces its output would otherwise move the key of
// every artifact built after it, and every comparison across that boundary would
// report different with an empty diff. Sorting is safe — Conda does not depend on
// explicit-file order — so the result is still a working `micromamba install
// --file` input, and `sha256sum explicit.txt` reproduces the identity by hand.
//
// A trailing #<md5> is dropped if the exporter emits one. It checksums the
// download rather than what was installed, and the URL already pins channel,
// subdirectory, name, version and build string — everything that tells one
// package from another.
func CanonicalExplicit(raw []byte) ([]byte, error) {
	var (
		urls   []string
		marked bool
	)
	for _, line := range strings.Split(string(raw), "\n") {
		line = strings.TrimSpace(line)
		if line == "" {
			continue
		}
		if !marked {
			// Everything before the marker is generated commentary.
			marked = line == "@EXPLICIT"
			continue
		}
		if strings.HasPrefix(line, "#") {
			continue
		}
		if i := strings.IndexByte(line, '#'); i >= 0 {
			line = strings.TrimSpace(line[:i])
		}
		if line != "" {
			urls = append(urls, line)
		}
	}
	if !marked {
		return nil, fmt.Errorf("conda: export has no @EXPLICIT marker; not an explicit spec")
	}

	sort.Strings(urls)
	var sb strings.Builder
	sb.WriteString("@EXPLICIT\n")
	for _, u := range urls {
		sb.WriteString(u)
		sb.WriteByte('\n')
	}
	return []byte(sb.String()), nil
}

// environmentExport is the subset of `micromamba env export --no-builds` output
// that describes packages rather than where they happened to be installed.
type environmentExport struct {
	Channels     []string `yaml:"channels"`
	Dependencies []any    `yaml:"dependencies"`
}

// CanonicalEnvironment re-emits `micromamba env export --no-builds` output as a
// deterministic environment.yml: channels in priority order, dependencies sorted,
// no name and no prefix.
//
// name and prefix describe a temporary local environment rather than its
// packages, so they are dropped. Channels are kept because the file has to stay
// usable — `micromamba create -f` without them resolves against whatever the next
// machine has configured, which is not a rebuild.
//
// The channel *set* is the export's, never the configuration's: those are the
// channels that actually provided packages. priority supplies only the order,
// because Micromamba alphabetizes on export and its sequence is not the solve
// order. A channel absent from priority is appended after the known ones, sorted,
// so the result is total whatever the input. Adding an unused channel to a site's
// config therefore changes no artifact's key.
//
// A pip sub-list is preserved and sorted as exported. Nothing chases it further:
// a centrally managed overlay does not carry pip installs.
func CanonicalEnvironment(raw []byte, priority []string) ([]byte, error) {
	var doc environmentExport
	if err := yaml.Unmarshal(raw, &doc); err != nil {
		return nil, fmt.Errorf("conda: cannot parse environment export: %w", err)
	}

	var deps, pip []string
	for _, entry := range doc.Dependencies {
		switch v := entry.(type) {
		case string:
			deps = append(deps, v)
		case map[string]any:
			list, ok := v["pip"].([]any)
			if !ok {
				continue
			}
			for _, item := range list {
				pip = append(pip, fmt.Sprint(item))
			}
		}
	}
	sort.Strings(deps)
	sort.Strings(pip)

	var sb strings.Builder
	sb.WriteString("channels:\n")
	for _, ch := range OrderChannels(doc.Channels, priority) {
		fmt.Fprintf(&sb, "  - %s\n", ch)
	}
	sb.WriteString("dependencies:\n")
	for _, d := range deps {
		fmt.Fprintf(&sb, "  - %s\n", d)
	}
	if len(pip) > 0 {
		sb.WriteString("  - pip:\n")
		for _, p := range pip {
			fmt.Fprintf(&sb, "      - %s\n", p)
		}
	}
	return []byte(sb.String()), nil
}

// OrderChannels returns the exported channels in the configured priority order.
//
// Only channels present in exported survive: those are the ones that provided a
// package. Anything configured but unused is left out, and anything used but
// unconfigured — a mirror, or a `bioconda::star` annotation — is appended after
// the known ones, sorted, rather than dropped.
func OrderChannels(exported, priority []string) []string {
	present := make(map[string]bool, len(exported))
	for _, ch := range exported {
		if ch = strings.TrimSpace(ch); ch != "" {
			present[ch] = true
		}
	}

	out := make([]string, 0, len(present))
	for _, ch := range priority {
		if present[ch] {
			out = append(out, ch)
			delete(present, ch)
		}
	}
	rest := make([]string, 0, len(present))
	for ch := range present {
		rest = append(rest, ch)
	}
	slices.Sort(rest)
	return append(out, rest...)
}
