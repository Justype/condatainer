package build

import (
	"bufio"
	"bytes"
	"context"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/registry"
)

// bootstrap is a definition's Bootstrap: and From: directives. Agent is what
// Apptainer's docs call the bootstrap agent — docker, oras, library, shub.
type bootstrap struct {
	Agent string
	From  string
}

// parseBootstrap reads the Bootstrap and From headers out of a definition,
// stopping at the first %section as Apptainer does.
func parseBootstrap(def []byte) bootstrap {
	var b bootstrap
	scanner := bufio.NewScanner(bytes.NewReader(def))
	scanner.Buffer(make([]byte, 0, 64*1024), 1024*1024)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		// A From: inside %post is shell text, not a directive.
		if strings.HasPrefix(line, "%") {
			break
		}
		key, value, ok := strings.Cut(line, ":")
		if !ok {
			continue
		}
		switch strings.ToLower(strings.TrimSpace(key)) {
		case "bootstrap":
			b.Agent = strings.ToLower(strings.TrimSpace(value))
		case "from":
			b.From = strings.TrimSpace(value)
		}
	}
	return b
}

// resolvable reports whether this bootstrap pulls from an OCI registry, the only
// kind oras can resolve to a digest.
func (b bootstrap) resolvable() bool {
	switch b.Agent {
	case "docker", "oras":
		return b.From != ""
	}
	return false
}

// remote reports whether this bootstrap names something that can change under a
// fixed reference, and so needs a from line even when it cannot be resolved.
// Everything else has no upstream, or one that is already a local file.
func (b bootstrap) remote() bool {
	switch b.Agent {
	// library and shub have their own protocols: unresolvable here, but moving
	// targets all the same.
	case "docker", "oras", "library", "shub":
		return b.From != ""
	}
	return false
}

// pinned returns From with its tag replaced by digest. A reference that already
// names one is returned as-is.
func (b bootstrap) pinned(digest string) string {
	if digest == "" || !strings.HasPrefix(digest, "sha256:") {
		return b.From
	}
	if strings.Contains(b.From, "@") {
		return b.From
	}
	name := b.From
	// A colon is only a tag separator after the last slash; before it, it is a
	// registry port.
	if i := strings.LastIndex(name, ":"); i > strings.LastIndex(name, "/") {
		name = name[:i]
	}
	return name + "@" + digest
}

// captureSynthesizedRecipe adopts a synthesized definition as the artifact's
// recipe when the build had no source of its own — the scheme:// path, where
// nothing was read from disk. Without it such an image would carry no keys.
//
// The generated header carries a build date and still hashes stably, since
// whole-line comments never reach a preimage.
func (b *BuildObject) captureSynthesizedRecipe(path string, data []byte) {
	if b.spec.Source.BuildType() != 0 {
		return
	}
	b.spec.Source.Definition = &DefinitionSource{
		File: SourceFile{Name: filepath.Base(path), Data: data},
	}
	b.embedSource(SourceFile{Name: meta.RecipeFileName, Data: data})
}

// pinDefinition rewrites a definition's From: line to name digest, and reports
// whether it changed anything.
func pinDefinition(def []byte, digest string) ([]byte, bool) {
	if digest == "" || !strings.HasPrefix(digest, "sha256:") {
		return def, false
	}
	boot := parseBootstrap(def)
	if !boot.resolvable() {
		return def, false
	}
	pinned := boot.pinned(digest)
	if pinned == boot.From {
		return def, false
	}

	lines := strings.Split(string(def), "\n")
	for i, line := range lines {
		if strings.HasPrefix(strings.TrimSpace(line), "%") {
			break
		}
		key, _, ok := strings.Cut(line, ":")
		if !ok || strings.ToLower(strings.TrimSpace(key)) != "from" {
			continue
		}
		lines[i] = "From: " + pinned
		return []byte(strings.Join(lines, "\n")), true
	}
	return def, false
}

// resolveUpstream records what the definition's bootstrap reference points at
// right now, before Apptainer pulls it.
//
// It never fails the build: an unreachable registry is recorded as
// meta.Unrecorded, so the gap is visible rather than silent.
func (b *BuildObject) resolveUpstream(ctx context.Context, def []byte) {
	if b.spec.Source.Definition == nil {
		return
	}
	boot := parseBootstrap(def)
	if !boot.remote() {
		return
	}

	from := &meta.From{Bootstrap: boot.Agent, Ref: boot.From, Digest: meta.Unrecorded}
	if boot.resolvable() {
		digest, err := registry.Resolve(ctx, boot.From)
		if err != nil {
			logging.FromContext(ctx).Warn("could not resolve the upstream image; building without a pinned base",
				"from", boot.From, "err", err)
		} else {
			from.Digest = digest
			logging.FromContext(ctx).Info("resolved upstream image", "from", boot.From, "digest", digest)
		}
	}
	b.spec.Source.Definition.From = from
}
