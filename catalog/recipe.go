package catalog

import (
	"bytes"
	"io"
	"strings"
)

// Entry is what the index records about a recipe: enough to resolve a name and
// walk the dependency graph without fetching anything.
type Entry struct {
	Name           string              `json:"-"` // the index key
	Path           string              `json:"path"`
	Type           Type                `json:"type"`
	Description    string              `json:"description,omitempty"`
	URL            string              `json:"url,omitempty"`
	Deps           []string            `json:"deps,omitempty"` // #DEP: order, as written
	IsTemplate     bool                `json:"is_template,omitempty"`
	TargetTemplate string              `json:"target_template,omitempty"`
	PH             map[string][]string `json:"ph,omitempty"` // value order is load-bearing, key order is not
}

// Recipe is an Entry plus everything only reading the file can tell.
type Recipe struct {
	Entry
	Env        []EnvVar
	Inputs     []string // #INPUT: prompts, in order
	Directives []string // #SBATCH / #PBS / #BSUB, verbatim
	Arch       Arch     // #ARCH:; empty means the default, ArchNative
	// Text is the recipe as fetched, tokens and all. It is what an artifact
	// embeds and what a rebuild starts from.
	Text []byte
	// Rendered is Text with a template's placeholders substituted — what the
	// build runs. Empty for a recipe that is not a template.
	Rendered []byte
}

// Script returns the bytes to execute: the rendered copy for an expanded
// template, the recipe itself otherwise.
func (r *Recipe) Script() []byte {
	if len(r.Rendered) > 0 {
		return r.Rendered
	}
	return r.Text
}

// Arch is what a recipe asserts about where its payload may run. Nothing can
// work this out by inspection — only the author knows whether a build produced
// machine code, a tool-specific binary dump, or bytes that mean the same
// everywhere — so the default is strict and portability must be declared.
type Arch string

const (
	ArchNative Arch = "native" // the default: runs only where it was built
	ArchNoarch Arch = "noarch" // runs anywhere
)

// EnvVar is one #ENV: contribution. Segments keep {prefix} intact, since the
// mount path is not known until the artifact is loaded.
type EnvVar struct {
	Key      string
	Segments []Segment
	Note     string // the ## note
}

// Segment is a run of literal text, or a {token} when IsToken.
type Segment struct {
	Text    string
	IsToken bool
}

// Value renders the segments with tokens substituted from vars. A token with no
// entry is left in braces, which is how {prefix} survives to the manifest.
func (e EnvVar) Value(vars map[string]string) string {
	var b strings.Builder
	for _, s := range e.Segments {
		switch v, ok := vars[s.Text]; {
		case !s.IsToken:
			b.WriteString(s.Text)
		case ok:
			b.WriteString(v)
		default:
			b.WriteByte('{')
			b.WriteString(s.Text)
			b.WriteByte('}')
		}
	}
	return b.String()
}

// ParseRecipe reads a recipe's headers and body. path is relative to the source
// root and carries what the file cannot: the module name, and whether it is a
// .def. A leading recipes/ is optional, so an index path and a bare module path
// both work.
func ParseRecipe(path string, r io.Reader) (*Recipe, error) {
	text, err := io.ReadAll(r)
	if err != nil {
		return nil, err
	}

	isDef := strings.HasSuffix(path, ".def")
	rec := &Recipe{Text: text}
	rec.Path = path
	rec.Name = strings.TrimPrefix(strings.TrimSuffix(path, ".def"), recipesDir+"/")

	var declared string
	for _, line := range headerBlock(text) {
		switch {
		case strings.HasPrefix(line, "#SBATCH"),
			strings.HasPrefix(line, "#PBS"),
			strings.HasPrefix(line, "#BSUB"):
			rec.Directives = append(rec.Directives, line)
			continue
		}

		key, rest, ok := strings.Cut(line, ":")
		if !ok {
			continue
		}
		value, note := splitNote(rest)

		switch key {
		case "#TYPE":
			declared = strings.ToLower(value)
		case "#ARCH":
			if rec.Arch == "" {
				rec.Arch = Arch(strings.ToLower(value))
			}
		case "#DESC":
			rec.Description = firstOf(rec.Description, value)
		case "#URL":
			rec.URL = firstOf(rec.URL, value)
		case "#TARGET":
			rec.TargetTemplate = firstOf(rec.TargetTemplate, value)
		case "#DEP":
			if value != "" {
				rec.Deps = append(rec.Deps, Normalize(value))
			}
		case "#ENV":
			if env, ok := parseEnv(value, note); ok {
				rec.Env = append(rec.Env, env)
			}
		case "#INPUT":
			if prompt := strings.TrimSpace(value); prompt != "" {
				rec.Inputs = append(rec.Inputs, prompt)
			}
		case "#PH":
			name, values, ok := parsePH(value)
			if !ok {
				continue
			}
			if rec.PH == nil {
				rec.PH = map[string][]string{}
			}
			if _, dup := rec.PH[name]; !dup {
				rec.PH[name] = values
			}
		}
	}

	rec.IsTemplate = len(rec.PH) > 0 && rec.TargetTemplate != ""
	rec.Type = DeriveType(rec.Name, rec.TargetTemplate, isDef, declared)
	return rec, nil
}

// headerBoundary returns the byte offset where a recipe's body begins.
//
// Scanning from offset 0, a line is header when its content — after removing the
// line terminator and any leading spaces or tabs — is empty or begins with '#'.
// The body starts at the first byte of the first line that fails that test, or at
// EOF when none does. Two consequences are worth stating: a leading #!/bin/bash
// is header, and so is a blank line before the first command.
//
// It bounds what ParseRecipe reads, which is what keeps a #DEP: in a heredoc — or
// a #SBATCH in a file the recipe writes — inert. Keys do not use it at all: both
// hash StripComments, which removes the header along with every other comment.
//
// Nothing here parses shell. Shell code carries '#' inside strings and ${x#y}
// expansions, so a comment-aware rule could move the boundary on a body that
// never changed.
func headerBoundary(text []byte) int {
	off := 0
	for off < len(text) {
		next := len(text)
		line := text[off:]
		if end := bytes.IndexByte(line, '\n'); end >= 0 {
			line, next = line[:end], off+end+1
		}
		trimmed := bytes.TrimLeft(bytes.TrimSuffix(line, []byte("\r")), " \t")
		if len(trimmed) > 0 && trimmed[0] != '#' {
			return off
		}
		off = next
	}
	return len(text)
}

// StripComments removes every whole-line comment from a recipe, keeping blank
// lines and every other byte as it stands. It is the recipe preimage both keys
// hash: what the recipe *does*, with everything that only describes it gone.
//
// The header is all comments, so it disappears with the rest — which is the
// point. Every header item that reaches a key already has its own record line
// (#TYPE: to type=, #ENV: to env=, #DEP: to dep=, a selected #PH: to ph=), and
// hashing the file too would re-admit the ones deliberately left out: a #PH:
// menu that autoupdate grows, a reworded #DESC:, a changed scheduler directive.
// None of those change what a build produces.
//
// What the two keys do not share is dependencies, not text: identity carries
// every direct one, equivalence only the projected subset.
//
// Whole-line only. A trailing comment cannot be removed without parsing shell,
// since `echo "a # b"` and `${v#pre}` both carry a '#' that is not one. Blank
// lines are kept, so a blank line added to the header does move both keys.
func StripComments(text []byte) []byte {
	out := make([]byte, 0, len(text))
	for off := 0; off < len(text); {
		next := len(text)
		line := text[off:]
		if end := bytes.IndexByte(line, '\n'); end >= 0 {
			line, next = line[:end+1], off+end+1
		}
		content := bytes.TrimSuffix(bytes.TrimSuffix(line, []byte("\n")), []byte("\r"))
		if trimmed := bytes.TrimLeft(content, " \t"); len(trimmed) == 0 || trimmed[0] != '#' {
			out = append(out, line...)
		}
		off = next
	}
	return out
}

// headerBlock returns the header's comment lines, trimmed and without its blank
// ones. It reads only up to headerBoundary, so a #DEP: in a heredoc — or a
// #SBATCH in a job wrapper the recipe writes — stays inert.
func headerBlock(text []byte) []string {
	var out []string
	for line := range strings.SplitSeq(string(text[:headerBoundary(text)]), "\n") {
		if line = strings.TrimSpace(line); line != "" {
			out = append(out, line)
		}
	}
	return out
}

// splitNote separates a header value from its ## note.
func splitNote(s string) (value, note string) {
	if before, after, ok := strings.Cut(s, "##"); ok {
		return strings.TrimSpace(before), strings.TrimSpace(after)
	}
	return strings.TrimSpace(s), ""
}

// parseEnv splits KEY=value and tokenizes the value on {name}.
func parseEnv(value, note string) (EnvVar, bool) {
	key, raw, ok := strings.Cut(value, "=")
	key = strings.TrimSpace(key)
	if !ok || key == "" {
		return EnvVar{}, false
	}
	return EnvVar{Key: key, Segments: tokenize(raw), Note: note}, true
}

// tokenize splits text into literal runs and {token} runs.
func tokenize(text string) []Segment {
	var out []Segment
	for text != "" {
		open := strings.IndexByte(text, '{')
		if open < 0 {
			return append(out, Segment{Text: text})
		}
		close := strings.IndexByte(text[open:], '}')
		if close < 0 {
			return append(out, Segment{Text: text})
		}
		close += open
		if open > 0 {
			out = append(out, Segment{Text: text[:open]})
		}
		out = append(out, Segment{Text: text[open+1 : close], IsToken: true})
		text = text[close+1:]
	}
	return out
}

// firstOf keeps the first value a repeated header set.
func firstOf(have, next string) string {
	if have != "" {
		return have
	}
	return next
}
