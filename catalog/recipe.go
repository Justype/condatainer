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
	Sources    []SourceURL // #SOURCE: declarations, in order
	Inputs     []string    // #INPUT: prompts, in order
	Directives []string    // #SBATCH / #PBS / #BSUB, verbatim
	Arch       Arch        // #ARCH:; empty means the default, ArchNative
	// License is #LICENSE: as written — an SPDX expression, kept verbatim and
	// never parsed. Deriving redistribution permission from a licence
	// expression is a judgement a tool gets wrong in the permissive direction,
	// and that direction cannot be taken back, so this documents and decides
	// nothing. Redistribute is the decision.
	License string
	// Redistribute is #REDISTRIBUTE: lower-cased, "" when the recipe did not
	// answer. Validate rejects anything but yes or no; Redistributable reads it
	// as the tri-state the answer actually is.
	Redistribute string
	// DeclaredType is #TYPE: lower-cased and verbatim, "" when the recipe did
	// not declare one. Type holds what DeriveType made of it, which for anything
	// but app or data is the path-based default — so this is kept for Validate
	// to reject the declaration rather than let it read as silence.
	DeclaredType string
	// Text is the recipe as fetched, tokens and all. It is what an artifact
	// embeds and what a rebuild starts from.
	Text []byte
	// Rendered is Text with a template's placeholders substituted — what the
	// build runs. Empty for a recipe that is not a template.
	Rendered []byte
}

// Redistributable reports the recipe's #REDISTRIBUTE: answer, or nil when it
// did not answer.
//
// Three states, not two: an absent declaration means the question was never put
// to the author, which is a different fact from a considered "no" and defaults
// differently depending on the artifact's type. A caller that flattens nil to
// false refuses everything nobody has annotated yet.
func (r *Recipe) Redistributable() *bool {
	switch r.Redistribute {
	case "yes":
		yes := true
		return &yes
	case "no":
		no := false
		return &no
	}
	return nil
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
	rec.Directives = scanDirectives(text)
	for _, annotation := range ScanAnnotations(text) {
		value, note := annotation.Value, annotation.Note

		switch annotation.Key {
		case "#TYPE":
			declared = strings.ToLower(value)
			rec.DeclaredType = declared
		case "#ARCH":
			if rec.Arch == "" {
				rec.Arch = Arch(strings.ToLower(value))
			}
		case "#DESC":
			rec.Description = firstOf(rec.Description, value)
		case "#URL":
			rec.URL = firstOf(rec.URL, value)
		case "#LICENSE":
			rec.License = firstOf(rec.License, value)
		case "#REDISTRIBUTE":
			rec.Redistribute = firstOf(rec.Redistribute, strings.ToLower(strings.TrimSpace(value)))
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
		case "#SOURCE":
			if src, ok := parseSource(value); ok {
				rec.Sources = append(rec.Sources, src)
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

// StripComments removes every whole-line comment from a recipe, keeping blank
// lines and every other byte as it stands. It is the recipe preimage both keys
// hash: what the recipe *does*, with everything that only describes it gone.
//
// Annotations are comments, so they disappear with the rest — which is the
// point. Every annotation that reaches a key already has its own record line
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
// lines are kept, so adding one does move both keys.
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

// scanDirectives returns every scheduler directive line, verbatim.
//
// These are not annotations: `#SBATCH --time=01:00:00` has no `#KEY: value`
// shape, and cutting it on the first colon would split the value. They are
// matched by prefix and handed to the scheduler packages to parse.
func scanDirectives(text []byte) []string {
	var out []string
	for line := range strings.SplitSeq(string(text), "\n") {
		line = strings.TrimSpace(strings.TrimSuffix(line, "\r"))
		switch {
		case strings.HasPrefix(line, "#SBATCH"),
			strings.HasPrefix(line, "#PBS"),
			strings.HasPrefix(line, "#BSUB"):
			out = append(out, line)
		}
	}
	return out
}

// splitNote separates an annotation value from its ## note.
func splitNote(s string) (value, note string) {
	if before, after, ok := strings.Cut(s, "##"); ok {
		return strings.TrimSpace(before), strings.TrimSpace(after)
	}
	return strings.TrimSpace(s), ""
}

// SourceURL is one #SOURCE: declaration: a name the recipe body refers to as
// $CNT_SRC_<name>, and where the file comes from.
//
// A link only the user holds — expiring, per-user, EULA-gated — is declared with
// a Prompt instead of a URL, and the tool asks for it. The answer is never
// recorded.
type SourceURL struct {
	Name   string
	URL    string
	Prompt string
}

// askPrefix marks a #SOURCE: whose link is asked for rather than written down.
const askPrefix = "ask:"

// parseSource splits "name url" or "name ask:prompt".
//
// The name becomes an environment variable, so it is restricted to what one can
// hold. A malformed declaration is skipped here and reported by Validate, which
// is where an author is told rather than left with an unset variable.
func parseSource(value string) (SourceURL, bool) {
	name, rest, _ := strings.Cut(strings.TrimSpace(value), " ")
	rest = strings.TrimSpace(rest)
	if !validSourceName(name) || rest == "" {
		return SourceURL{}, false
	}
	if prompt, ok := strings.CutPrefix(rest, askPrefix); ok {
		prompt = strings.TrimSpace(prompt)
		if prompt == "" {
			return SourceURL{}, false
		}
		return SourceURL{Name: name, Prompt: prompt}, true
	}
	if strings.ContainsAny(rest, " \t") {
		return SourceURL{}, false
	}
	return SourceURL{Name: name, URL: rest}, true
}

// Prompts returns every question a build must put to the user, in the order the
// answers are supplied: the #INPUT: prompts, then each #SOURCE: ask: in
// declaration order. #INPUT: answers reach the recipe; the rest are consumed by
// the fetch and never reach it.
func (r *Recipe) Prompts() []string {
	out := append([]string(nil), r.Inputs...)
	for _, src := range r.Sources {
		if src.Prompt != "" {
			out = append(out, src.Prompt)
		}
	}
	return out
}

// scanMalformedSources returns the #SOURCE: values parseSource refused, as
// written. Parsing skips them so one bad line cannot take a collection out of a
// listing; Validate reads them back here to tell the author.
func scanMalformedSources(text []byte) []string {
	var out []string
	for _, annotation := range ScanAnnotations(text) {
		if annotation.Key != "#SOURCE" {
			continue
		}
		if _, ok := parseSource(annotation.Value); !ok {
			out = append(out, annotation.Value)
		}
	}
	return out
}

// validSourceName reports whether name can be the tail of $CNT_SRC_<name>.
func validSourceName(name string) bool {
	if name == "" {
		return false
	}
	for _, r := range name {
		switch {
		case r >= 'a' && r <= 'z', r >= 'A' && r <= 'Z', r >= '0' && r <= '9', r == '_':
		default:
			return false
		}
	}
	return name[0] < '0' || name[0] > '9'
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
