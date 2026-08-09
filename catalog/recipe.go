package catalog

import (
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
	Text       []byte
}

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

// headerBlock returns the leading run of comment and blank lines, trimmed.
// It stops at the first line that is neither, so a #DEP: in a heredoc — or a
// #SBATCH in a job wrapper the recipe writes — stays inert.
func headerBlock(text []byte) []string {
	var out []string
	for line := range strings.SplitSeq(string(text), "\n") {
		line = strings.TrimSpace(line)
		if line == "" {
			continue
		}
		if !strings.HasPrefix(line, "#") {
			break
		}
		out = append(out, line)
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
