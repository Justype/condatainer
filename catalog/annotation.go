package catalog

import "strings"

// Annotation is one `#KEY: value ## note` declaration.
//
// Every consumer of script metadata reads these — recipes, user scripts,
// project scanning — so a key is tokenized in exactly one place and means the
// same thing everywhere it appears.
type Annotation struct {
	// Key is the notation including its '#', e.g. "#DEP".
	Key string
	// Value is the text after the colon, with any trailing `# comment` removed.
	Value string
	// Note is the text after `##`, unparsed. Empty when there is none.
	Note string
	// Line is the 1-based line the annotation was written on.
	Line int
}

// ScanAnnotations returns every annotation in text, in the order written.
//
// A line qualifies when, after leading blanks and tabs, it begins with '#',
// then a key of upper-case letters, digits or underscores, then ':'. Position
// carries no meaning: a declaration is wherever its author put it, so a header
// block is a convention rather than a rule the parser enforces.
//
// The requirement that a key be upper-case is what keeps ordinary prose out.
// A comment like `# note: rerun this weekly` has a lower-case key and is not an
// annotation; `#DEP: star/2.7.11b` is.
func ScanAnnotations(text []byte) []Annotation {
	var out []Annotation
	for i, raw := range strings.Split(string(text), "\n") {
		line := strings.TrimLeft(strings.TrimSuffix(raw, "\r"), " \t")
		if !strings.HasPrefix(line, "#") {
			continue
		}
		key, rest, ok := cutKey(line)
		if !ok {
			continue
		}
		value, note := splitNote(rest)
		out = append(out, Annotation{Key: key, Value: stripInlineComment(value), Note: note, Line: i + 1})
	}
	return out
}

// Find returns the value of the first annotation with key, and whether one was
// present. Most keys are single-valued, and the first wins.
func Find(annotations []Annotation, key string) (string, bool) {
	for _, annotation := range annotations {
		if annotation.Key == key {
			return annotation.Value, true
		}
	}
	return "", false
}

// Select returns every annotation with key, in order.
func Select(annotations []Annotation, key string) []Annotation {
	var out []Annotation
	for _, annotation := range annotations {
		if annotation.Key == key {
			out = append(out, annotation)
		}
	}
	return out
}

// cutKey splits `#KEY:rest`, requiring an upper-case key immediately after the
// '#'. It reports false for a line that only looks like one, such as a shebang
// or a prose comment containing a colon.
func cutKey(line string) (key, rest string, ok bool) {
	body, rest, ok := strings.Cut(line[1:], ":")
	if !ok || body == "" {
		return "", "", false
	}
	for _, r := range body {
		if (r < 'A' || r > 'Z') && (r < '0' || r > '9') && r != '_' {
			return "", "", false
		}
	}
	return "#" + body, rest, true
}

// stripInlineComment removes a trailing `# comment` from an annotation value.
//
// utils.StripInlineComment is the same three lines, kept separate because
// catalog sits below utils and the schedulers reach it through utils. Neither
// belongs in the other's package for the sake of one string operation.
func stripInlineComment(value string) string {
	if i := strings.IndexByte(value, '#'); i >= 0 {
		value = value[:i]
	}
	return strings.TrimSpace(value)
}
