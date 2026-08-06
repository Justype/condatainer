package catalog

import (
	"fmt"
	"maps"
	"regexp"
	"slices"
	"strconv"
	"strings"
)

// Template is a #TARGET: pattern: the module path a recipe builds, with {name}
// placeholders standing in for #PH: values.
type Template struct {
	Raw   string
	parts []Segment
}

// Variant is one concrete member of a template's family.
type Variant struct {
	Name string
	Vars map[string]string
}

// NewTemplate parses a #TARGET: pattern.
func NewTemplate(raw string) *Template {
	return &Template{Raw: raw, parts: tokenize(raw)}
}

// Names lists the placeholders in the order they appear, which is the order the
// built name reads in. Nothing stores a separate key order.
func (t *Template) Names() []string {
	var out []string
	for _, p := range t.parts {
		if p.IsToken && !slices.Contains(out, p.Text) {
			out = append(out, p.Text)
		}
	}
	return out
}

// Fill substitutes vars into the pattern. Every placeholder must have a value:
// a half-filled name would name nothing.
func (t *Template) Fill(vars map[string]string) (string, error) {
	var b strings.Builder
	for _, p := range t.parts {
		if !p.IsToken {
			b.WriteString(p.Text)
			continue
		}
		v, ok := vars[p.Text]
		if !ok {
			return "", fmt.Errorf("catalog: %s has no value for {%s}", t.Raw, p.Text)
		}
		b.WriteString(v)
	}
	return b.String(), nil
}

// Match recovers the values that would produce name, matching against the #PH:
// value sets rather than a wildcard: a closed set is its own pattern, and only
// a * branch is permissive.
func (t *Template) Match(name string, ph map[string][]string) (map[string]string, bool) {
	if len(t.parts) == 0 || !strings.HasPrefix(name, t.leadingLiteral()) {
		return nil, false
	}

	var pattern strings.Builder
	var order []string
	pattern.WriteByte('^')
	for _, p := range t.parts {
		if !p.IsToken {
			pattern.WriteString(regexp.QuoteMeta(p.Text))
			continue
		}
		pattern.WriteByte('(')
		pattern.WriteString(valuePattern(ph[p.Text]))
		pattern.WriteByte(')')
		order = append(order, p.Text)
	}
	pattern.WriteByte('$')

	re, err := regexp.Compile(pattern.String())
	if err != nil {
		return nil, false
	}
	m := re.FindStringSubmatch(name)
	if m == nil {
		return nil, false
	}

	vars := make(map[string]string, len(order))
	for i, key := range order {
		// A placeholder used twice must capture the same value both times;
		// RE2 has no backreference to enforce it in the pattern.
		if prev, seen := vars[key]; seen && prev != m[i+1] {
			return nil, false
		}
		vars[key] = m[i+1]
	}
	return vars, true
}

// Enumerate returns every concrete variant, vars included so callers do not
// have to run each generated name back through Match. Open-ended values cannot
// be enumerated and are skipped.
func (t *Template) Enumerate(ph map[string][]string) []Variant {
	names := t.Names()
	combos := []map[string]string{{}}
	for _, key := range names {
		var concrete []string
		for _, v := range ph[key] {
			if v != "*" {
				concrete = append(concrete, v)
			}
		}
		if len(concrete) == 0 {
			continue
		}
		next := make([]map[string]string, 0, len(combos)*len(concrete))
		for _, base := range combos {
			for _, v := range concrete {
				combo := make(map[string]string, len(base)+1)
				for k, bv := range base {
					combo[k] = bv
				}
				combo[key] = v
				next = append(next, combo)
			}
		}
		combos = next
	}

	var out []Variant
	seen := map[string]bool{}
	for _, combo := range combos {
		name, err := t.Fill(combo)
		if err != nil || name == "" || seen[name] {
			continue
		}
		seen[name] = true
		out = append(out, Variant{Name: name, Vars: combo})
	}
	return out
}

// leadingLiteral is the fixed prefix every match must start with, used to skip
// compiling a pattern that cannot match.
func (t *Template) leadingLiteral() string {
	if len(t.parts) > 0 && !t.parts[0].IsToken {
		return t.parts[0].Text
	}
	return ""
}

// valuePattern turns a #PH: value set into an alternation. Order does not
// matter: RE2 finds a match under an anchored pattern whenever one exists.
func valuePattern(values []string) string {
	var alts []string
	for _, v := range values {
		if v == "*" {
			alts = append(alts, `[^/]+`) // must never swallow a separator
			continue
		}
		alts = append(alts, regexp.QuoteMeta(v))
	}
	if len(alts) == 0 {
		return `[^/]+`
	}
	return strings.Join(alts, "|")
}

// ValidateTemplate reports disagreements between a #TARGET: and its #PH:
// declarations. Returns one message per problem, or nil when the template is
// sound.
func ValidateTemplate(target string, ph map[string][]string) []string {
	if len(ph) == 0 {
		return nil
	}
	if target == "" {
		return []string{"#PH: is declared but there is no #TARGET:"}
	}

	t := NewTemplate(target)
	inTarget := t.Names()

	var problems []string
	for _, name := range slices.Sorted(maps.Keys(ph)) {
		if !slices.Contains(inTarget, name) {
			problems = append(problems,
				fmt.Sprintf("#PH:%s is declared but {%s} is missing from #TARGET:", name, name))
		}
	}
	for _, name := range inTarget {
		if _, ok := ph[name]; !ok {
			problems = append(problems,
				fmt.Sprintf("#TARGET: uses {%s} but there is no #PH:%s declaration", name, name))
		}
	}

	// Two open-ended placeholders with no literal between them can match the
	// same name in more than one way, so nothing can recover the values.
	for i := 0; i+1 < len(t.parts); i++ {
		a, b := t.parts[i], t.parts[i+1]
		if a.IsToken && b.IsToken && isOpen(ph[a.Text]) && isOpen(ph[b.Text]) {
			problems = append(problems,
				fmt.Sprintf("#TARGET: has {%s}{%s} adjacent and both are open-ended", a.Text, b.Text))
		}
	}
	return problems
}

func isOpen(values []string) bool { return slices.Contains(values, "*") }

// parsePH parses a "#PH:name:values" body into its name and concrete values.
func parsePH(value string) (string, []string, bool) {
	name, raw, ok := strings.Cut(value, ":")
	name = strings.TrimSpace(name)
	if !ok || name == "" {
		return "", nil, false
	}
	values := ParseValues(raw)
	if len(values) == 0 {
		return "", nil, false
	}
	return name, values, true
}

var phRange = regexp.MustCompile(`^(\d+)-(\d+)$`)

// ParseValues expands a value list into concrete values, ranges included.
//
// The separator sets the order, because values[0] is the default: a comma
// list sorts newest first, since for versions that is the sensible default,
// while a pipe list keeps the order written, since for labels "larger" means
// nothing. A * is always last, being a fallback rather than a value. Where both
// separators appear the pipe wins; validate_recipes.py rejects the mixture.
//
// Exported because the dialect is shared: recipe #PH: and helper #VALUE: lists
// are written the same way, and a second implementation would let them drift.
func ParseValues(raw string) []string {
	raw = strings.TrimSpace(raw)
	sep, sorted := ",", true
	if strings.Contains(raw, "|") {
		sep, sorted = "|", false
	}

	var values []string
	open := false
	seen := map[string]bool{}
	for _, tok := range strings.Split(raw, sep) {
		tok = strings.TrimSpace(tok)
		switch {
		case tok == "":
			continue
		case tok == "*":
			open = true
			continue
		}
		for _, v := range expandRange(tok) {
			if !seen[v] {
				seen[v] = true
				values = append(values, v)
			}
		}
	}

	if sorted {
		slices.SortStableFunc(values, func(a, b string) int { return CompareVersions(b, a) })
	}
	if open {
		values = append(values, "*")
	}
	return values
}

// expandRange turns "22-49" into every integer between, or returns the token.
func expandRange(tok string) []string {
	m := phRange.FindStringSubmatch(tok)
	if m == nil {
		return []string{tok}
	}
	lo, err1 := strconv.Atoi(m[1])
	hi, err2 := strconv.Atoi(m[2])
	if err1 != nil || err2 != nil {
		return []string{tok}
	}
	if lo > hi {
		lo, hi = hi, lo
	}
	out := make([]string, 0, hi-lo+1)
	for v := lo; v <= hi; v++ {
		out = append(out, strconv.Itoa(v))
	}
	return out
}
