package record

import (
	"fmt"
	"strings"

	"github.com/Justype/condatainer/catalog"
)

// group is a field's position in the fixed order a record is written in.
type group int

// The source section is split in two so recipe-then-from is enforced by the same
// comparison every other group uses.
const (
	groupType group = iota
	groupEnv
	groupRecipe
	groupFrom
	groupPH
	groupDep
)

// groupOf maps a line's key to its group, and reports whether the key is known.
// An unknown key is rejected rather than ignored: a v1 record is fully specified,
// so a key this build does not know means a format that should have changed its
// tag.
func groupOf(key string) (group, bool) {
	switch key {
	case "type":
		return groupType, true
	case "env":
		return groupEnv, true
	case "recipe":
		return groupRecipe, true
	case "from":
		return groupFrom, true
	case "ph":
		return groupPH, true
	case "dep":
		return groupDep, true
	}
	return 0, false
}

// Parse reads a record back. It accepts only canonical form: the exact bytes
// Marshal produces. Groups out of order, an unsorted group, a blank line, or
// trailing whitespace are all rejected, because a record that parses but does not
// re-marshal identically would hash to something other than the key it carries.
func Parse(data []byte) (Record, error) {
	text := string(data)
	if text == "" {
		return Record{}, fmt.Errorf("%w: empty", ErrInvalid)
	}
	if !strings.HasSuffix(text, "\n") {
		return Record{}, fmt.Errorf("%w: no trailing newline", ErrInvalid)
	}
	lines := strings.Split(strings.TrimSuffix(text, "\n"), "\n")

	var r Record
	switch lines[0] {
	case KindIdentity.Tag():
		r.Kind = KindIdentity
	case KindEquiv.Tag():
		r.Kind = KindEquiv
	default:
		return Record{}, fmt.Errorf("%w: %q is not a record format line", ErrInvalid, lines[0])
	}

	var (
		last     = groupType
		haveType bool
	)
	for i, line := range lines[1:] {
		no := i + 2
		if line == "" {
			return Record{}, fmt.Errorf("%w: line %d is blank", ErrInvalid, no)
		}
		if line != strings.TrimRight(line, " \t") {
			return Record{}, fmt.Errorf("%w: line %d has trailing whitespace", ErrInvalid, no)
		}
		key, value, ok := strings.Cut(line, "=")
		if !ok {
			return Record{}, fmt.Errorf("%w: line %d is not key=value: %q", ErrInvalid, no, line)
		}
		g, known := groupOf(key)
		if !known {
			return Record{}, fmt.Errorf("%w: line %d has unknown key %q", ErrInvalid, no, key)
		}
		if g < last {
			return Record{}, fmt.Errorf("%w: line %d: %s appears after a later field", ErrInvalid, no, key)
		}
		last = g

		switch g {
		case groupType:
			if haveType {
				return Record{}, fmt.Errorf("%w: line %d: type appears twice", ErrInvalid, no)
			}
			haveType = true
			r.Type = catalog.Type(value)
		case groupEnv:
			k, v, ok := strings.Cut(value, "=")
			if !ok {
				return Record{}, fmt.Errorf("%w: line %d: env is not KEY=value", ErrInvalid, no)
			}
			if n := len(r.Env); n > 0 && r.Env[n-1].Key >= k {
				return Record{}, fmt.Errorf("%w: line %d: env %s is out of order", ErrInvalid, no, k)
			}
			r.Env = append(r.Env, Env{Key: k, Value: v})
		case groupRecipe:
			if r.Recipe != "" {
				return Record{}, fmt.Errorf("%w: line %d: recipe appears twice", ErrInvalid, no)
			}
			r.Recipe = value
		case groupFrom:
			if r.From != "" {
				return Record{}, fmt.Errorf("%w: line %d: from appears twice", ErrInvalid, no)
			}
			r.From = value
		case groupPH:
			n, v, ok := strings.Cut(value, "=")
			if !ok {
				return Record{}, fmt.Errorf("%w: line %d: ph is not name=value", ErrInvalid, no)
			}
			if c := len(r.Placeholders); c > 0 && r.Placeholders[c-1].Name >= n {
				return Record{}, fmt.Errorf("%w: line %d: placeholder %s is out of order", ErrInvalid, no, n)
			}
			r.Placeholders = append(r.Placeholders, Placeholder{Name: n, Value: v})
		case groupDep:
			// Repeated identical lines are legal, so this compares for a strict
			// decrease rather than requiring every line to differ.
			if c := len(r.Deps); c > 0 && r.Deps[c-1].text() > value {
				return Record{}, fmt.Errorf("%w: line %d: dependency %q is out of order", ErrInvalid, no, value)
			}
			fields := strings.Split(value, " ")
			r.Deps = append(r.Deps, Dep{Type: catalog.Type(fields[0]), Fields: fields[1:]})
		}
	}

	if !haveType {
		return Record{}, fmt.Errorf("%w: no type line", ErrInvalid)
	}
	if err := validate(r); err != nil {
		return Record{}, err
	}
	return r, nil
}
