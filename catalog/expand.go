package catalog

import (
	"fmt"
	"slices"
	"strings"
)

// Expand resolves a template recipe to one concrete variant.
//
// Every placeholder must have a value. A missing one is an error rather than a
// token left standing, which would run with {gencode_version} in a URL.
//
// The body is substituted too, because the recipe hash is taken over the
// expanded text — an unexpanded template would give every variant one key.
func Expand(r *Recipe, vars map[string]string) (*Recipe, error) {
	if !r.IsTemplate {
		return nil, fmt.Errorf("catalog: %s is not a template", r.Name)
	}
	t := NewTemplate(r.TargetTemplate)
	name, err := t.Fill(vars)
	if err != nil {
		return nil, err
	}

	out := *r
	out.Name = name
	out.IsTemplate = false
	out.TargetTemplate = ""
	out.Description = replaceVars(r.Description, vars)
	out.Text = []byte(replaceVars(string(r.Text), vars))

	// PH keeps the chosen value per placeholder: what the artifact records as
	// the vars it was built with.
	out.PH = make(map[string][]string, len(t.Names()))
	for _, key := range t.Names() {
		out.PH[key] = []string{vars[key]}
	}

	out.Deps = make([]string, len(r.Deps))
	for i, dep := range r.Deps {
		out.Deps[i] = replaceVars(dep, vars)
	}

	out.Env = make([]EnvVar, len(r.Env))
	for i, env := range r.Env {
		env.Segments = expandSegments(env.Segments, vars)
		env.Note = replaceVars(env.Note, vars)
		out.Env[i] = env
	}

	out.Inputs = slices.Clone(r.Inputs)
	for i := range out.Inputs {
		out.Inputs[i] = replaceVars(out.Inputs[i], vars)
	}
	out.Directives = slices.Clone(r.Directives)
	return &out, nil
}

// replaceVars substitutes {name} for each var. Only declared names are touched,
// so a shell ${VAR} in the body is left alone.
func replaceVars(s string, vars map[string]string) string {
	if s == "" || len(vars) == 0 {
		return s
	}
	for name, value := range vars {
		s = strings.ReplaceAll(s, "{"+name+"}", value)
	}
	return s
}

// expandSegments turns resolved tokens into literals. {prefix} has no var and
// survives, since where the payload mounts is not known until it is loaded.
func expandSegments(in []Segment, vars map[string]string) []Segment {
	out := make([]Segment, 0, len(in))
	for _, s := range in {
		if v, ok := vars[s.Text]; s.IsToken && ok {
			s = Segment{Text: v}
		}
		out = append(out, s)
	}
	return out
}
