package catalog

import (
	"context"
	"fmt"
	"slices"
	"strings"
)

// Have reports the versions of a name a caller already has: an image directory
// for one tool, a modulefile tree for another.
//
// It answers with versions, not with a verdict. Deciding whether 1.23.1>=1.10
// is satisfied is the resolver's, so that both callers apply one rule.
type Have func(name string) []string

// Node is one module in a build plan.
type Node struct {
	Dep    Dep
	Entry  *Entry // nil when no source provides it
	Source *Source
	Vars   map[string]string // set when the name came from a template

	// Installed is the version Have reported, when one satisfies the dep.
	Installed string
}

// Name is the module path this node resolves to. The constraint is not part of
// it: ">=1.16" selected the version, it does not name the module.
func (n Node) Name() string {
	if n.Entry != nil && !n.Entry.IsTemplate {
		return n.Entry.Name
	}
	if n.Dep.Version == "" {
		return n.Dep.Name
	}
	return n.Dep.Name + "/" + n.Dep.Version
}

// Plan is a build order and the subset of it still to act on.
type Plan struct {
	Order   []Node // topological, dependencies first
	Missing []Node // the nodes the caller must build, pull or hand to a fallback
}

// Resolve walks the dependency graph of roots and returns it in build order.
//
// It runs over Entry, so a whole tree is computed from the index without
// fetching a recipe. Nodes carry their Type, which is how a caller filters what
// it cannot build, and a node no source provides carries a nil Entry rather
// than failing the walk — that is the normal path for a conda package.
func (c Catalog) Resolve(ctx context.Context, roots []string, have Have) (*Plan, error) {
	r := &resolver{cat: c, have: have, seen: map[string]bool{}, onStack: map[string]bool{}}
	for _, root := range roots {
		dep, err := ParseDep(root)
		if err != nil {
			return nil, err
		}
		if err := r.visit(ctx, dep, nil); err != nil {
			return nil, err
		}
	}

	plan := &Plan{Order: r.order}
	for _, n := range r.order {
		if n.Installed == "" {
			plan.Missing = append(plan.Missing, n)
		}
	}
	return plan, nil
}

type resolver struct {
	cat     Catalog
	have    Have
	seen    map[string]bool
	onStack map[string]bool
	order   []Node
	stack   []string
}

// visit resolves one dep and everything below it, appending in dependency-first
// order. vars carries the parent's expansion, since a template dep is
// substituted before it has a name at all.
func (r *resolver) visit(ctx context.Context, dep Dep, vars map[string]string) error {
	name := replaceVars(dep.String(), vars)
	dep, err := ParseDep(name)
	if err != nil {
		return err
	}
	key := dep.Name + "/" + dep.Version

	if r.onStack[key] {
		return fmt.Errorf("catalog: dependency cycle: %s", strings.Join(append(r.stack, key), " -> "))
	}
	if r.seen[key] {
		return nil // diamond: already placed, and placed before its dependents
	}

	node, err := r.resolveOne(ctx, dep)
	if err != nil {
		return err
	}

	r.onStack[key], r.stack = true, append(r.stack, key)
	if node.Entry != nil && node.Installed == "" {
		for _, raw := range node.Entry.Deps {
			child, err := ParseDep(raw)
			if err != nil {
				return err
			}
			if err := r.visit(ctx, child, node.Vars); err != nil {
				return err
			}
		}
	}
	r.onStack[key], r.stack = false, r.stack[:len(r.stack)-1]

	r.seen[key] = true
	r.order = append(r.order, node)
	return nil
}

// resolveOne picks the version a dep gets and finds what provides it.
//
// The preferred version wins when it resolves. Otherwise the newest installed
// in range, so a bare #DEP: does not rebuild the moment upstream moves; then
// the newest available; then nothing provides it, and the caller's fallback
// picks its own.
func (r *resolver) resolveOne(ctx context.Context, dep Dep) (Node, error) {
	node := Node{Dep: dep}

	if v := r.installed(dep); v != "" {
		node.Installed = v
		if dep.Version == "" {
			node.Dep.Version = v
		}
	}

	m, found, err := r.cat.Lookup(ctx, dep.Name+"/"+node.Dep.Version)
	if err != nil {
		return node, err
	}
	if !found && node.Dep.Version == "" {
		if v := r.newestAvailable(ctx, dep); v != "" {
			node.Dep.Version = v
			m, found, err = r.cat.Lookup(ctx, dep.Name+"/"+v)
			if err != nil {
				return node, err
			}
		}
	}
	if found {
		node.Entry, node.Source, node.Vars = m.Entry, m.Source, m.Vars
	}
	return node, nil
}

// installed returns the newest version the caller has that the dep admits.
func (r *resolver) installed(dep Dep) string {
	if r.have == nil {
		return ""
	}
	candidates := slices.Clone(r.have(dep.Name))
	slices.SortStableFunc(candidates, func(a, b string) int { return CompareVersions(b, a) })
	for _, v := range candidates {
		if dep.Satisfies(v) {
			return v
		}
	}
	return ""
}

// newestAvailable returns the newest version the catalog offers in range.
func (r *resolver) newestAvailable(ctx context.Context, dep Dep) string {
	versions, err := r.cat.Versions(ctx, dep.Name)
	if err != nil {
		return ""
	}
	for _, v := range versions {
		if dep.Satisfies(v) {
			return v
		}
	}
	return ""
}
