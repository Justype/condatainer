// Package restore makes a project's locked identities locally available.
//
// Planning is separated from acquisition on purpose. A plan is computed from
// the checkout and what is already installed, writes nothing, contacts nothing,
// and acquires no lock — so `--dry-run` can print exactly what would happen
// before anything expensive or irreversible starts.
package restore

import (
	"fmt"
	"os"
	"path/filepath"
	"slices"
	"sort"
	"strings"

	"github.com/Justype/condatainer/internal/artifact/capsule"
	"github.com/Justype/condatainer/internal/artifact/compare"
	"github.com/Justype/condatainer/internal/artifact/key"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/project"
	"github.com/Justype/condatainer/internal/project/lock"
	"github.com/Justype/condatainer/internal/store"
)

// Match is which of an artifact's two keys a restored copy has to agree with.
// Defined by internal/project so restore and project-aware execution cannot
// disagree about what a lock means.
type Match = project.Match

const (
	// MatchEquivalence accepts anything that can substitute for what was locked.
	MatchEquivalence = project.MatchEquivalence
	// MatchIdentity accepts only the exact build the lock names.
	MatchIdentity = project.MatchIdentity
)

// Action is how one artifact will be made available.
type Action string

const (
	// ActionAdopt reuses a verified artifact that is already here.
	ActionAdopt Action = "adopt"
	// ActionFetch downloads a recorded registry remote.
	ActionFetch Action = "fetch"
	// ActionBuild rebuilds from the vendored sources.
	ActionBuild Action = "build"
	// ActionUnavailable is an artifact that is neither here nor obtainable. A
	// plan carrying one is never complete, so nothing runs on it.
	ActionUnavailable Action = "unavailable"
)

// Step is one artifact's place in a restore.
type Step struct {
	// Artifact is the vendored directory, relative to the lock directory.
	Artifact string `json:"artifact"`
	Name     string `json:"name"`
	// Identity is what the lock records.
	Identity string `json:"identity"`
	// Found is the identity actually present, set only when an equivalent
	// artifact stands in for the locked one. Empty for an exact hit, so a reader
	// can tell a substitution from a match without comparing two digests.
	Found  string `json:"found,omitempty"`
	Action Action `json:"action"`
	// Direct reports that a pin mounts this artifact. A closure-only step
	// is a build dependency, materialized only because something above it must be
	// rebuilt.
	Direct bool `json:"direct"`
	// Requests are the pin keys this artifact answers, sorted. Empty for
	// a closure-only step.
	Requests []string `json:"requests,omitempty"`
	// Destination is the project-relative path this artifact must end up at,
	// set only for a `path:` pin. Empty means the artifact is addressed
	// by identity and goes wherever the store's destination rule puts it.
	//
	// A project path is an output, not somewhere to look: a step carrying one
	// is never satisfied by, and never installed into, a shared images root.
	Destination string `json:"destination,omitempty"`
	// Replaces describes a different artifact already sitting at Destination
	// that this step will overwrite — its identity, or "an unreadable file" when
	// the bytes there are not a readable artifact. Empty when the path is free
	// or holds something that answers the lock.
	//
	// A restore renames over that file, so a plan that did not say this would
	// preview an overwrite as an ordinary build.
	Replaces string `json:"replaces,omitempty"`
	// Path and Layout are set when the artifact was found locally. For a
	// destination step Path is that destination, already present and verified.
	Path   string       `json:"path,omitempty"`
	Layout store.Layout `json:"layout,omitempty"`
	// Remotes are the recorded fetch locations, in retry order.
	Remotes []lock.Remote `json:"remotes,omitempty"`
	// DependsOn are the artifacts that must exist first, sorted.
	DependsOn []string `json:"depends_on,omitempty"`
	// RequiresInput reports that a rebuild would prompt, for #INPUT: or #SOURCE: ask: answers.
	RequiresInput bool `json:"requires_input,omitempty"`
	// Arch is the artifact's recorded architecture, or "noarch".
	Arch string `json:"arch,omitempty"`
	// StoreOnly keeps this step out of the bare name because a pin in the
	// same plan answers to it. Set only on a build dependency being installed.
	StoreOnly bool `json:"store_only,omitempty"`
	// Unbuildable reports an artifact with no sources to rebuild from, so a
	// registry copy is the only thing that can produce it.
	Unbuildable bool `json:"unbuildable,omitempty"`
}

// verdict is how a step's result relates to the lock, decided by whether an
// equivalent artifact stood in for the locked identity.
func (s Step) verdict() compare.Verdict {
	if s.Found != "" {
		return compare.Equivalent
	}
	return compare.Exact
}

// Plan is a dependency-first restore, with everything that would prevent it.
type Plan struct {
	Root  string `json:"root"`
	Match Match  `json:"match"`
	// Steps are ordered so every step's DependsOn precede it.
	Steps    []Step   `json:"steps"`
	Problems []string `json:"problems,omitempty"`
}

// Complete reports whether the plan can run as computed.
func (p *Plan) Complete() bool { return len(p.Problems) == 0 }

// Work reports how many steps would acquire something. A plan of nothing but
// adoptions is already satisfied.
func (p *Plan) Work() int {
	var n int
	for _, step := range p.Steps {
		if step.Action != ActionAdopt {
			n++
		}
	}
	return n
}

// Options tunes planning.
type Options struct {
	// Match is which key a local copy has to agree with. Empty means the lock's.
	Match Match
	// SkipPrebuilt ignores every recorded remote and plans a local build
	// instead. It is not a claim that the machine is offline: a build needs the
	// network too, and generally more of it — a Conda replay downloads every
	// pinned package URL, a recipe fetches its own sources, and resolving an
	// absent base bootstraps one from a registry. On a compute node without
	// egress the answer is the proxy (internal/runtime/proxy), not this flag.
	SkipPrebuilt bool
	// KeepBuildDeps installs a newly produced build dependency through the ordinary
	// destination rule instead of discarding it with the restore. An input that
	// was already installed is adopted in place either way.
	KeepBuildDeps bool
	// Replace allows a project path holding something that does not answer the
	// lock to be overwritten. Without it that is a planning problem, refused
	// before anything is acquired.
	//
	// A project path is the only place a restore can destroy a file: an artifact
	// addressed by name goes wherever the store's destination rule puts it, which
	// never takes a name another identity already holds.
	Replace bool
	// SearchDirs overrides the configured image roots.
	SearchDirs []string
	// HostArch is the architecture to plan for. Empty means this machine.
	HostArch string
	// SubmitJobs allows a rebuild carrying scheduler directives to be handed to
	// the scheduler instead of run here. Ignored inside a job, which is already
	// the machine the work was sent to.
	SubmitJobs bool
	// Only narrows the restore to one vendored artifact and its closure. A
	// submitted job carries it so the job produces exactly what it was sent for
	// and leaves every other artifact to whoever asked for that one.
	Only string
	// Answers are the answers to a recipe's prompts per artifact path, #INPUT: first and then
	// each #SOURCE: ask: in declaration order.
	// Supplied per invocation because they are never recorded, and their
	// presence is what makes an interactive rebuild plannable at all.
	Answers map[string][]string
	// lookup finds a local copy of one locked artifact. Injected for tests.
	lookup lookupFunc
	// lookupAt finds the artifact at one exact path. Injected for tests.
	lookupAt lookupAtFunc
	// lookupInput finds something that may satisfy one dependency edge.
	// Injected for tests.
	lookupInput lookupInputFunc
}

// lookupFunc finds a local artifact answering to keys under a match mode.
type lookupFunc = project.LookupFunc

// lookupAtFunc reports whether the artifact at one exact path satisfies keys.
type lookupAtFunc = project.LookupAtFunc

// lookupInputFunc finds a local artifact that may satisfy one dependency edge.
type lookupInputFunc func(dep meta.Dependency, match Match, dirs []string) (store.Candidate, bool)

func (o Options) resolver() lookupFunc {
	if o.lookup != nil {
		return o.lookup
	}
	return project.LookupLocal
}

func (o Options) pathResolver() lookupAtFunc {
	if o.lookupAt != nil {
		return o.lookupAt
	}
	return project.LookupAt
}

func (o Options) inputResolver() lookupInputFunc {
	if o.lookupInput != nil {
		return o.lookupInput
	}
	if o.lookup != nil {
		// A caller that injected only a name lookup answers for build dependencies the
		// same way, so it can never fall through to the real image roots.
		return func(dep meta.Dependency, match Match, dirs []string) (store.Candidate, bool) {
			return o.lookup(dep.Name, meta.Keys{Identity: dep.Identity, Equiv: dep.Equiv}, match, dirs)
		}
	}
	return project.LookupInput
}

// Compute plans a restore from a verified lock.
//
// It runs the same checkout-local validation `project validate` does and stops
// there when the lock is not sound: planning an acquisition against an
// inconsistent specification would only fail later, further from the cause.
//
// Nothing here writes, fetches, or locks. Every artifact that is already
// present is adopted, and only a genuine miss becomes a fetch or a rebuild.
func Compute(root string, l *lock.Lock, opts Options) *Plan {
	match := opts.Match
	if match == "" {
		match = l.Match
	}
	match = match.Normalize()
	plan := &Plan{Root: root, Match: match}

	verified, problems := lock.Verify(root, l)
	if len(problems) > 0 {
		for _, problem := range problems {
			plan.Problems = append(plan.Problems, problem.String())
		}
		return plan
	}

	hostArch := opts.HostArch
	if hostArch == "" {
		hostArch = meta.NativeArch()
	}
	resolve, resolveAt, resolveInput := opts.resolver(), opts.pathResolver(), opts.inputResolver()

	// How each artifact is depended on, which decides what may stand in for it
	// when it is only a build dependency.
	edges := inboundEdges(verified)

	named, destinations, conflicts := split(l)
	plan.Problems = append(plan.Problems, conflicts...)

	order, cycle := topological(verified)
	if cycle != "" {
		plan.Problems = append(plan.Problems, cycle)
		return plan
	}

	// Steps and their would-be problems are collected first: a closure-only step
	// may be pruned below, and a pruned step's problems are not the plan's.
	var steps []Step
	pending := map[string][]string{}

	for _, artifact := range order {
		entry := verified.Entries[artifact]
		base := Step{
			Artifact: artifact,
			Name:     entry.Manifest.Name,
			Identity: entry.Identity.Digest(),
			Remotes:  l.Remotes[artifact],
			Arch:     entry.Manifest.Platform.Arch,

			RequiresInput: entry.Manifest.Source.RequiresInput,
			Unbuildable:   key.IsSnapshot(entry.Manifest),
			DependsOn:     dependsOn(verified, entry),
		}
		keys := meta.Keys{Identity: entry.Identity, Equiv: entry.Equiv}

		// One step per project destination: two paths selecting one artifact are
		// two files to produce, not one artifact to place twice.
		if placements := destinations[artifact]; len(placements) > 0 {
			for _, destination := range sortedKeys(placements) {
				step := base
				step.Direct, step.Destination = true, destination
				step.Requests = sorted(placements[destination])
				at := filepath.Join(root, filepath.FromSlash(destination))
				step, reasons := plan.classify(step, keys, hostArch, opts, func() (store.Candidate, bool) {
					return resolveAt(at, entry.Manifest.Name, keys, match)
				})
				if step.Action != ActionAdopt {
					step.Replaces = occupantAt(at)
				}
				if step.Replaces != "" && !opts.Replace {
					reasons = append(reasons, fmt.Sprintf(
						"%s already holds %s, which is not what %s pins; pass --replace to overwrite it",
						destination, step.Replaces, entry.Manifest.Name))
				}
				steps = append(steps, step)
				pending[artifact] = append(pending[artifact], reasons...)
			}
			continue
		}

		step := base
		step.Direct = len(named[artifact]) > 0
		step.Requests = sorted(named[artifact])
		// A pin has to answer for its own keys — someone asked for it by
		// name. A build dependency only has to leave its dependent's equivalence
		// unchanged, which is what its edge's role states.
		find := func() (store.Candidate, bool) {
			if step.Direct {
				return resolve(entry.Manifest.Name, keys, match, opts.SearchDirs)
			}
			dep, ok := edges[artifact]
			if !ok {
				return resolve(entry.Manifest.Name, keys, match, opts.SearchDirs)
			}
			return resolveInput(dep, match, opts.SearchDirs)
		}
		step, reasons := plan.classify(step, keys, hostArch, opts, find)
		steps = append(steps, step)
		pending[artifact] = append(pending[artifact], reasons...)
	}

	baseArtifact := l.Pins[lock.BaseKey].Artifact
	kept, reported := prune(steps, pending)
	yieldSharedNames(kept, opts)
	if opts.Only != "" {
		var problem string
		if kept, problem = restrict(kept, opts.Only, baseArtifact); problem != "" {
			plan.Problems = append(plan.Problems, problem)
			return plan
		}
		reported = problemsFor(kept, pending)
	}
	plan.Steps = reorderBaseFirst(kept, baseArtifact)
	plan.Problems = append(plan.Problems, reported...)
	return plan
}

// reorderBaseFirst moves the project's root to the front of an otherwise
// dependency-ordered plan.
//
// Nothing in the manifest graph points at it — no #DEP: may name an os
// artifact — so topological order leaves it wherever it falls alphabetically.
// Every conda or script build that resolves its own root from the lock
// (build.LockedSpec.Base, set from Run's `available` map) needs it to have
// already run, and an os build never has dependencies of its own to violate
// by moving first.
func reorderBaseFirst(steps []Step, baseArtifact string) []Step {
	if baseArtifact == "" {
		return steps
	}
	for i, step := range steps {
		if step.Artifact != baseArtifact {
			continue
		}
		if i == 0 {
			return steps
		}
		out := make([]Step, 0, len(steps))
		out = append(out, step)
		out = append(out, steps[:i]...)
		out = append(out, steps[i+1:]...)
		return out
	}
	return steps
}

// restrict drops every step outside one artifact's closure, keeping the order.
//
// Directness is not rewritten. An artifact inside the closure that some
// pin also names by itself stays a pin and is still installed where
// the destination rule puts it, because that is what it is — narrowing what a
// restore covers says nothing about what its members are.
//
// The base artifact is always kept, whether or not only's closure reaches it:
// nothing in the manifest graph ever depends on it (#DEP: cannot name an os
// artifact), yet a submitted job re-enters Compute with --only set and still
// has to resolve its own root the same way an unrestricted restore does — see
// reorderBaseFirst and internal/project/restore/restore.go's use of Run's
// `available` map.
func restrict(steps []Step, only, baseArtifact string) ([]Step, string) {
	keep := map[string]bool{only: true}
	if baseArtifact != "" {
		keep[baseArtifact] = true
	}
	found := false
	for i := len(steps) - 1; i >= 0; i-- {
		step := steps[i]
		if step.Artifact == only {
			found = true
		}
		if !keep[step.Artifact] {
			continue
		}
		for _, dependency := range step.DependsOn {
			keep[dependency] = true
		}
	}
	if !found {
		return nil, fmt.Sprintf("--only %s is not an artifact this restore covers", only)
	}
	kept := make([]Step, 0, len(steps))
	for _, step := range steps {
		if keep[step.Artifact] {
			kept = append(kept, step)
		}
	}
	return kept, ""
}

// yieldSharedNames gives the bare name to the pin whenever a kept build
// dependency shares it, marking the dependency store-only.
//
// The flat name answers `exec -o`, `list`, and every other checkout, so it
// belongs to what someone asked for by name rather than to what a rebuild
// happened to need. Left to the store alone the winner is whoever installs
// first, and dependencies are built first — exactly backwards. Deciding it here
// makes it a property of the plan, so a dry run shows the same answer the
// restore will reach.
//
// Only two identities of one name ever contend. One identity is one step, and a
// discarded build dependency is never installed at all.
func yieldSharedNames(steps []Step, opts Options) {
	pinned := map[string]bool{}
	for _, step := range steps {
		if step.Direct && step.Destination == "" {
			pinned[step.Name] = true
		}
	}
	for i, step := range steps {
		if step.Direct || step.Destination != "" || isTransient(step, opts) {
			continue
		}
		steps[i].StoreOnly = pinned[step.Name]
	}
}

// occupantAt describes the artifact already at a project destination that does
// not answer the lock, or "" when the path is free.
//
// Read only for a step that is going to write there: a restore renames over
// whatever regular .sqf it finds, so this is the plan's one chance to say what
// is about to be lost. A file whose metadata cannot be read is still reported —
// it is still about to be replaced.
func occupantAt(path string) string {
	info, err := os.Lstat(path)
	if err != nil || !info.Mode().IsRegular() {
		return ""
	}
	artifact, err := compare.Read(path)
	if err != nil {
		return "an unreadable file"
	}
	return artifact.IdentityRef().Digest()
}

// problemsFor collects the problems the given steps carry, in step order.
func problemsFor(steps []Step, problems map[string][]string) []string {
	var out []string
	for _, step := range steps {
		out = append(out, problems[step.Artifact]...)
	}
	return out
}

// split separates pins by where the artifact has to end up: named ones go
// wherever the store's destination rule puts them, `path:` ones go to a declared
// place inside the project. It also reports any artifact pinned both ways.
//
// The two cannot be planned together. A named pin is satisfied by a copy
// in any readable images root; a project path is a file that must exist at one
// exact location. An artifact asked for both has two destinations and no way to
// choose, so it is a lock to fix rather than a case to resolve.
//
// The returned slices hold at most one request each in practice. A name key must
// equal the artifact's own name, so two name keys cannot select one artifact,
// and a destination key *is* its destination. They stay slices because that is
// what reports the both-ways conflict above without a second pass.
func split(l *lock.Lock) (named map[string][]string, destinations map[string]map[string][]string, problems []string) {
	named = map[string][]string{}
	destinations = map[string]map[string][]string{}
	for _, request := range l.Requests() {
		artifact := l.Pins[request].Artifact
		destination, isPath := strings.CutPrefix(request, lock.PathPrefix)
		if !isPath {
			named[artifact] = append(named[artifact], request)
			continue
		}
		if destinations[artifact] == nil {
			destinations[artifact] = map[string][]string{}
		}
		destinations[artifact][destination] = append(destinations[artifact][destination], request)
	}
	for artifact := range destinations {
		if len(named[artifact]) > 0 {
			problems = append(problems, fmt.Sprintf(
				"%s is pinned both by name (%s) and at a project path; it has one payload and cannot have two destinations",
				artifact, strings.Join(sorted(named[artifact]), ", ")))
		}
	}
	sort.Strings(problems)
	return named, destinations, problems
}

// classify decides how one step is satisfied and records it.
//
// found is the only thing that differs between a named artifact and a project
// path: the first is looked for across the images roots, the second only at its
// declared path. Everything after — architecture, remotes, interactive input —
// applies to a miss either way.
// classify decides how one step is satisfied, reporting the problems that would
// follow if it is kept.
//
// Problems are returned rather than recorded because a closure-only step may be
// pruned: a platform mismatch or a missing #INPUT: answer for an artifact
// nothing needs is not a reason to refuse the plan.
func (p *Plan) classify(step Step, keys meta.Keys, hostArch string, opts Options,
	found func() (store.Candidate, bool)) (Step, []string) {

	if candidate, ok := found(); ok {
		step.Action, step.Path, step.Layout = ActionAdopt, candidate.Path, candidate.Layout
		if candidate.Identity != keys.Identity {
			step.Found = candidate.Identity.Digest()
		}
		return step, nil
	}

	// Only a miss has to be acquired, so only a miss has to run here.
	var problems []string
	if reason := runsOn(step.Arch, hostArch); reason != "" {
		problems = append(problems, fmt.Sprintf("%s: %s", step.Name, reason))
	}
	switch {
	// --no-prebuilt chooses building over fetching, so it has nothing to say
	// about an artifact that cannot be built: skipping the fetch there would
	// refuse the restore rather than take the other route.
	case len(step.Remotes) > 0 && (!opts.SkipPrebuilt || step.Unbuildable):
		step.Action = ActionFetch
	case step.Unbuildable:
		// Refused while planning rather than attempted: a frozen environment was
		// captured, not built, so no source it could be rebuilt from exists
		// anywhere.
		step.Action = ActionUnavailable
		problems = append(problems, fmt.Sprintf(
			"%s is a frozen environment and cannot be rebuilt; it is only obtainable from a registry, so publish it with `condatainer project push` from a checkout that has it",
			step.Name))
	default:
		step.Action = ActionBuild
		if step.RequiresInput && len(opts.Answers[step.Artifact]) == 0 {
			problems = append(problems, fmt.Sprintf(
				"%s must be rebuilt and its recipe asks for input (#INPUT: or a #SOURCE: ask:), which needs a terminal", step.Name))
		}
	}
	return step, problems
}

// prune drops closure-only steps nothing needs, and returns the problems the
// surviving steps carry.
//
// A closure entry is a *build dependency*: it exists so its dependent can be rebuilt
// to the exact identity the lock records. When that dependent is adopted or
// fetched instead, the input is never opened, so acquiring it would be pure
// waste — and its problems would refuse a plan that is in fact fine.
//
// Need flows backwards along dependency edges from steps that will actually
// build, which is why this walks the topological order in reverse. A fetch does
// not propagate need: it downloads a finished artifact and opens none of its
// inputs.
func prune(steps []Step, problems map[string][]string) ([]Step, []string) {
	needed := make(map[string]bool, len(steps))
	for i := len(steps) - 1; i >= 0; i-- {
		step := steps[i]
		if step.Direct {
			needed[step.Artifact] = true
		}
		if needed[step.Artifact] && step.Action == ActionBuild {
			for _, dependency := range step.DependsOn {
				needed[dependency] = true
			}
		}
	}
	kept := make([]Step, 0, len(steps))
	var reported []string
	for _, step := range steps {
		if !needed[step.Artifact] {
			continue
		}
		kept = append(kept, step)
		reported = append(reported, problems[step.Artifact]...)
	}
	return kept, reported
}

func sorted(in []string) []string {
	out := slices.Clone(in)
	sort.Strings(out)
	return out
}

func sortedKeys[V any](in map[string]V) []string {
	out := make([]string, 0, len(in))
	for key := range in {
		out = append(out, key)
	}
	sort.Strings(out)
	return out
}

// runsOn reports why an artifact cannot be materialized on this host, or "".
// noarch runs anywhere; anything else runs only where it was built, because
// nothing infers portability — a recipe declares it or it is absent.
func runsOn(arch, hostArch string) string {
	switch arch {
	case "":
		return "the artifact records no architecture"
	case meta.ArchNone, hostArch:
		return ""
	default:
		return fmt.Sprintf("built for %s, this host is %s", arch, hostArch)
	}
}

// inboundEdges maps each artifact to the edge that depends on it, keeping the
// strictest role when several do.
//
// Strictness is the order of the roles themselves: data pins content, app pins a
// version, history pins nothing. One artifact reached as data by one dependent
// and as history by another must satisfy the data claim, or that dependent's
// equivalence would move.
func inboundEdges(verified *lock.Verified) map[string]meta.Dependency {
	rank := map[string]int{meta.RoleData: 2, meta.RoleApp: 1, meta.RoleHistory: 0}
	out := map[string]meta.Dependency{}
	for _, entry := range verified.Entries {
		for _, dep := range entry.Manifest.Dependencies {
			artifact := lock.EntryPath(capsule.EntryName(dep.Name, dep.Identity.Digest()))
			if held, ok := out[artifact]; !ok || rank[dep.Role] > rank[held.Role] {
				out[artifact] = dep
			}
		}
	}
	return out
}

// dependsOn is the artifact paths one entry's edges point at, sorted. Verify has
// already established that each resolves, so this only has to name them.
func dependsOn(verified *lock.Verified, entry *lock.Entry) []string {
	var out []string
	for _, dep := range entry.Manifest.Dependencies {
		for path, candidate := range verified.Entries {
			if candidate.Manifest.Name == dep.Name && candidate.Identity == dep.Identity {
				out = append(out, path)
				break
			}
		}
	}
	sort.Strings(out)
	return out
}

// topological orders the reachable closure so every artifact follows the ones
// it depends on, and reports a cycle rather than looping forever.
//
// A capsule cannot contain a cycle — a dependency image existed before the
// artifact that mounted it — but a hand-edited lock is untrusted input and gets
// checked anyway.
func topological(verified *lock.Verified) ([]string, string) {
	const (
		unvisited = 0
		active    = 1
		done      = 2
	)
	state := map[string]int{}
	var order []string
	var cycle string

	var visit func(artifact string, trail []string)
	visit = func(artifact string, trail []string) {
		switch state[artifact] {
		case done:
			return
		case active:
			if cycle == "" {
				cycle = fmt.Sprintf("dependency cycle: %s", cycleText(append(trail, artifact)))
			}
			return
		}
		state[artifact] = active
		entry, ok := verified.Entries[artifact]
		if ok {
			for _, dep := range dependsOn(verified, entry) {
				visit(dep, append(trail, artifact))
			}
		}
		state[artifact] = done
		if ok {
			order = append(order, artifact)
		}
	}

	reachable := make([]string, 0, len(verified.Reachable))
	for artifact := range verified.Reachable {
		reachable = append(reachable, artifact)
	}
	sort.Strings(reachable)
	for _, artifact := range reachable {
		visit(artifact, nil)
	}
	return order, cycle
}

func cycleText(trail []string) string {
	out := trail[0]
	for _, artifact := range trail[1:] {
		out += " → " + artifact
	}
	return out
}
