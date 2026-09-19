package conda

import (
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strings"
)

// Shell selects which syntax ReactivateScript emits.
type Shell string

const (
	ShellBash Shell = "bash"
	ShellZsh  Shell = "zsh"
	ShellFish Shell = "fish"
)

// ReactivateScript returns a shell snippet that re-sources root's
// etc/conda/deactivate.d (reverse filename order) then etc/conda/
// activate.d (forward order), for a caller to eval after installing,
// updating, or removing packages. Globs fresh on each call, so a script
// for a package just removed is silently skipped rather than erroring.
// Only scripts ending in shell's own hook extension are sourced, so a
// bash/zsh-only *.sh hook is silently skipped under ShellFish rather than
// fed to fish's parser as invalid syntax.
func ReactivateScript(root string, shell Shell) string {
	ext := ".sh"
	if shell == ShellFish {
		ext = ".fish"
	}
	var b strings.Builder
	for _, f := range sortedHookScripts(filepath.Join(root, "etc", "conda", "deactivate.d"), ext, true) {
		writeSourceLine(&b, shell, f)
	}
	for _, f := range sortedHookScripts(filepath.Join(root, "etc", "conda", "activate.d"), ext, false) {
		writeScopedSourceLine(&b, shell, root, f)
	}
	return b.String()
}

// writeSourceLine sources file with no CONDA_PREFIX scoping (deactivate
// scripts don't read it).
func writeSourceLine(b *strings.Builder, shell Shell, file string) {
	if shell == ShellFish {
		fmt.Fprintf(b, "source %s\n", shellQuote(file))
		return
	}
	fmt.Fprintf(b, ". %s\n", shellQuote(file))
}

// writeScopedSourceLine sources file with CONDA_PREFIX=root visible only to
// that one call, matching what a real activate.d hook (e.g. libxml2's,
// which reads CONDA_PREFIX for XML_CATALOG_FILES) expects to see. bash/zsh
// get this for free from the "VAR=val cmd" prefix-assignment form; fish has
// no equivalent, so a begin/end block scopes a local "set -lx" instead.
func writeScopedSourceLine(b *strings.Builder, shell Shell, root, file string) {
	if shell == ShellFish {
		fmt.Fprintf(b, "begin; set -lx CONDA_PREFIX %s; source %s; end\n", shellQuote(root), shellQuote(file))
		return
	}
	fmt.Fprintf(b, "CONDA_PREFIX=%s . %s\n", shellQuote(root), shellQuote(file))
}

// sortedHookScripts lists dir's files ending in ext, sorted by name
// (reverse when asked). A missing dir is the common case — most packages
// ship neither hook — so it yields nil rather than an error.
func sortedHookScripts(dir, ext string, reverse bool) []string {
	entries, err := os.ReadDir(dir)
	if err != nil {
		return nil
	}
	var scripts []string
	for _, entry := range entries {
		if !entry.IsDir() && strings.HasSuffix(entry.Name(), ext) {
			scripts = append(scripts, filepath.Join(dir, entry.Name()))
		}
	}
	sort.Strings(scripts)
	if reverse {
		sort.Sort(sort.Reverse(sort.StringSlice(scripts)))
	}
	return scripts
}

// shellQuote renders a path as a single-quoted shell word. The \' escape
// works in both POSIX shells (close-quote, escaped literal quote, reopen)
// and fish (its \' escape inside single quotes), so one implementation
// covers every Shell.
func shellQuote(s string) string {
	return "'" + strings.ReplaceAll(s, "'", `'\''`) + "'"
}
