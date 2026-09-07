package container

import (
	"fmt"
	"strings"

	"github.com/Justype/condatainer/catalog"
)

// ActivationScript returns a bash preamble that sources every mounted app
// overlay's own etc/conda/activate.d/*.sh, and /cnt_env's when a conda
// env/.img is mounted (lastImg != ""), in overlay order — the same order
// BuildPathEnv iterates before its own prepend reverses it, so the
// last-listed overlay's hook still runs last and so still wins a name
// collision, consistent with how PATH already resolves one. Each block scopes
// CONDA_PREFIX to that overlay's own prefix while its scripts run: an
// activate.d script reads CONDA_PREFIX to build its own exports (e.g.
// libxml2's own XML_CATALOG_FILES), so sourcing it under some other overlay's
// CONDA_PREFIX would compute the wrong value.
//
// This replays only conda's own "first activation" branch (activate.py's
// build_activate): a container run mounts, executes once, and exits — never
// nested activate/deactivate within one invocation — so there is no
// CONDA_SHLVL stack, no deactivate.d, and nothing to restore.
//
// Returns "" when no mounted overlay could contribute an activate.d
// directory, so a caller can skip wrapping the command at all in that case.
func ActivationScript(overlays []string, lastImg string) string {
	var b strings.Builder
	for _, ov := range overlays {
		contribution, _ := resolveImage(cleanOverlayPath(ov))
		if contribution.Type != catalog.TypeApp || contribution.Prefix == "" {
			continue
		}
		writeActivateBlock(&b, contribution.Prefix)
	}
	if lastImg != "" {
		writeActivateBlock(&b, "/cnt_env")
	}
	return b.String()
}

func writeActivateBlock(b *strings.Builder, prefix string) {
	q := shellQuote(prefix)
	fmt.Fprintf(b, "if [ -d %s/etc/conda/activate.d ]; then\n", q)
	fmt.Fprintf(b, "  CONDA_PREFIX=%s\n", q)
	b.WriteString("  export CONDA_PREFIX\n")
	fmt.Fprintf(b, "  for __cnt_f in %s/etc/conda/activate.d/*.sh; do\n", q)
	b.WriteString("    [ -e \"$__cnt_f\" ] || continue\n")
	b.WriteString("    . \"$__cnt_f\"\n")
	b.WriteString("  done\n")
	b.WriteString("fi\n")
}

// shellQuote renders a path as a single-quoted shell word.
func shellQuote(s string) string {
	return "'" + strings.ReplaceAll(s, "'", `'\''`) + "'"
}
