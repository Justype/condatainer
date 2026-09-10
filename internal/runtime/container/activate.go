package container

import (
	"fmt"
	"strings"
)

// ActivationMode selects which activate.d scripts ActivationScript sources.
type ActivationMode string

const (
	// ActivationAll sources every mounted app overlay's own activate.d plus
	// the mounted conda environment's.
	ActivationAll ActivationMode = "all"
	// ActivationEnv sources only the mounted conda environment's activate.d,
	// skipping every app overlay's own.
	ActivationEnv ActivationMode = "env"
	// ActivationNone sources neither.
	ActivationNone ActivationMode = "none"
)

// ActivationScript returns a bash preamble sourcing activate.d scripts for
// mode: ActivationAll sources every ContributesBin() overlay (each app,
// plus the mounted conda environment's own) deduplicated by Prefix, so a
// writable .img and its paired env.sqf snapshot are sourced once, not
// twice; ActivationEnv sources only /cnt_env; ActivationNone sources
// nothing. Returns "" when there is nothing to source, so a caller can
// skip wrapping the command at all.
func ActivationScript(overlays []string, envMounted bool, mode ActivationMode) string {
	var b strings.Builder
	switch mode {
	case ActivationAll:
		seen := map[string]bool{}
		for _, ov := range overlays {
			contribution, _ := resolveImage(cleanOverlayPath(ov))
			if !contribution.ContributesBin() || seen[contribution.Prefix] {
				continue
			}
			seen[contribution.Prefix] = true
			writeActivateBlock(&b, contribution.Prefix)
		}
	case ActivationEnv:
		if envMounted {
			writeActivateBlock(&b, EnvPrefix)
		}
	}
	return b.String()
}

// MMHelperScript returns a bash snippet defining and exporting an mm shell
// function that runs `condatainer env "$@"` and, after install/update/
// remove, re-sources activate.d/deactivate.d via `condatainer env
// reactivate --shell bash` (hardcoded, not detected: mm only ever runs
// inside bash, regardless of the ambient $SHELL reactivate would otherwise
// detect). export -f is what lets the function survive the subsequent
// exec into an interactive bash shell; it does not survive exec into zsh
// or fish. Returns "" when no environment is mounted.
func MMHelperScript(envMounted bool) string {
	if !envMounted {
		return ""
	}
	return `mm() {
    condatainer env "$@" || return
    case "$1" in
        install|update|remove) eval "$(condatainer env reactivate --shell bash)" ;;
    esac
}
export -f mm
`
}

func writeActivateBlock(b *strings.Builder, prefix string) {
	q := shellQuote(prefix)
	fmt.Fprintf(b, "if [ -d %s/etc/conda/activate.d ]; then\n", q)
	fmt.Fprintf(b, "  for __cnt_f in %s/etc/conda/activate.d/*.sh; do\n", q)
	b.WriteString("    [ -e \"$__cnt_f\" ] || continue\n")
	fmt.Fprintf(b, "    CONDA_PREFIX=%s . \"$__cnt_f\"\n", q)
	b.WriteString("  done\n")
	b.WriteString("fi\n")
}

// shellQuote renders a path as a single-quoted shell word.
func shellQuote(s string) string {
	return "'" + strings.ReplaceAll(s, "'", `'\''`) + "'"
}
