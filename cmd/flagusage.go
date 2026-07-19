package cmd

import (
	"strings"

	"github.com/spf13/pflag"
)

// flagUsages renders a FlagSet's help like pflag's FlagUsages, but for sections
// where no flag has a shorthand it reclaims the 4 extra columns pflag reserves
// for the "-x, " shorthand slot — giving a 2-space indent instead of 6 while
// keeping the aligned description column.
//
// Sections that do contain shorthands (e.g. the general "Flags:" block) are
// returned unchanged, so "-n, --name" style alignment is preserved.
func flagUsages(fs *pflag.FlagSet) string {
	hasShorthand := false
	fs.VisitAll(func(f *pflag.Flag) {
		if f.Shorthand != "" && f.ShorthandDeprecated == "" {
			hasShorthand = true
		}
	})

	usage := fs.FlagUsages()
	if hasShorthand {
		return usage
	}

	// Long-only section: every line is indented by pflag's fixed 6 spaces
	// (2 base + 4 reserved shorthand columns). Drop 4 to land at a 2-space indent.
	// Shifting every line equally preserves the description alignment.
	lines := strings.Split(usage, "\n")
	for i, ln := range lines {
		lines[i] = strings.TrimPrefix(ln, "    ")
	}
	return strings.Join(lines, "\n")
}
