package helper

import (
	"fmt"
	"strings"

	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
)

// ApplyParamFlags scans raw CLI args for flags defined by params.
// Standard resource flags (-c/-m/-t/-g/-b/-e/-o/-w) are already parsed by cobra
// and must NOT be in args — this handles only helper-specific #PARAM: flags.
// Returns (values map, unrecognised args, error).
func ApplyParamFlags(params []HelperParam, args []string) (map[string]string, []string, error) {
	values := make(map[string]string)
	var extra []string

	// Build lookup maps: long flag → param, short flag → param
	byLong := make(map[string]*HelperParam, len(params))
	byShort := make(map[string]*HelperParam, len(params))
	for i := range params {
		p := &params[i]
		if p.LongFlag != "" {
			byLong[p.LongFlag] = p
		}
		if p.ShortFlag != "" {
			byShort[p.ShortFlag] = p
		}
	}

	i := 0
	for i < len(args) {
		arg := args[i]
		if p := byLong[arg]; p != nil {
			i++
			if i >= len(args) {
				return nil, nil, fmt.Errorf("flag %s requires a value", arg)
			}
			values[p.Key] = args[i]
			i++
			continue
		}
		if p := byShort[arg]; p != nil {
			i++
			if i >= len(args) {
				return nil, nil, fmt.Errorf("flag %s requires a value", arg)
			}
			values[p.Key] = args[i]
			i++
			continue
		}
		// Handle --flag=value form
		if strings.HasPrefix(arg, "--") {
			flagName, val, hasVal := strings.Cut(arg, "=")
			if p := byLong[flagName]; p != nil {
				if !hasVal {
					return nil, nil, fmt.Errorf("flag %s requires a value", flagName)
				}
				values[p.Key] = val
				i++
				continue
			}
		}
		extra = append(extra, arg)
		i++
	}
	return values, extra, nil
}

// ReservedShortFlags is the set of short flags consumed by parsePostScriptHelperFlags.
// -l/-u are cobra-level only and never reach script arg parsing, so they are excluded.
var ReservedShortFlags = map[string]string{
	"-c": "--cpus", "-m": "--mem", "-t": "--time", "-g": "--gpu",
	"-b": "--base", "-e": "--env", "-o": "--overlay", "-w": "--cwd",
}

// ValidateHelperParamConflicts returns an error if any #PARAM: short flag clashes
// with a short flag reserved by condatainer's script-level flag parser.
func ValidateHelperParamConflicts(params []HelperParam) error {
	for _, p := range params {
		if p.ShortFlag != "" {
			if builtin, conflict := ReservedShortFlags[p.ShortFlag]; conflict {
				return fmt.Errorf(
					"helper param %q uses short flag %s which conflicts with built-in %s\n"+
						"  Fix: change the short flag in the #PARAM: header to a different letter",
					p.Key, p.ShortFlag, builtin)
			}
		}
	}
	return nil
}

// PrintHelperUsage prints the combined --help text for a helper.
// versions is the map from ParseHelperScriptMeta().ParamValues; pass nil if unavailable.
func PrintHelperUsage(scriptName string, params []HelperParam, spec *scheduler.ResourceSpec, versions map[string][]string) {
	cpus := 0
	mem := ""
	timeStr := ""
	if spec != nil {
		cpus = spec.CpusPerTask
		if mb := spec.GetMemPerNodeMB(); mb > 0 {
			mem = utils.FormatMemoryMB(mb)
		}
		if spec.Time > 0 {
			timeStr = utils.FormatDuration(spec.Time)
		}
	}
	fmt.Printf("Usage: condatainer helper %s [flags]\n\n", scriptName)
	fmt.Println("Resource flags:")
	fmt.Printf("  -c, --cpus        CPUs per job   (default: %d)\n", cpus)
	fmt.Printf("  -m, --mem         Memory         (default: %s)\n", mem)
	fmt.Printf("  -t, --time        Walltime       (default: %s)\n", timeStr)
	fmt.Println("  -g, --gpu         GPU spec (e.g. a100:1)")
	fmt.Println("  -b, --base        Override base image")
	fmt.Println("  -e, --env         Writable overlay (.img)")
	fmt.Println("  -o, --overlay     Additional read-only overlay (repeatable)")
	fmt.Println("  -w, --cwd         Use current directory as working dir")
	fmt.Println("      --new         Skip reuse prompt, force new instance")
	if len(params) > 0 {
		fmt.Println("\nHelper-specific flags:")
		for _, p := range params {
			flags := ""
			if p.ShortFlag != "" && p.LongFlag != "" {
				flags = fmt.Sprintf("%s, %s", p.ShortFlag, p.LongFlag)
			} else if p.LongFlag != "" {
				flags = p.LongFlag
			} else if p.ShortFlag != "" {
				flags = p.ShortFlag
			}
			defStr := ""
			if p.Default != "" {
				defStr = fmt.Sprintf(" (default: %s)", p.Default)
			}
			choicesStr := ""
			if versions != nil {
				if vlist := versions[p.Key]; len(vlist) > 0 {
					choicesStr = "  " + utils.StyleDebug(utils.FormatChoicesInline(vlist))
				}
			}
			if flags != "" {
				fmt.Printf("  %-20s  %s%s%s\n", flags, p.Desc, defStr, choicesStr)
			} else {
				fmt.Printf("  %-20s  %s%s%s  [prompt-only]\n", p.Key, p.Desc, defStr, choicesStr)
			}
		}
	}
	fmt.Printf("\nSaved defaults: condatainer helper %s config [show | set KEY VALUE]\n", scriptName)
}
