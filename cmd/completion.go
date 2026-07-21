package cmd

import (
	"bytes"
	"fmt"
	"os"
	"path/filepath"
	"regexp"
	"strings"

	"github.com/spf13/cobra"
)

// completionAlias is an extra command name to register completion for (--alias).
var completionAlias string

// aliasNamePattern restricts --alias to a plain shell word: the value is
// interpolated into the generated script, so spaces or quotes would silently
// produce a broken completion file.
var aliasNamePattern = regexp.MustCompile(`^[A-Za-z0-9_-]+$`)

// shellFromPath maps an executable path to a supported shell name, matching on
// the base name so a path like /home/zsh-fan/bin/foo does not match.
// Returns "" for anything else.
func shellFromPath(path string) string {
	base := strings.ToLower(filepath.Base(path))
	switch {
	case strings.Contains(base, "fish"):
		return "fish"
	case strings.Contains(base, "zsh"):
		return "zsh"
	case strings.Contains(base, "bash"):
		return "bash"
	}
	return ""
}

// detectCompletionShell detects which shell is requesting completions from the
// parent process via /proc, falling back to $SHELL only when /proc is unreadable.
//
// An unrecognized parent is reported as unsupported rather than guessed at: the
// output is meant to be sourced, so emitting the wrong shell's script turns a
// clear error into a screenful of syntax errors in the user's login.
func detectCompletionShell() (string, error) {
	if exe, err := os.Readlink(fmt.Sprintf("/proc/%d/exe", os.Getppid())); err == nil {
		if shell := shellFromPath(exe); shell != "" {
			return shell, nil
		}
		return "", fmt.Errorf("no supported shell detected (parent process: %s)", filepath.Base(exe))
	}
	if shell := shellFromPath(os.Getenv("SHELL")); shell != "" {
		return shell, nil
	}
	return "", fmt.Errorf("could not detect the current shell")
}

var completionCmd = &cobra.Command{
	Use:   "completion [bash|zsh|fish]",
	Short: "Generate shell completion script",
	Long: `Generate shell completion script for condatainer.

If no shell is specified, the current shell will be auto-detected.
Use --alias to define a shorter name, e.g. --alias cnt lets you type 'cnt avail'.
An alias only takes effect when the script is sourced from shell rc file.`,
	Example: `  Bash:
  # Current session
  source <(condatainer completion bash)

  Zsh:
  # Current session
  source <(condatainer completion zsh)
  # All sessions (install once)
  condatainer completion zsh > "${fpath[1]}/_condatainer"

  Fish:
  # Current session
  condatainer completion fish | source
  # All sessions (install once)
  condatainer completion fish > ~/.config/fish/completions/condatainer.fish`,
	DisableFlagsInUseLine: true,
	ValidArgs:             []string{"bash", "zsh", "fish"},
	Args:                  cobra.MatchAll(cobra.MaximumNArgs(1), cobra.OnlyValidArgs),
	Run: func(cmd *cobra.Command, args []string) {
		var shell string
		if len(args) > 0 {
			shell = args[0]
		} else {
			detected, err := detectCompletionShell()
			if err != nil {
				ExitWithError("%v.\nSpecify it explicitly: condatainer completion bash|zsh|fish", err)
			}
			shell = detected
		}

		if completionAlias != "" && !aliasNamePattern.MatchString(completionAlias) {
			ExitWithError("Invalid --alias %q: use letters, digits, '-' or '_' only.", completionAlias)
		}
		name := cmd.Root().Name()

		switch shell {
		case "bash":
			// Generate to buffer so we can post-process
			var buf bytes.Buffer
			if err := cmd.Root().GenBashCompletionV2(&buf, true); err != nil {
				_ = cmd.Root().GenBashCompletion(&buf)
			}
			// Post-process to handle -- and add overlay fzf support
			script := postProcessBashCompletion(buf.String())
			os.Stdout.WriteString(script + bashAliasBinding(name, completionAlias))
		case "zsh":
			// Generate to buffer so we can post-process
			var buf bytes.Buffer
			cmd.Root().GenZshCompletion(&buf)
			script := postProcessZshCompletion(buf.String())
			os.Stdout.WriteString(script + zshAliasBinding(name, completionAlias))
		case "fish":
			// Generate to buffer so we can post-process
			var buf bytes.Buffer
			cmd.Root().GenFishCompletion(&buf, true)
			script := postProcessFishCompletion(buf.String())
			os.Stdout.WriteString(script + fishAliasBinding(script, name, completionAlias))
		}
	},
}

func init() {
	rootCmd.AddCommand(completionCmd)
	completionCmd.Flags().StringVar(&completionAlias, "alias", "",
		"Define a short alias (e.g. cnt) and complete it too")
}

// bashAliasBinding defines the alias and binds it to the generated completion
// function. Empty when no alias was requested.
func bashAliasBinding(name, alias string) string {
	if alias == "" {
		return ""
	}
	return fmt.Sprintf(`
alias %[2]s=%[1]s
if [[ $(type -t compopt) = "builtin" ]]; then
    complete -o default -F __start_%[1]s %[2]s
else
    complete -o default -o nospace -F __start_%[1]s %[2]s
fi
`, name, alias)
}

// zshAliasBinding defines the alias and binds it to the generated completion
// function. Empty when no alias was requested.
func zshAliasBinding(name, alias string) string {
	if alias == "" {
		return ""
	}
	return fmt.Sprintf("\nalias %[2]s=%[1]s\ncompdef _%[1]s %[2]s\n", name, alias)
}

// fishAliasBinding mirrors the script's 'complete -c <name>' lines onto alias.
// Fish binds each command separately, so the lines are duplicated verbatim.
// Empty when no alias was requested.
func fishAliasBinding(script, name, alias string) string {
	if alias == "" {
		return ""
	}
	var b strings.Builder
	fmt.Fprintf(&b, "\nalias %s %s\n", alias, name)
	prefix := "complete -c " + name + " "
	for _, line := range strings.Split(script, "\n") {
		if strings.HasPrefix(line, prefix) {
			fmt.Fprintf(&b, "complete -c %s %s\n", alias, strings.TrimPrefix(line, prefix))
		}
	}
	return b.String()
}

// bashCompletionFallback defines minimal stubs for _get_comp_words_by_ref and
// _init_completion so that the generated script works on systems that do not
// have the bash-completion package installed.
const bashCompletionFallback = `
# Minimal fallback for systems without the bash-completion package
if ! type _get_comp_words_by_ref &>/dev/null; then
    _get_comp_words_by_ref() {
        while [[ $# -gt 0 ]]; do
            case $1 in
                -n) shift 2 ;;
                cur)   cur=${COMP_WORDS[COMP_CWORD]} ;;
                prev)  prev=${COMP_WORDS[COMP_CWORD-1]} ;;
                words) words=("${COMP_WORDS[@]}") ;;
                cword) cword=$COMP_CWORD ;;
            esac
            shift
        done
    }
fi
if ! type _init_completion &>/dev/null; then
    _init_completion() {
        COMPREPLY=()
        while [[ $1 == -* ]]; do shift; [[ $1 == -n ]] && shift; done
        cur=${COMP_WORDS[COMP_CWORD]}
        prev=${COMP_WORDS[COMP_CWORD-1]}
        words=("${COMP_WORDS[@]}")
        cword=$COMP_CWORD
        return 0
    }
fi
`

// postProcessBashCompletion modifies the generated bash completion script
func postProcessBashCompletion(script string) string {
	// 0. Prepend fallback stubs for systems without bash-completion installed
	script = bashCompletionFallback + script

	// 1. Handle -- properly (use file completion after --)
	oldCode := `args=("${words[@]:1}")
    requestComp="${words[0]} __complete ${args[*]}"`

	newCode := `args=("${words[@]:1}")
    for word in "${words[@]}"; do
        if [[ "$word" == "--" ]]; then
            return
        fi
    done
    requestComp="${words[0]} __complete ${args[*]}"`

	script = strings.Replace(script, oldCode, newCode, 1)

	// 2. Hide short flags, then add fzf support if available for 'e/exec'.
	// Both are injected at the same anchor so there is only one replacement to keep
	// in sync: the line where 'out' holds the candidates with the directive removed.
	fzfInject := `# Bash lists short and long flags on separate lines, doubling the candidate
    # list. Drop the short forms when completing a flag.
    if [[ -n ${out} && ${cur} == -* ]]; then
        out=$(printf '%s\n' "${out}" | grep -vE '^-[a-zA-Z0-9]([[:space:]]|$)')
    fi

    __condatainer_debug "The completions are: ${out}"

    # Use fzf if available and we're completing for 'e/exec'
    if command -v fzf >/dev/null 2>&1 && [[ -n "$out" ]]; then
        local is_overlay_cmd=false
        for i in "${!words[@]}"; do
            word="${words[$i]}"
            if [[ "$word" == "e" || "$word" == "exec" ]]; then
                is_overlay_cmd=true
                break
            fi
        done

        if [[ "$cur" == -* ]]; then
            is_overlay_cmd=false
        fi

        if [[ "$is_overlay_cmd" == "true" ]]; then
            if [[ $(echo "$out" | wc -l) -gt 1 ]]; then
                 local selection
                 selection=$(echo "$out" | fzf --reverse --header="Select overlay" --query="$cur" --select-1 --exit-0)
                 if [[ -n "$selection" ]]; then
                     out="$selection"
                 fi
            fi
        fi
    fi`

	script = strings.Replace(script, `__condatainer_debug "The completions are: ${out}"`, fzfInject, 1)

	return script
}

// postProcessZshCompletion modifies the generated zsh completion script
// EXPERIMENTAL: zsh completion is not tested and may have edge cases. Feedback welcome.
func postProcessZshCompletion(script string) string {
	// Add fzf support if available for 'e/exec'
	fzfInject := `    __condatainer_debug "completions: ${out}"

    # Use fzf if available and we're completing for 'e/exec'
    if command -v fzf >/dev/null 2>&1 && [[ -n "$out" ]]; then
        local is_overlay_cmd=false
        local _word
        for _word in "${words[@]}"; do
            if [[ "$_word" == "e" || "$_word" == "exec" ]]; then
                is_overlay_cmd=true
                break
            fi
        done

        if [[ "${words[$CURRENT]}" == -* ]]; then
            is_overlay_cmd=false
        fi

        if [[ "$is_overlay_cmd" == "true" ]]; then
             local -a lines
             lines=("${(@f)out}")
             if [[ ${#lines} -gt 1 ]]; then
                 zle -I
                 local selection
                 selection=$(echo "$out" | fzf --height 40% --reverse --header="Select overlay" --query="$lastParam" --select-1 --exit-0)
                 if [[ -n "$selection" ]]; then
                     out="$selection"
                 fi
             fi
        fi
    fi`

	script = strings.Replace(script, `__condatainer_debug "completions: ${out}"`, fzfInject, 1)

	return script
}

// postProcessFishCompletion modifies the generated fish completion script
// EXPERIMENTAL: fish completion is not tested and may have edge cases. Feedback welcome.
func postProcessFishCompletion(script string) string {
	// Add fzf support if available for 'e/exec'
	// Note: We try to match the line where results are captured.
	// Cobra sometimes generates 'eval' and sometimes not depending on version.
	// We handle the standard pattern found in recent versions.
	// Match cobra v1.10.x generated format (eval, no extra redirect)
	target := `set -l results (eval $requestComp 2> /dev/null)`
	if !strings.Contains(script, target) {
		// Older cobra without 'eval'
		target = `set -l results ($requestComp 2> /dev/null)`
	}

	fzfInject := target + `

    # Use fzf if available and we're completing for 'e/exec'
    if type -q fzf
        set -l _is_overlay_cmd false
        if contains "e" $words; or contains "exec" $words
            set _is_overlay_cmd true
        end
        if test $_is_overlay_cmd = true
            if not string match -q -- "-*" (commandline -t)
                # Check if -- is in the command line; if so, skip fzf and use default results
                if not contains -- "--" $words
                    set -l candidates $results[1..-2]
                    if test (count $candidates) -gt 1
                        set -l selection (string join \n $candidates | fzf --height 40% --reverse --select-1 --exit-0 --query (commandline -t) --header "Select overlay")
                        if test -n "$selection"
                            set results $selection $results[-1]
                        end
                    end
                end
            end
        end
    end`

	return strings.Replace(script, target, fzfInject, 1)
}
