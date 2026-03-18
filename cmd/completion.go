package cmd

import (
	"bytes"
	"fmt"
	"os"
	"strings"

	"github.com/spf13/cobra"
)

// detectLoginShell reads $SHELL (the user's configured login shell).
func detectLoginShell() string {
	shell := os.Getenv("SHELL")
	shellLower := strings.ToLower(shell)

	// Check for specific shells
	if strings.Contains(shellLower, "fish") {
		return "fish"
	}
	if strings.Contains(shellLower, "zsh") {
		return "zsh"
	}

	// Default to bash
	return "bash"
}

// detectCompletionShell detects which shell is requesting completions by
// inspecting the parent process via /proc (Linux/HPC). Falls back to $SHELL.
func detectCompletionShell() string {
	ppid := os.Getppid()
	if exe, err := os.Readlink(fmt.Sprintf("/proc/%d/exe", ppid)); err == nil {
		exe = strings.ToLower(exe)
		switch {
		case strings.Contains(exe, "fish"):
			return "fish"
		case strings.Contains(exe, "zsh"):
			return "zsh"
		case strings.Contains(exe, "bash"):
			return "bash"
		}
	}
	return detectLoginShell() // fallback to $SHELL on non-Linux or unrecognized parent
}

var completionCmd = &cobra.Command{
	Use:   "completion [bash|zsh|fish]",
	Short: "Generate shell completion script",
	Long: `Generate shell completion script for condatainer.

If no shell is specified, the current shell will be auto-detected.`,
	Example: `  Bash:
  # Current session
  source <(condatainer completion bash)
  # All sessions (install once)
  condatainer completion bash > /etc/bash_completion.d/condatainer

  Zsh:
  # Current session
  source <(condatainer completion zsh)
  # All sessions (install once)
  condatainer completion zsh > "${fpath[1]}/_condatainer"
  # Note: If compinit is not enabled, add to ~/.zshrc:
  # autoload -U compinit; compinit

  Fish:
  # Current session
  condatainer completion fish | source
  # All sessions (install once)
  condatainer completion fish > ~/.config/fish/completions/condatainer.fish`,
	DisableFlagsInUseLine: true,
	ValidArgs:             []string{"bash", "zsh", "fish"},
	Args:                  cobra.MatchAll(cobra.MaximumNArgs(1), cobra.OnlyValidArgs),
	Run: func(cmd *cobra.Command, args []string) {
		shell := detectCompletionShell()
		if len(args) > 0 {
			shell = args[0]
		}

		switch shell {
		case "bash":
			// Generate to buffer so we can post-process
			var buf bytes.Buffer
			if err := cmd.Root().GenBashCompletionV2(&buf, true); err != nil {
				_ = cmd.Root().GenBashCompletion(&buf)
			}
			// Post-process to handle -- and add overlay fzf support
			script := postProcessBashCompletion(buf.String())
			os.Stdout.WriteString(script)
		case "zsh":
			// Generate to buffer so we can post-process
			var buf bytes.Buffer
			cmd.Root().GenZshCompletion(&buf)
			script := postProcessZshCompletion(buf.String())
			os.Stdout.WriteString(script)
		case "fish":
			// Generate to buffer so we can post-process
			var buf bytes.Buffer
			cmd.Root().GenFishCompletion(&buf, true)
			script := postProcessFishCompletion(buf.String())
			os.Stdout.WriteString(script)
		}
	},
}

func init() {
	rootCmd.AddCommand(completionCmd)
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

	// 2. Add fzf support if available for 'e/exec/instance start'
	fzfInject := `    __condatainer_debug "The completions are: ${out}"

    # Use fzf if available and we're completing for 'e/exec/instance start'
    if command -v fzf >/dev/null 2>&1 && [[ -n "$out" ]]; then
        local is_overlay_cmd=false
        for i in "${!words[@]}"; do
            word="${words[$i]}"
            if [[ "$word" == "e" || "$word" == "exec" ]]; then
                is_overlay_cmd=true
                break
            elif [[ "$word" == "instance" && "${words[$((i+1))]}" == "start" ]]; then
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
	// Add fzf support if available for 'e/exec/instance start'
	fzfInject := `    __condatainer_debug "completions: ${out}"

    # Use fzf if available and we're completing for 'e/exec/instance start'
    if command -v fzf >/dev/null 2>&1 && [[ -n "$out" ]]; then
        local is_overlay_cmd=false
        local _prev="" _word
        for _word in "${words[@]}"; do
            if [[ "$_word" == "e" || "$_word" == "exec" ]]; then
                is_overlay_cmd=true
                break
            elif [[ "$_prev" == "instance" && "$_word" == "start" ]]; then
                is_overlay_cmd=true
                break
            fi
            _prev="$_word"
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
	// Add fzf support if available for 'e/exec/instance start'
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

    # Use fzf if available and we're completing for 'e/exec/instance start'
    if type -q fzf
        set -l _is_overlay_cmd false
        if contains "e" $words; or contains "exec" $words
            set _is_overlay_cmd true
        else if contains "instance" $words; and contains "start" $words
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
