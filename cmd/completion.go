package cmd

import (
	"bytes"
	"os"
	"strings"

	"github.com/spf13/cobra"
)

// detectShell auto-detects the current shell from environment
func detectShell() string {
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

var completionCmd = &cobra.Command{
	Use:   "completion [bash|zsh|fish]",
	Short: "Generate shell completion script",
	Long: func() string {
		detected := detectShell()
		return `Generate shell completion script for condatainer.

If no shell is specified, ` + detected + ` will be used (auto-detected from $SHELL).

To load completions:

Bash:
  # Current session:
  $ source <(condatainer completion bash)

  # All sessions (install once):
  $ condatainer completion bash > /etc/bash_completion.d/condatainer

Zsh:
  # Current session:
  $ source <(condatainer completion zsh)

  # All sessions (install once):
  $ condatainer completion zsh > "${fpath[1]}/_condatainer"
  
  # Note: If compinit is not enabled, add to ~/.zshrc:
  # autoload -U compinit; compinit

Fish:
  # Current session:
  $ condatainer completion fish | source

  # All sessions (install once):
  $ condatainer completion fish > ~/.config/fish/completions/condatainer.fish
`
	}(),
	DisableFlagsInUseLine: true,
	ValidArgs:             []string{"bash", "zsh", "fish"},
	Args:                  cobra.MatchAll(cobra.MaximumNArgs(1), cobra.OnlyValidArgs),
	Run: func(cmd *cobra.Command, args []string) {
		shell := detectShell()
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

// postProcessBashCompletion modifies the generated bash completion script
func postProcessBashCompletion(script string) string {
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

	// 2. Add fzf support if available for 'e' and 'exec'
	fzfInject := `    __condatainer_debug "The completions are: ${out}"

    # Use fzf if available and we're completing for 'e' or 'exec'
    if command -v fzf >/dev/null 2>&1 && [[ -n "$out" ]]; then
        local is_overlay_cmd=false
        for word in "${words[@]}"; do
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
	// Add fzf support if available for 'e' and 'exec'
	fzfInject := `    __condatainer_debug "completions: ${out}"

    # Use fzf if available and we're completing for 'e' or 'exec'
    if command -v fzf >/dev/null 2>&1 && [[ -n "$out" ]]; then
        local is_overlay_cmd=false
        for word in "${words[@]}"; do
            if [[ "$word" == "e" || "$word" == "exec" ]]; then
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
	// Add fzf support if available for 'e' and 'exec'
	// Note: We try to match the line where results are captured.
	// Cobra sometimes generates 'eval' and sometimes not depending on version.
	// We handle the standard pattern found in recent versions.
	target := `set -l results ($requestComp 2> /dev/null)`
	// Fallback if the generator uses eval
	if !strings.Contains(script, target) {
		target = `set -l results (eval $requestComp 2> /dev/null)`
	}

	fzfInject := target + `

    # Use fzf if available and we're completing for 'e' or 'exec'
    if type -q fzf
        if contains "e" $words; or contains "exec" $words
            if not string match -q -- "-*" (commandline -t)
                # Check if -- is in the command line; if so, skip fzf and use default results
                if not contains -- "--" $words
                    set -l candidates $results[1..-2]
                    if test (count $candidates) -gt 1
                        set -l selection (string join \n $candidates | fzf --height 40% --reverse --select-1 --exit-0 --query (commandline -t) < /dev/tty)
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
