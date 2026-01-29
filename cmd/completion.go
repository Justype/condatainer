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
			// Post-process to handle -- for file completion
			script := postProcessBashCompletion(buf.String())
			os.Stdout.WriteString(script)
		case "zsh":
			cmd.Root().GenZshCompletion(os.Stdout)
		case "fish":
			cmd.Root().GenFishCompletion(os.Stdout, true)
		}
	},
}

func init() {
	rootCmd.AddCommand(completionCmd)
}

// postProcessBashCompletion modifies the generated bash completion script
// to handle -- properly (use file completion after --)
func postProcessBashCompletion(script string) string {
	// Find the __condatainer_get_completion_results function and add -- handling
	// We inject code to check if -- is in the words array, and if so, use default file completion
	oldCode := `args=("${words[@]:1}")
    requestComp="${words[0]} __complete ${args[*]}"`

	newCode := `args=("${words[@]:1}")
    # Check if -- is in the command line; if so, use default file completion
    for word in "${words[@]}"; do
        if [[ "$word" == "--" ]]; then
            return
        fi
    done
    requestComp="${words[0]} __complete ${args[*]}"`

	return strings.Replace(script, oldCode, newCode, 1)
}
