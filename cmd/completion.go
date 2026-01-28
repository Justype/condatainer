package cmd

import (
	"bytes"
	"os"
	"strings"

	"github.com/spf13/cobra"
	"github.com/spf13/pflag"
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
	if strings.Contains(shellLower, "pwsh") || strings.Contains(shellLower, "powershell") {
		return "powershell"
	}

	// Default to bash
	return "bash"
}

var completionCmd = &cobra.Command{
	Use:   "completion [bash|zsh|fish|powershell]",
	Short: "Generate shell completion script",
	Long: func() string {
		detected := detectShell()
		return `Generate shell completion script for condatainer.

If no shell is specified, ` + detected + ` will be used (auto-detected from $SHELL).

To load completions:

Bash:
  $ source <(condatainer completion bash)

  # To load completions for each session, execute once:
  # Linux:
  $ condatainer completion bash > /etc/bash_completion.d/condatainer
  # macOS:
  $ condatainer completion bash > $(brew --prefix)/etc/bash_completion.d/condatainer

Zsh:
  # If shell completion is not already enabled in your environment,
  # you will need to enable it.  You can execute the following once:
  $ echo "autoload -U compinit; compinit" >> ~/.zshrc

  # To load completions for each session, execute once:
  $ condatainer completion zsh > "${fpath[1]}/_condatainer"

  # You will need to start a new shell for this setup to take effect.

Fish:
  $ condatainer completion fish | source

  # To load completions for each session, execute once:
  $ condatainer completion fish > ~/.config/fish/completions/condatainer.fish

PowerShell:
  PS> condatainer completion powershell | Out-String | Invoke-Expression

  # To load completions for every new session, run:
  PS> condatainer completion powershell > condatainer.ps1
  # and source this file from your PowerShell profile.
`
	}(),
	DisableFlagsInUseLine: true,
	ValidArgs:             []string{"bash", "zsh", "fish", "powershell"},
	Args:                  cobra.MatchAll(cobra.MaximumNArgs(1), cobra.OnlyValidArgs),
	Run: func(cmd *cobra.Command, args []string) {
		shell := detectShell()
		if len(args) > 0 {
			shell = args[0]
		}

		// Temporarily strip short shorthands (-x) from flags so completion shows
		// only long options (e.g., --long-arg). We restore them after generation.
		saved := stripShortFlagShorthands(cmd.Root())
		defer restoreShortFlagShorthands(cmd.Root(), saved)

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
		case "powershell":
			cmd.Root().GenPowerShellCompletionWithDesc(os.Stdout)
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

// stripShortFlagShorthands walks the command tree and clears the Shorthand
// field for any flag that has one, returning a map of saved values so they
// can be restored later.
func stripShortFlagShorthands(root *cobra.Command) map[string]string {
	saved := make(map[string]string)

	// Helper to strip shorthand from a flag and save it
	stripFlag := func(f *pflag.Flag) {
		if f.Shorthand != "" {
			saved[f.Name] = f.Shorthand
			f.Shorthand = ""
		}
	}

	var walk func(c *cobra.Command)
	walk = func(c *cobra.Command) {
		// Strip from local flags
		c.LocalFlags().VisitAll(stripFlag)
		// Strip from persistent flags defined on this command
		c.PersistentFlags().VisitAll(stripFlag)
		// Strip from inherited flags (persistent flags from parent commands)
		c.InheritedFlags().VisitAll(stripFlag)

		for _, child := range c.Commands() {
			walk(child)
		}
	}
	walk(root)
	return saved
}

// restoreShortFlagShorthands restores previously-saved shorthand values.
func restoreShortFlagShorthands(root *cobra.Command, saved map[string]string) {
	// Helper to restore shorthand for a flag
	restoreFlag := func(f *pflag.Flag) {
		if old, ok := saved[f.Name]; ok {
			f.Shorthand = old
		}
	}

	var walk func(c *cobra.Command)
	walk = func(c *cobra.Command) {
		c.LocalFlags().VisitAll(restoreFlag)
		c.PersistentFlags().VisitAll(restoreFlag)
		c.InheritedFlags().VisitAll(restoreFlag)

		for _, child := range c.Commands() {
			walk(child)
		}
	}
	walk(root)
}
