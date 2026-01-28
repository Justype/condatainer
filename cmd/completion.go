package cmd

import (
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
		switch shell {
		case "bash":
			cmd.Root().GenBashCompletion(os.Stdout)
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
