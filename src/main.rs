use anyhow::Result;
use clap::{CommandFactory, Parser, Subcommand};
use clap_complete;

mod apptainer;
mod build;
mod commands;
mod config;
mod overlay;
mod scheduler;
mod utils;

use commands::avail::{handle_avail, AvailArgs};
use commands::config::{handle_config, ConfigArgs};
use commands::create::{handle_create, CreateArgs};
use commands::e::{handle_e, EArgs};
use commands::exec::{handle_exec, ExecArgs};
use commands::list::{handle_list, ListArgs};
use commands::o::{handle_o, OArgs};
use commands::overlay::{handle_overlay, OverlayArgs};
use commands::remove::{handle_remove, RemoveArgs};
use commands::scheduler::{handle_scheduler, SchedulerArgs};
use config::AppConfig;

fn main() -> Result<()> {
    let Cli {
        debug,
        local,
        command,
    } = Cli::parse();
    
    // Set debug mode in utils
    utils::set_debug_mode(debug);
    
    let config = AppConfig::new(debug, !local)?;

    if debug {
        utils::print_debug(&format!("CondaTainer Version: {}", config.version));
        utils::print_debug(&format!("Program base: {}", config.program_base.display()));
        utils::print_debug(&format!("Base directory: {}", config.base_dir.display()));
        utils::print_debug(&format!("Apptainer binary: {}", config.apptainer_bin.display()));
        println!("{config:#?}");
    }

    run_command(command, &config)
}

fn run_command(command: Command, config: &AppConfig) -> Result<()> {
    match command {
        Command::Avail(args) => handle_avail(&args, config),
        Command::Completion { shell } => {
            let detected_shell = if let Some(s) = shell {
                s
            } else {
                // Try to auto-detect from $SHELL
                match std::env::var("SHELL") {
                    Ok(shell_path) => {
                        let shell_name = std::path::Path::new(&shell_path)
                            .file_name()
                            .and_then(|n| n.to_str())
                            .unwrap_or("");
                        
                        match shell_name {
                            "bash" => clap_complete::Shell::Bash,
                            "zsh" => clap_complete::Shell::Zsh,
                            "fish" => clap_complete::Shell::Fish,
                            "elvish" => clap_complete::Shell::Elvish,
                            "powershell" | "pwsh" => clap_complete::Shell::PowerShell,
                            _ => {
                                utils::print_error(&format!("Could not auto-detect shell from SHELL={}. Supported shells: bash, zsh, fish, elvish, powershell", shell_path));
                                utils::print_message("\nUsage: condatainer completion <SHELL>");
                                utils::print_message("\nExample: condatainer completion bash");
                                return Ok(());
                            }
                        }
                    }
                    Err(_) => {
                        utils::print_error("Could not auto-detect shell. SHELL environment variable not set.");
                        utils::print_message("\nUsage: condatainer completion <SHELL>");
                        utils::print_message("\nSupported shells: bash, zsh, fish, elvish, powershell");
                        utils::print_message("Example: condatainer completion bash");
                        return Ok(());
                    }
                }
            };
            
            // Generate enhanced completion for bash, zsh, and fish
            match detected_shell {
                clap_complete::Shell::Bash => {
                    commands::completion::generate_enhanced_bash_completion();
                }
                clap_complete::Shell::Zsh => {
                    commands::completion::generate_enhanced_zsh_completion();
                }
                clap_complete::Shell::Fish => {
                    commands::completion::generate_enhanced_fish_completion();
                }
                _ => {
                    let mut cmd = Cli::command();
                    let bin_name = cmd.get_name().to_string();
                    clap_complete::generate(detected_shell, &mut cmd, bin_name, &mut std::io::stdout());
                }
            }
            Ok(())
        }
        Command::CompleteOverlay { include_data, to_complete } => {
            let include_data_bool = include_data == "true";
            let suggestions = commands::completion::get_overlay_suggestions(include_data_bool, &to_complete);
            for suggestion in suggestions {
                println!("{}", suggestion);
            }
            Ok(())
        }
        Command::CompleteBaseImage { to_complete } => {
            let suggestions = commands::completion::get_base_image_suggestions(&to_complete);
            for suggestion in suggestions {
                println!("{}", suggestion);
            }
            Ok(())
        }
        Command::Config(args) => handle_config(&args, config),
        Command::Create(args) => handle_create(args, config),
        Command::E(args) => handle_e(&args, config),
        Command::Exec(args) => handle_exec(&args, config),
        Command::List(args) => handle_list(&args, config),
        Command::O(args) => handle_o(&args, config),
        Command::Overlay(args) => handle_overlay(&args, config),
        Command::Remove(args) => handle_remove(&args, config),
        Command::Scheduler(args) => handle_scheduler(&args, config),
    }
}

#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "Use Apptainer/Conda/SquashFS to manage tools for HPC users.",
    long_about = "For the full manual, see https://github.com/Justype/condatainer/blob/main/docs/manuals/condatainer.md"
)]
struct Cli {
    #[arg(short, long, help = "Enable debug mode with verbose output")]
    debug: bool,

    #[arg(
        long,
        help = "Skip submitting jobs to a scheduler and keep everything local"
    )]
    local: bool,

    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    #[command(about = "Check available build scripts", alias = "av")]
    Avail(AvailArgs),
    #[command(
        about = "Generate shell completion scripts",
        long_about = "Generate shell completion scripts for various shells with dynamic overlay completion.

Auto-detects your shell from $SHELL if no argument is provided.

Enhanced completion features:
  - Tab completion for overlay names (e and exec -o)
  - Tab completion for base images (e -b and exec -b)
  - Real-time filesystem queries (always up-to-date after create/remove)

Installation:
  Bash:
    source <(condatainer completion)
    # Or permanently: echo 'source <(condatainer completion)' >> ~/.bashrc

  Zsh:
    source <(condatainer completion)
    # Or permanently: echo 'source <(condatainer completion)' >> ~/.zshrc

  Fish:
    condatainer completion fish > ~/.config/fish/completions/condatainer.fish

  PowerShell:
    condatainer completion powershell | Out-String | Invoke-Expression

Examples:
  condatainer completion              # Auto-detect shell
  condatainer completion bash         # Generate bash completion
  condatainer completion zsh          # Generate zsh completion
  condatainer completion fish         # Generate fish completion"
    )]
    Completion {
        #[arg(value_enum, help = "Shell to generate completion for (auto-detected if omitted)")]
        shell: Option<clap_complete::Shell>,
    },
    #[command(hide = true, about = "Internal: complete overlay names")]
    CompleteOverlay {
        #[arg(default_value = "true")]
        include_data: String,
        #[arg(default_value = "")]
        to_complete: String,
    },
    #[command(hide = true, about = "Internal: complete base image names")]
    CompleteBaseImage {
        #[arg(default_value = "")]
        to_complete: String,
    },
    #[command(about = "Manage condatainer configuration")]
    Config(ConfigArgs),
    #[command(about = "Create a new SquashFS overlay", alias = "install", alias = "i")]
    Create(CreateArgs),
    #[command(about = "Run bash using writable overlays")]
    E(EArgs),
    #[command(about = "Execute a command using overlays")]
    Exec(ExecArgs),
    #[command(about = "List installed overlays", alias = "ls")]
    List(ListArgs),
    #[command(about = "Shortcut for 'overlay create'")]
    O(OArgs),
    #[command(about = "Manage persistent overlay images (create, resize, check, info)")]
    Overlay(OverlayArgs),
    #[command(about = "Remove installed overlays matching search terms", alias = "rm", alias = "delete", alias = "uninstall")]
    Remove(RemoveArgs),
    #[command(about = "Display scheduler information", alias = "sched")]
    Scheduler(SchedulerArgs),
}
