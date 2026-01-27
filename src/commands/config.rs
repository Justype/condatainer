use anyhow::{Context, Result};
use clap::{Args, Subcommand};
use std::env;
use std::io::{self, Write};
use std::process::Command;

use crate::config::{
    self, detect_scheduler_bin, force_detect_and_save, parse_duration, validate_binary,
    AppConfig,
};
use crate::utils::{self, style_error, style_success, style_warning, style_hint, style_info};

#[derive(Args, Debug)]
pub struct ConfigArgs {
    #[command(subcommand)]
    pub subcommand: ConfigSubcommand,
}

#[derive(Subcommand, Debug)]
pub enum ConfigSubcommand {
    #[command(about = "Show current configuration")]
    Show {
        #[arg(long, help = "Show only the config file path")]
        path: bool,
    },
    #[command(about = "Get a configuration value")]
    Get {
        #[arg(help = "Configuration key (e.g., apptainer_bin, build.default_cpus)")]
        key: String,
    },
    #[command(about = "Set a configuration value")]
    Set {
        #[arg(help = "Configuration key")]
        key: String,
        #[arg(help = "Configuration value")]
        value: String,
    },
    #[command(about = "Create a config file with defaults")]
    Init,
    #[command(about = "Edit config file in default editor")]
    Edit,
    #[command(about = "Validate configuration (exit code 0=valid, 1=invalid)")]
    Valid {
        #[arg(long, short, help = "Quiet mode - only exit codes, no output")]
        quiet: bool,
    },
}

pub fn handle_config(args: &ConfigArgs, config: &AppConfig) -> Result<()> {
    match &args.subcommand {
        ConfigSubcommand::Show { path } => handle_show(*path, config),
        ConfigSubcommand::Get { key } => handle_get(key, config),
        ConfigSubcommand::Set { key, value } => handle_set(key, value, config),
        ConfigSubcommand::Init => handle_init(config),
        ConfigSubcommand::Edit => handle_edit(config),
        ConfigSubcommand::Valid { quiet } => handle_valid(config, *quiet),
    }
}

fn handle_show(show_path: bool, config: &AppConfig) -> Result<()> {
    let config_path = AppConfig::get_user_config_path()?;

    if show_path {
        println!("{}", config_path.display());
        return Ok(());
    }

    // Show config file location
    println!("Config file: {}", utils::style_path(&config_path.display().to_string()));
    if config_path.exists() {
        println!("  {}\n", style_success("exists"));
    } else {
        println!(
            "  {} (use 'condatainer config init' to create)\n",
            style_warning("not found")
        );
    }

    // Show all settings
    println!("{}", utils::style_title("Current Configuration:"));
    println!();

    // Directories
    println!("{}", utils::style_title("Directories:"));
    println!("  program_base:   {}", config.program_base.display());
    println!("  base_dir:       {}", config.base_dir.display());
    println!("  logs_dir:       {}", config.logs_dir.display());
    println!();

    // Binary paths
    println!("{}", utils::style_title("Binaries:"));
    println!("  apptainer_bin:  {}", config.apptainer_bin.display());
    match &config.scheduler_bin {
        Some(bin) => println!("  scheduler_bin:  {}", bin.display()),
        None => println!("  scheduler_bin:  {}", style_warning("not configured")),
    }
    println!();

    // Runtime settings
    println!("{}", utils::style_title("Runtime:"));
    println!("  submit_job:     {}", config.submit_job);
    println!();

    // Build settings
    println!("{}", utils::style_title("Build Configuration:"));
    println!("  default_cpus:    {}", config.build.default_cpus);
    println!("  default_mem_mb:  {}", config.build.default_memory_mb);
    println!(
        "  default_time:    {}s",
        config.build.default_time.as_secs()
    );
    println!("  tmp_size_mb:     {}", config.build.tmp_size_mb);
    println!("  compress_args:   {}", config.build.compress_args);
    println!("  overlay_type:    {}", config.build.overlay_type);
    println!();

    // Show environment variable overrides
    println!("{}", utils::style_title("Environment Variable Overrides:"));
    let env_vars = [
        "CONDATAINER_BASE_DIR",
        "CONDATAINER_LOGS_DIR",
        "CONDATAINER_APPTAINER_BIN",
        "CONDATAINER_SCHEDULER_BIN",
        "CONDATAINER_SUBMIT_JOB",
        "CONDATAINER_BUILD_DEFAULT_CPUS",
        "CONDATAINER_BUILD_DEFAULT_MEM_MB",
        "CONDATAINER_BUILD_DEFAULT_TIME",
    ];
    let mut has_overrides = false;
    for var in &env_vars {
        if let Ok(val) = env::var(var) {
            println!("  {}={}", var, val);
            has_overrides = true;
        }
    }
    if !has_overrides {
        println!("  {}", style_info("none"));
    }

    Ok(())
}

fn handle_get(_key: &str, _config: &AppConfig) -> Result<()> {
    // For now, we'll just return an error since we don't have a full viper-like system
    // In a full implementation, you'd want to load the config file and extract the value
    anyhow::bail!("Config get is not yet fully implemented. Use 'config show' to see all values.");
}

fn handle_set(key: &str, value: &str, _config: &AppConfig) -> Result<()> {
    let known_keys = [
        "base_dir",
        "logs_dir",
        "apptainer_bin",
        "scheduler_bin",
        "scheduler_type",
        "submit_job",
        "build.default_cpus",
        "build.default_mem_mb",
        "build.default_time",
        "build.tmp_size_mb",
        "build.compress_args",
        "build.overlay_type",
    ];

    if !known_keys.contains(&key) {
        utils::print_warning(&format!("'{}' is not a standard config key", key));
    }

    // Validate value based on key type
    if key == "build.default_time" {
        parse_duration(value).context("Invalid duration format")?;
        utils::print_success("Valid duration format");
    }

    // For now, we'll just show a message
    // In a full implementation, you'd want to load, modify, and save the config file
    utils::print_note("This command is not yet fully implemented");
    println!("To set config values, edit the config file directly:");
    println!("  {}", style_hint("condatainer config edit"));

    Ok(())
}

fn handle_init(_config: &AppConfig) -> Result<()> {
    let config_path = AppConfig::get_user_config_path()?;

    // Check if config already exists
    if config_path.exists() {
        print!(
            "{} Config file already exists: {}\nOverwrite? [y/N]: ",
            style_warning("Warning:"),
            config_path.display()
        );
        io::stdout().flush()?;

        let mut response = String::new();
        io::stdin().read_line(&mut response)?;
        let response = response.trim().to_lowercase();

        if response != "y" && response != "yes" {
            utils::print_note("Cancelled");
            return Ok(());
        }
    }

    // Warn if inside container
    if config::is_inside_container() {
        utils::print_warning("Running inside a container - config changes may not persist on the host");
    }

    // Create a new config with defaults
    let mut new_config = AppConfig::new(false, true)?;

    // Force re-detect binaries from current environment
    let updated = force_detect_and_save(&mut new_config)?;

    // Check if apptainer was found
    if !validate_binary(&new_config.apptainer_bin.to_string_lossy()) {
        utils::print_warning("Cannot find apptainer. Please load apptainer module:");
        utils::print_hint("ml apptainer");
    }

    if updated {
        utils::print_success(&format!(
            "Config file created with auto-detected settings: {}",
            config_path.display()
        ));
    } else {
        utils::print_success(&format!(
            "Config file created: {}",
            config_path.display()
        ));
    }

    // Show what was detected
    println!();
    println!("{}", utils::style_title("Detected settings:"));
    println!("  Program base: {}", new_config.program_base.display());
    println!("  Data base:    {}", new_config.base_dir.display());
    if new_config.program_base != new_config.base_dir {
        println!("    {}", style_info("(using separate writable directory for data)"));
    }
    println!("  Apptainer:    {}", new_config.apptainer_bin.display());
    match &new_config.scheduler_bin {
        Some(bin) => {
            let scheduler_type = detect_scheduler_bin()
                .map(|(_, t)| t)
                .unwrap_or_else(|| "unknown".to_string());
            println!("  Scheduler:    {} ({})", bin.display(), scheduler_type);
        }
        None => println!("  Scheduler:    {}", style_warning("not found")),
    }
    println!("  Compression:  {}", new_config.build.compress_args);

    Ok(())
}

fn handle_edit(_config: &AppConfig) -> Result<()> {
    let config_path = AppConfig::get_user_config_path()?;

    // Create config if it doesn't exist
    if !config_path.exists() {
        utils::print_note(&format!(
            "Config file doesn't exist, creating it first..."
        ));
        let config = AppConfig::new(false, true)?;
        config.save_config()?;
    }

    // Get editor from environment
    let editor = env::var("EDITOR").unwrap_or_else(|_| "vi".to_string());

    // Open editor
    let status = Command::new(&editor)
        .arg(&config_path)
        .status()
        .context(format!("Failed to open editor: {}", editor))?;

    if !status.success() {
        anyhow::bail!("Editor exited with non-zero status");
    }

    Ok(())
}

fn handle_valid(config: &AppConfig, quiet: bool) -> Result<()> {
    use std::process;
    
    let mut valid = true;
    
    if !quiet {
        println!("{}", utils::style_title("Validating configuration..."));
        println!();
    }

    // Check apptainer binary
    if validate_binary(&config.apptainer_bin.to_string_lossy()) {
        if !quiet {
            println!(
                "{} Apptainer binary: {}",
                style_success("✓"),
                config.apptainer_bin.display()
            );
        }
    } else {
        if !quiet {
            println!(
                "{} Apptainer binary not found: {}",
                style_error("✗"),
                config.apptainer_bin.display()
            );
        } else {
            utils::print_error(&format!(
                "Apptainer binary not found or not executable: {}",
                config.apptainer_bin.display()
            ));
        }
        valid = false;
    }

    // Check scheduler binary
    if !quiet {
        match &config.scheduler_bin {
            Some(bin) => {
                if validate_binary(&bin.to_string_lossy()) {
                    println!("{} Scheduler binary: {}", style_success("✓"), bin.display());
                } else {
                    println!("{} Scheduler binary not found: {}", style_error("✗"), bin.display());
                    valid = false;
                }
            }
            None => {
                println!(
                    "{} Scheduler binary: {}",
                    style_warning("⚠"),
                    style_info("not configured")
                );
            }
        }
    }

    // Check base_dir writability
    if config::is_writable(&config.base_dir) {
        if !quiet {
            println!("{} Base directory writable: {}", style_success("✓"), config.base_dir.display());
        }
    } else {
        if !quiet {
            println!(
                "{} Base directory not writable: {}",
                style_error("✗"),
                config.base_dir.display()
            );
        } else {
            utils::print_error(&format!(
                "Base directory is not writable: {}",
                config.base_dir.display()
            ));
        }
        valid = false;
    }

    // Check build config
    if !quiet {
        if config.build.default_cpus > 0 {
            println!("{} Build CPUs: {}", style_success("✓"), config.build.default_cpus);
        } else {
            println!(
                "{} Build CPUs must be > 0: {}",
                style_error("✗"),
                config.build.default_cpus
            );
            valid = false;
        }

        if config.build.default_memory_mb > 0 {
            println!(
                "{} Build Memory: {} MB",
                style_success("✓"),
                config.build.default_memory_mb
            );
        } else {
            println!(
                "{} Build Memory must be > 0: {}",
                style_error("✗"),
                config.build.default_memory_mb
            );
            valid = false;
        }

        println!();
    }

    if valid {
        if !quiet {
            utils::print_success("Configuration is valid");
        }
        Ok(())
    } else {
        if !quiet {
            utils::print_error("Configuration has errors");
        }
        process::exit(1);
    }
}
