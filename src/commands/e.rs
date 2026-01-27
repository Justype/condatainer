use anyhow::{anyhow, Context, Result};
use clap::Parser;
use std::collections::HashMap;
use std::env;
use std::fs;
use std::io::IsTerminal;
use std::path::{Path, PathBuf};
use std::process::Command;

use crate::apptainer::{Apptainer, ExecOptions};
use crate::config::AppConfig;
use crate::utils::{is_img, is_overlay, is_sqf, print_debug, print_error, print_message, print_note, print_warning, style_highlight, style_path};

#[derive(Parser, Debug)]
pub struct EArgs {
    #[arg(short = 'r', long = "read-only", help = "Mount .img overlays as read-only (default: writable)")]
    pub read_only: bool,
    
    #[arg(short = 'b', long = "base-image", help = "Base image to use instead of default")]
    pub base_image: Option<String>,
    
    #[arg(short = 'n', long = "no-autoload", help = "Disable autoloading 'env.img' from current directory")]
    pub no_autoload: bool,
    
    #[arg(short = 'f', long = "fakeroot", help = "Run container with fakeroot privileges")]
    pub fakeroot: bool,
    
    #[arg(trailing_var_arg = true, allow_hyphen_values = true, help = "Overlay files to mount (defaults to env.img)")]
    pub overlays: Vec<String>,
}

pub fn handle_e(args: &EArgs, config: &AppConfig) -> Result<()> {
    let apptainer = Apptainer::from_config(config);
    
    // Determine base image
    let base_image = if let Some(ref custom_base) = args.base_image {
        PathBuf::from(custom_base)
    } else {
        config.base_image.clone()
    };
    
    if !base_image.exists() {
        return Err(anyhow!("Base image not found at {}", base_image.display()));
    }
    
    // Get installed overlays map
    let installed_overlays = get_installed_overlays_map(config)?;
    
    // Parse overlays and commands (look for -- separator)
    let (overlay_list, command_list) = parse_overlays_and_commands(&args.overlays);
    
    // Resolve overlay paths
    let mut overlay_paths = resolve_overlay_paths(&overlay_list, &installed_overlays, config)?;
    
    // Autoload env.img if no .img overlay and not disabled
    if !args.no_autoload && !overlay_paths.iter().any(|p| is_img(&p.to_string_lossy())) {
        let pwd = env::current_dir()?;
        let local_env_path = pwd.join("env.img");
        if local_env_path.exists() {
            print_message(&format!(
                "Autoload env.img at {}",
                style_path(&local_env_path.to_string_lossy())
            ));
            overlay_paths.push(local_env_path);
        }
    }
    
    // Validate only one .img overlay
    validate_single_img(&overlay_paths)?;
    
    // Reorder overlays: put .img at the end
    let ordered_overlays = reorder_overlays(overlay_paths);
    
    // Check if we have a writable .img overlay
    let writable_img = !args.read_only;
    let has_writable_img = writable_img && ordered_overlays.iter().any(|p| is_img(&p.to_string_lossy()));
    
    // Determine fakeroot requirement
    let mut use_fakeroot = args.fakeroot;
    if has_writable_img && !use_fakeroot {
        if let Some(img_path) = ordered_overlays.iter().find(|p| is_img(&p.to_string_lossy())) {
            let uid_status = inspect_image_uid_status(img_path)?;
            if uid_status == 0 {
                print_note(&format!(
                    "Root overlay {}. --fakeroot added automatically",
                    style_path(&img_path.file_name().unwrap().to_string_lossy())
                ));
                use_fakeroot = true;
            } else if uid_status == 2 {
                print_warning(&format!(
                    "{}'s inner UID is different. --fakeroot added automatically",
                    style_path(&img_path.file_name().unwrap().to_string_lossy())
                ));
                use_fakeroot = true;
            }
        }
    }
    
    // Build environment variables
    let mut env_vars = build_env_vars(&ordered_overlays, has_writable_img, config)?;
    
    // Track container layer
    let layer = if let Ok(in_cnt) = env::var("IN_CONDATAINER") {
        in_cnt.parse::<i32>().unwrap_or(1) + 1
    } else {
        1
    };
    
    if layer > 1 {
        print_warning("You are trying to mount an .img overlay inside an existing CondaTainer environment. This may lead to unexpected behavior.");
    }
    
    let ps1_prefix = if layer == 1 {
        "CNT".to_string()
    } else {
        format!("CNT_{}", layer)
    };
    
    env_vars.push(format!("IN_CONDATAINER={}", layer));
    env_vars.push(format!("PS1={} \\[\\e[0;34m\\]\\w\\[\\e[0m\\]> ", ps1_prefix));
    
    // Display environment notes if in interactive shell
    if is_interactive_shell() && command_list.is_empty() {
        display_env_notes(&ordered_overlays, has_writable_img)?;
    }
    
    // Build exec options
    let mut opts = ExecOptions::new();
    
    // Add overlays
    for overlay in &ordered_overlays {
        let overlay_str = overlay.to_string_lossy().to_string();
        if writable_img && is_img(&overlay_str) {
            print_debug(&format!("[EXEC] Making overlay {} writable.", style_path(&overlay_str)));
            opts = opts.overlay(overlay_str);
        } else {
            opts = opts.overlay(format!("{}:ro", overlay_str));
        }
    }
    
    // Add environment variables
    opts = opts.envs(env_vars);
    
    // Add fakeroot
    if use_fakeroot {
        opts = opts.fakeroot();
    }
    
    // Use provided command or default to bash
    let final_command = if command_list.is_empty() {
        vec!["bash".to_string()]
    } else {
        command_list
    };
    
    // Execute
    let result = apptainer.exec(&base_image, &final_command, Some(opts));
    
    match result {
        Ok(_) => Ok(()),
        Err(e) => {
            if has_writable_img && is_interactive_shell() {
                print_error("Command execution failed inside the container with writable image.");
                print_note("Please check the troubleshooting documentation for exec issues.");
            }
            Err(e)
        }
    }
}

/// Parse overlays and commands, looking for -- separator
fn parse_overlays_and_commands(args: &[String]) -> (Vec<String>, Vec<String>) {
    if let Some(pos) = args.iter().position(|s| s == "--") {
        let overlays = args[..pos].to_vec();
        let commands = args[pos + 1..].to_vec();
        (overlays, commands)
    } else {
        // All are overlays, no command
        (args.to_vec(), Vec::new())
    }
}

/// Get installed overlays as name/version -> path mapping
fn get_installed_overlays_map(config: &AppConfig) -> Result<HashMap<String, PathBuf>> {
    let mut overlays = HashMap::new();
    
    if !config.images_dir.exists() {
        return Ok(overlays);
    }
    
    let entries = fs::read_dir(&config.images_dir).with_context(|| {
        format!("Unable to list installed overlays in {}", config.images_dir.display())
    })?;
    
    for entry in entries {
        let entry = entry?;
        if entry.file_type()?.is_dir() {
            continue;
        }
        
        let file_name = entry.file_name();
        let file_name_str = file_name.to_string_lossy();
        
        if !is_overlay(&file_name_str) {
            continue;
        }
        
        // Convert filename to name/version format
        let name_version = file_name_str
            .trim_end_matches(".sqf")
            .trim_end_matches(".sqsh")
            .trim_end_matches(".squashfs")
            .trim_end_matches(".img")
            .trim_end_matches(".sif");
        
        // Convert samtools--1.21 to samtools/1.21
        let normalized = name_version.replace("--", "/");
        overlays.insert(normalized, entry.path());
    }
    
    Ok(overlays)
}

/// Normalize package name/version format
fn normalize_name_version(s: &str) -> String {
    s.trim()
        .replace('=', "/")
        .replace('@', "/")
        .replace("--", "/")
}

/// Resolve overlay names/paths to absolute paths
fn resolve_overlay_paths(
    overlays: &[String],
    installed_overlays: &HashMap<String, PathBuf>,
    config: &AppConfig,
) -> Result<Vec<PathBuf>> {
    let mut paths = Vec::new();
    
    for overlay in overlays {
        // Check if it's a file path
        if is_overlay(overlay) {
            let path = PathBuf::from(overlay);
            if path.exists() {
                paths.push(path.canonicalize()?);
            } else {
                return Err(anyhow!("Overlay file {} not found.", overlay));
            }
        } else {
            // Treat as package name
            let normalized = normalize_name_version(overlay);
            if let Some(path) = installed_overlays.get(&normalized) {
                paths.push(path.clone());
            } else {
                // Try as filename in images_dir
                let formatted = normalized.replace('/', "--");
                let overlay_path = config.images_dir.join(format!("{}.sqf", formatted));
                if overlay_path.exists() {
                    paths.push(overlay_path);
                } else {
                    return Err(anyhow!(
                        "Overlay {} not found. Try 'condatainer list' to see installed overlays.",
                        style_highlight(overlay)
                    ));
                }
            }
        }
    }
    
    Ok(paths)
}

/// Validate only one .img overlay is used
fn validate_single_img(overlays: &[PathBuf]) -> Result<()> {
    let img_count = overlays
        .iter()
        .filter(|p| is_img(&p.to_string_lossy()))
        .count();
    
    if img_count > 1 {
        return Err(anyhow!("Only one .img overlay can be used at a time."));
    }
    
    Ok(())
}

/// Reorder overlays to put .img at the end
fn reorder_overlays(overlays: Vec<PathBuf>) -> Vec<PathBuf> {
    let mut img_overlays = Vec::new();
    let mut other_overlays = Vec::new();
    
    for overlay in overlays {
        if is_img(&overlay.to_string_lossy()) {
            img_overlays.push(overlay);
        } else {
            other_overlays.push(overlay);
        }
    }
    
    other_overlays.extend(img_overlays);
    other_overlays
}

/// Build environment variables for the container
fn build_env_vars(
    overlays: &[PathBuf],
    has_writable_img: bool,
    config: &AppConfig,
) -> Result<Vec<String>> {
    let mut env_vars = Vec::new();
    
    // Build PATH environment variable
    let mut paths = vec!["/usr/sbin".to_string(), "/usr/bin".to_string()];
    
    for overlay in overlays {
        let file_name = overlay
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("");
        let normalized = normalize_name_version(file_name);
        let slash_count = normalized.matches('/').count();
        
        // Skip reference overlays (more than 1 slash)
        if slash_count > 1 {
            continue;
        }
        
        let overlay_str = overlay.to_string_lossy();
        if is_img(&overlay_str) {
            paths.insert(0, "/ext3/env/bin".to_string());
        } else if is_sqf(&overlay_str) {
            paths.insert(0, format!("/cnt/{}/bin", normalized));
        }
    }
    
    env_vars.push(format!("PATH={}", paths.join(":")));
    env_vars.push("LC_ALL=C.UTF-8".to_string());
    env_vars.push("LANG=C.UTF-8".to_string());
    
    // Add .env file configurations
    let env_configs = load_env_configs(overlays)?;
    for (key, value) in env_configs {
        env_vars.push(format!("{}={}", key, value));
    }
    
    // Check if there's any img overlay
    let has_img = overlays.iter().any(|p| is_img(&p.to_string_lossy()));
    
    // Add conda/mamba environment variables if there's an img overlay
    if has_img {
        env_vars.push("MAMBA_ROOT_PREFIX=/ext3/env".to_string());
        env_vars.push("CONDA_PREFIX=/ext3/env".to_string());
        env_vars.push("CONDA_DEFINE_ENV=env".to_string());
        env_vars.push("RETICULATE_PYTHON=/ext3/env/bin/python".to_string());
        
        // CNT_CONDA_PREFIX is only set when the img is writable
        if has_writable_img {
            env_vars.push("CNT_CONDA_PREFIX=/ext3/env".to_string());
        }
    }
    
    // Add CONDATINER_DIR bind
    env_vars.push(format!("CONDATAINER_DIR={}", config.base_dir.display()));
    
    Ok(env_vars)
}

/// Load environment configurations from .env files
fn load_env_configs(overlays: &[PathBuf]) -> Result<HashMap<String, String>> {
    let mut configs = HashMap::new();
    
    for overlay_path in overlays {
        let env_path = overlay_path.with_extension(
            format!("{}.env", overlay_path.extension().and_then(|s| s.to_str()).unwrap_or(""))
        );
        
        if !env_path.exists() {
            continue;
        }
        
        let content = fs::read_to_string(&env_path)?;
        for line in content.lines() {
            let line = line.trim();
            if line.starts_with('#') || line.is_empty() {
                continue;
            }
            
            if let Some((key, value)) = line.split_once('=') {
                if configs.contains_key(key) {
                    print_message(&format!(
                        "Environment variable {} is defined in multiple overlays. Using the value from {}.",
                        key,
                        overlay_path.file_name().unwrap().to_string_lossy()
                    ));
                }
                configs.insert(key.to_string(), value.to_string());
            }
        }
    }
    
    Ok(configs)
}

/// Display environment variable notes
fn display_env_notes(overlays: &[PathBuf], has_writable_img: bool) -> Result<()> {
    let mut env_notes = HashMap::new();
    
    // Load .env files and extract notes
    for overlay_path in overlays {
        let env_path = overlay_path.with_extension(
            format!("{}.env", overlay_path.extension().and_then(|s| s.to_str()).unwrap_or(""))
        );
        
        if !env_path.exists() {
            continue;
        }
        
        let content = fs::read_to_string(&env_path)?;
        for line in content.lines() {
            if let Some(note_content) = line.strip_prefix("#ENVNOTE:") {
                if let Some((key, note)) = note_content.split_once('=') {
                    env_notes.insert(key.trim().to_string(), note.trim().to_string());
                }
            } else if !line.starts_with('#') && line.contains('=') {
                if let Some((key, value)) = line.split_once('=') {
                    if !env_notes.contains_key(key.trim()) {
                        env_notes.insert(key.trim().to_string(), value.trim().to_string());
                    }
                }
            }
        }
    }
    
    if has_writable_img {
        env_notes.insert("CNT_CONDA_PREFIX".to_string(), "/ext3/env".to_string());
    }
    
    if !env_notes.is_empty() {
        let max_len = env_notes.keys().map(|k| k.len()).max().unwrap_or(0);
        print_message("Overlay envs:");
        for (key, note) in env_notes {
            println!("  {}: {}", style_highlight(&format!("{:width$}", key, width = max_len)), note);
        }
        println!();
    }
    
    Ok(())
}

/// Check if running in an interactive shell
fn is_interactive_shell() -> bool {
    std::io::stdin().is_terminal()
}

/// Inspect image UID status using debugfs
/// Returns: 0=root owned, 1=same UID, 2=different UID, 9=failed to check
fn inspect_image_uid_status(img_path: &Path) -> Result<i32> {
    if !is_img(&img_path.to_string_lossy()) || !img_path.exists() {
        return Ok(9);
    }
    
    let output = Command::new("debugfs")
        .args(&["-R", "stat upper", &img_path.to_string_lossy()])
        .output();
    
    match output {
        Ok(output) => {
            let stdout = String::from_utf8_lossy(&output.stdout);
            if let Some(caps) = regex::Regex::new(r"User:\s+(\d+)")?.captures(&stdout) {
                if let Some(uid_str) = caps.get(1) {
                    let uid: u32 = uid_str.as_str().parse()?;
                    let current_uid = unsafe { libc::getuid() };
                    
                    if uid == 0 {
                        return Ok(0);
                    } else if uid == current_uid {
                        return Ok(1);
                    } else {
                        return Ok(2);
                    }
                }
            }
            Ok(9)
        }
        Err(_) => {
            print_warning(&format!(
                "Failed to check fakeroot status of image {}. Assuming fakeroot.",
                img_path.display()
            ));
            Ok(9)
        }
    }
}
