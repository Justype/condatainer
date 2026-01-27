//! Completion helpers for shell autocomplete

use anyhow::Result;
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

use crate::config::AppConfig;
use crate::utils;

/// Get installed overlays as name -> path mapping
pub fn get_installed_overlays() -> Result<HashMap<String, PathBuf>> {
    let config = AppConfig::new(false, false)?;
    let mut overlays = HashMap::new();

    if !config.images_dir.exists() {
        return Ok(overlays);
    }

    for entry in fs::read_dir(&config.images_dir)? {
        let entry = entry?;
        let path = entry.path();
        
        if let Some(ext) = path.extension() {
            if ext == "sqf" || ext == "sqsh" || ext == "squashfs" {
                if let Some(stem) = path.file_stem() {
                    let name = stem.to_string_lossy().replace("--", "/");
                    overlays.insert(name, path);
                }
            }
        }
    }

    Ok(overlays)
}

/// Get local overlay files in current directory
pub fn get_local_overlays(to_complete: &str) -> Vec<String> {
    let path_dir = if let Some(idx) = to_complete.rfind('/') {
        &to_complete[..=idx]
    } else {
        ""
    };

    let dir_for_read = if path_dir.is_empty() { "." } else { path_dir };

    let Ok(entries) = fs::read_dir(dir_for_read) else {
        return Vec::new();
    };

    let mut suggestions = Vec::new();
    for entry in entries.flatten() {
        let path = entry.path();
        if path.is_dir() {
            continue;
        }

        if let Some(name) = path.file_name() {
            let name_str = name.to_string_lossy();
            if !utils::is_overlay(&name_str) {
                continue;
            }
            
            // Exclude .sif files - only show overlay files (.img, .sqf, .sqsh)
            if name_str.ends_with(".sif") {
                continue;
            }

            let candidate = if path_dir.is_empty() {
                name_str.to_string()
            } else {
                format!("{}{}", path_dir, name_str)
            };

            if to_complete.is_empty() || candidate.starts_with(to_complete) {
                suggestions.push(candidate);
            }
        }
    }

    suggestions
}

/// Get overlay suggestions for completion (includes data overlays)
pub fn get_overlay_suggestions(include_data: bool, to_complete: &str) -> Vec<String> {
    let mut suggestions = Vec::new();

    // Add installed overlays
    if let Ok(installed) = get_installed_overlays() {
        let config = AppConfig::new(false, false).ok();
        
        for (name, path) in installed {
            let is_app = if let Some(ref cfg) = config {
                is_app_overlay(&path, &cfg.images_dir)
            } else {
                true
            };

            if include_data || is_app {
                if to_complete.is_empty() || name.starts_with(to_complete) {
                    suggestions.push(name);
                }
            }
        }
    }

    // Add local overlay files
    suggestions.extend(get_local_overlays(to_complete));

    suggestions.sort();
    suggestions.dedup();
    suggestions
}

/// Get base image suggestions for completion (app overlays only, slash count < 2)
pub fn get_base_image_suggestions(to_complete: &str) -> Vec<String> {
    let mut suggestions = Vec::new();

    // Add app overlays (slash count < 2)
    if let Ok(installed) = get_installed_overlays() {
        let config = AppConfig::new(false, false).ok();
        
        for (name, path) in installed {
            let slash_count = name.matches('/').count();
            if slash_count < 2 {
                let is_app = if let Some(ref cfg) = config {
                    is_app_overlay(&path, &cfg.images_dir)
                } else {
                    true
                };

                if is_app && (to_complete.is_empty() || name.starts_with(to_complete)) {
                    suggestions.push(name);
                }
            }
        }
    }

    // Add local .sif and .sqf files
    if let Ok(entries) = fs::read_dir(".") {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_file() {
                if let Some(ext) = path.extension() {
                    if ext == "sif" || ext == "sqf" {
                        if let Some(name) = path.file_name() {
                            let candidate = name.to_string_lossy().to_string();
                            if to_complete.is_empty() || candidate.starts_with(to_complete) {
                                suggestions.push(candidate);
                            }
                        }
                    }
                }
            }
        }
    }

    suggestions.sort();
    suggestions.dedup();
    suggestions
}

/// Check if path is an app overlay (not a data overlay)
fn is_app_overlay(path: &Path, images_dir: &Path) -> bool {
    path_within_dir(path, images_dir)
}

/// Generate enhanced bash completion script with dynamic suggestions
pub fn generate_enhanced_bash_completion() {
    use clap::CommandFactory;
    use crate::Cli;
    
    // First generate the base clap completion
    let mut cmd = Cli::command();
    let mut buf = Vec::new();
    clap_complete::generate(clap_complete::Shell::Bash, &mut cmd, "condatainer", &mut buf);
    let base_completion = String::from_utf8_lossy(&buf);
    
    // Print the enhanced wrapper
    print!(r#"#!/usr/bin/env bash
# Enhanced bash completion for condatainer with dynamic overlay suggestions

# Load the base clap-generated completion
{}

# Enhanced completion function that adds dynamic overlay suggestions
_condatainer_enhanced() {{
    local cur prev words cword
    
    # Simple variable setup (works without bash-completion package)
    COMPREPLY=()
    cur="${{COMP_WORDS[COMP_CWORD]}}"
    prev="${{COMP_WORDS[COMP_CWORD-1]}}"
    words=("${{COMP_WORDS[@]}}")
    cword=$COMP_CWORD

    local cmd=""
    local i
    
    # Find the subcommand
    for ((i=1; i < cword; i++)); do
        if [[ "${{words[i]}}" != -* ]]; then
            cmd="${{words[i]}}"
            break
        fi
    done

    # Handle different commands
    case "$cmd" in
        e)
            # Check if completing after -b flag
            case "$prev" in
                -b|--base-image)
                    local images
                    images=$(condatainer complete-base-image "$cur" 2>/dev/null)
                    if [[ -n "$images" ]]; then
                        COMPREPLY=($(compgen -W "$images" -- "$cur"))
                        return 0
                    fi
                    ;;
                *)
                    # For 'e' command, suggest overlays for positional args
                    if [[ "$cur" != -* ]]; then
                        local overlays
                        overlays=$(condatainer complete-overlay true "$cur" 2>/dev/null)
                        if [[ -n "$overlays" ]]; then
                            COMPREPLY=($(compgen -W "$overlays" -- "$cur"))
                            return 0
                        fi
                    fi
                    ;;
            esac
            ;;
        exec)
            # For 'exec' command, check if completing after -o or -b
            case "$prev" in
                -o|--overlay)
                    local overlays
                    overlays=$(condatainer complete-overlay true "$cur" 2>/dev/null)
                    if [[ -n "$overlays" ]]; then
                        COMPREPLY=($(compgen -W "$overlays" -- "$cur"))
                        return 0
                    fi
                    ;;
                -b|--base-image)
                    local images
                    images=$(condatainer complete-base-image "$cur" 2>/dev/null)
                    if [[ -n "$images" ]]; then
                        COMPREPLY=($(compgen -W "$images" -- "$cur"))
                        return 0
                    fi
                    ;;
                *)
                    # Positional args - could be overlays or commands
                    if [[ "$cur" != -* ]]; then
                        local overlays
                        overlays=$(condatainer complete-overlay true "$cur" 2>/dev/null)
                        if [[ -n "$overlays" ]]; then
                            COMPREPLY=($(compgen -W "$overlays" -- "$cur"))
                            return 0
                        fi
                    fi
                    ;;
            esac
            ;;
    esac

    # Fall back to the default clap completion
    _condatainer "$@"
}}

# Override the completion function
complete -F _condatainer_enhanced condatainer
"#, base_completion);
}

/// Generate enhanced zsh completion script with dynamic suggestions
pub fn generate_enhanced_zsh_completion() {
    use clap::CommandFactory;
    use crate::Cli;
    
    // First generate the base clap completion
    let mut cmd = Cli::command();
    let mut buf = Vec::new();
    clap_complete::generate(clap_complete::Shell::Zsh, &mut cmd, "condatainer", &mut buf);
    let base_completion = String::from_utf8_lossy(&buf);
    
    // Print the enhanced wrapper
    print!(r#"#compdef condatainer

# Enhanced zsh completion for condatainer with dynamic overlay suggestions

{}

# Custom completion function for overlay arguments
_condatainer_overlays() {{
    local overlays
    overlays=($(condatainer complete-overlay true "" 2>/dev/null))
    _describe 'overlays' overlays
}}

# Custom completion function for base image arguments
_condatainer_base_images() {{
    local images
    images=($(condatainer complete-base-image "" 2>/dev/null))
    _describe 'base images' images
}}

# Override completions for specific commands
_condatainer_e_overlays() {{
    case $words[$CURRENT-1] in
        -b|--base-image)
            _condatainer_base_images
            return 0
            ;;
        *)
            if [[ $words[$CURRENT] != -* ]]; then
                _condatainer_overlays
                return 0
            fi
            ;;
    esac
    return 1
}}

_condatainer_exec_overlays() {{
    case $words[$CURRENT-1] in
        -o|--overlay)
            _condatainer_overlays
            return 0
            ;;
        -b|--base-image)
            _condatainer_base_images
            return 0
            ;;
        *)
            if [[ $words[$CURRENT] != -* ]]; then
                _condatainer_overlays
                return 0
            fi
            ;;
    esac
    return 1
}}

# Hook into the main completion function
(( $+functions[_condatainer_commands] )) && {{
    _condatainer_e_original=${{functions[_condatainer__e_cmd]}}
    _condatainer__e_cmd() {{
        _condatainer_e_overlays || ${{_condatainer_e_original}}
    }}
    
    _condatainer_exec_original=${{functions[_condatainer__exec_cmd]}}
    _condatainer__exec_cmd() {{
        _condatainer_exec_overlays || ${{_condatainer_exec_original}}
    }}
}}
"#, base_completion);
}

/// Generate enhanced fish completion script with dynamic suggestions
pub fn generate_enhanced_fish_completion() {
    use clap::CommandFactory;
    use crate::Cli;
    
    // First generate the base clap completion
    let mut cmd = Cli::command();
    let mut buf = Vec::new();
    clap_complete::generate(clap_complete::Shell::Fish, &mut cmd, "condatainer", &mut buf);
    let base_completion = String::from_utf8_lossy(&buf);
    
    // Print the enhanced wrapper
    print!(r#"# Enhanced fish completion for condatainer with dynamic overlay suggestions

{}

# Dynamic overlay completion for 'e' command
complete -c condatainer -n '__fish_seen_subcommand_from e; and not __fish_seen_argument -s b -l base-image' -a "(condatainer complete-overlay true '' 2>/dev/null)"

# Dynamic base image completion for 'e' command with -b flag
complete -c condatainer -n '__fish_seen_subcommand_from e; and __fish_seen_argument -s b -l base-image' -a "(condatainer complete-base-image '' 2>/dev/null)"

# Dynamic overlay completion for 'exec' command with -o flag
complete -c condatainer -n '__fish_seen_subcommand_from exec; and __fish_seen_argument -s o -l overlay' -a "(condatainer complete-overlay true '' 2>/dev/null)"

# Dynamic base image completion for 'exec' command with -b flag
complete -c condatainer -n '__fish_seen_subcommand_from exec; and __fish_seen_argument -s b -l base-image' -a "(condatainer complete-base-image '' 2>/dev/null)"

# Dynamic overlay completion for 'exec' positional args (when no flags present)
complete -c condatainer -n '__fish_seen_subcommand_from exec; and not __fish_seen_argument -s o -l overlay -s b -l base-image' -a "(condatainer complete-overlay true '' 2>/dev/null)"
"#, base_completion);
}

/// Check if path is within a directory
fn path_within_dir(path: &Path, dir: &Path) -> bool {
    if let Ok(rel) = path.strip_prefix(dir) {
        !rel.starts_with("..")
    } else {
        false
    }
}
