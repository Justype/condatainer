use anyhow::{Context, Result};
use clap::Parser;
use colored::*;
use regex::Regex;
use std::collections::HashMap;
use std::fs;
use std::io::{self, Write};
use std::path::PathBuf;

use crate::config::AppConfig;
use crate::utils::{print_error, print_success, print_warning, style_highlight};

#[derive(Parser, Debug)]
pub struct RemoveArgs {
    #[arg(
        help = "Names or search terms (AND logic applied)",
        required = true
    )]
    pub terms: Vec<String>,
}

pub fn handle_remove(args: &RemoveArgs, config: &AppConfig) -> Result<()> {
    // Normalize terms (replace =, @, -- with /)
    let normalized_terms: Vec<String> = args
        .terms
        .iter()
        .map(|term| normalize_name_version(term))
        .collect();

    // Get installed overlays
    let installed_overlays = get_installed_overlays_with_paths(config)?;

    if installed_overlays.is_empty() {
        print_warning("No installed overlays found.");
        return Ok(());
    }

    // Filter overlays based on search terms
    let filtered_overlays = filter_overlays(&installed_overlays, &normalized_terms);

    if filtered_overlays.is_empty() {
        print_warning("No matching installed overlays found.");
        return Ok(());
    }

    // Display overlays to be removed with highlighted search terms
    println!("Overlays to be removed:");
    for overlay_name in &filtered_overlays {
        let highlighted = highlight_terms(overlay_name, &normalized_terms);
        println!(" - {}", highlighted);
    }

    // Ask for confirmation
    print!("{}", "[CNT] Are you sure? Cannot be undone. [y/N]: ");
    io::stdout().flush()?;
    
    let mut input = String::new();
    io::stdin().read_line(&mut input)?;
    let choice = input.trim().to_lowercase();
    
    if choice != "y" {
        println!("Removal cancelled.");
        return Ok(());
    }

    // Remove overlays
    for overlay_name in &filtered_overlays {
        if let Some(overlay_path) = installed_overlays.get(overlay_name) {
            match fs::remove_file(overlay_path) {
                Ok(_) => {
                    print_success(&format!(
                        "Overlay {} removed.",
                        style_highlight(overlay_name)
                    ));
                    
                    // Also remove .env file if it exists
                    let env_path = format!("{}.env", overlay_path.display());
                    if PathBuf::from(&env_path).exists() {
                        let _ = fs::remove_file(&env_path);
                    }
                }
                Err(e) => {
                    print_error(&format!(
                        "Failed to remove overlay {}: {}",
                        style_highlight(overlay_name),
                        e
                    ));
                }
            }
        }
    }

    Ok(())
}

/// Normalize package spec formats so that `name/version`, `name=version`, `name@version`
/// are treated the same. This converts them all to use forward slashes.
fn normalize_name_version(name_version: &str) -> String {
    name_version
        .trim()
        .replace('=', "/")
        .replace('@', "/")
        .replace("--", "/")
}

/// Get installed overlays as a map of name/version -> file path
fn get_installed_overlays_with_paths(config: &AppConfig) -> Result<HashMap<String, PathBuf>> {
    let mut overlays = HashMap::new();

    if !config.images_dir.exists() {
        return Ok(overlays);
    }

    let entries = fs::read_dir(&config.images_dir).with_context(|| {
        format!(
            "Unable to list installed overlays in {}",
            config.images_dir.display()
        )
    })?;

    for entry in entries {
        let entry = entry?;
        let file_name = entry.file_name();
        let file_name_str = file_name.to_string_lossy();

        if entry.file_type()?.is_dir() {
            continue;
        }

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

/// Check if a filename is an overlay file
fn is_overlay(filename: &str) -> bool {
    filename.ends_with(".sqf")
        || filename.ends_with(".sqsh")
        || filename.ends_with(".squashfs")
        || filename.ends_with(".img")
        || filename.ends_with(".sif")
}

/// Filter overlays based on search terms using AND logic
fn filter_overlays(
    installed_overlays: &HashMap<String, PathBuf>,
    terms: &[String],
) -> Vec<String> {
    let mut filtered = Vec::new();

    // Check if first term is an exact match - if so, use exact match mode
    if let Some(first_term) = terms.first() {
        if installed_overlays.contains_key(first_term) {
            // Exact match mode - look for each term exactly
            for term in terms {
                if installed_overlays.contains_key(term) {
                    filtered.push(term.clone());
                } else {
                    print_warning(&format!(
                        "Overlay {} not found among installed overlays.",
                        style_highlight(term)
                    ));
                }
            }
            return filtered;
        }
    }

    // Search mode - filter using AND logic with case-insensitive regex matching
    for (overlay_name, _) in installed_overlays {
        let matches_all_terms = terms.iter().all(|term| {
            // Use case-insensitive matching
            let pattern = regex::escape(term);
            if let Ok(re) = Regex::new(&format!("(?i){}", pattern)) {
                re.is_match(overlay_name)
            } else {
                // If regex fails, fall back to simple contains
                overlay_name.to_lowercase().contains(&term.to_lowercase())
            }
        });

        if matches_all_terms {
            filtered.push(overlay_name.clone());
        }
    }

    // Sort for consistent output
    filtered.sort();
    filtered
}

/// Highlight search terms in the overlay name
fn highlight_terms(overlay_name: &str, terms: &[String]) -> String {
    let mut result = overlay_name.to_string();

    for term in terms {
        let pattern = regex::escape(term);
        if let Ok(re) = Regex::new(&format!("(?i){}", pattern)) {
            result = re
                .replace_all(&result, |caps: &regex::Captures| {
                    caps[0].yellow().to_string()
                })
                .to_string();
        }
    }

    result
}
