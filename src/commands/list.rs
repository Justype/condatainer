use anyhow::{Context, Result};
use clap::Args;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::Path;

use crate::config::AppConfig;
use crate::utils::{self, is_img, is_sqf, is_sif};

#[derive(Args, Debug)]
pub struct ListArgs {
    #[arg(help = "Search terms to filter overlays (AND logic)")]
    pub terms: Vec<String>,

    #[arg(short = 'd', long = "delete", visible_alias = "remove", visible_short_alias = 'r', help = "Delete listed overlays after confirmation (used with terms)")]
    pub delete: bool,
}

pub fn handle_list(args: &ListArgs, config: &AppConfig) -> Result<()> {
    let filters = normalize_filters(&args.terms);

    let app_overlays = collect_app_overlays(&config.images_dir, &filters)?;
    let ref_overlays = collect_reference_overlays(&config.images_dir, &filters)?;

    let mut printed = false;

    if !app_overlays.is_empty() {
        printed = true;
        println!("Available app overlays:");
        
        let mut names: Vec<&String> = app_overlays.keys().collect();
        names.sort();
        
        let name_width = names.iter().map(|n| n.len()).max().unwrap_or(0);
        
        for name in names {
            let versions = &app_overlays[name];
            let colored_versions: Vec<String> = versions
                .iter()
                .map(|v| {
                    if v == "(system app overlay)" {
                        v.clone()
                    } else {
                        utils::style_info(v)
                    }
                })
                .collect();
            
            let values = colored_versions.join(", ");
            let name_field = format!("{:<width$}", name, width = name_width);
            println!(" {}: {}", utils::style_name(&name_field), values);
        }
    }

    if !ref_overlays.is_empty() {
        if printed {
            println!();
        }
        println!("Available reference overlays:");
        for ref_name in &ref_overlays {
            println!(" {}", utils::style_name(ref_name));
        }
        printed = true;
    }

    if !printed {
        utils::print_warning("No installed overlays match the provided search terms.");
    }

    // Handle deletion if --delete flag is set
    if args.delete && !args.terms.is_empty() {
        // Collect all matching overlay names for deletion
        let mut to_delete = Vec::new();
        
        // Collect from app overlays
        for (name, versions) in &app_overlays {
            for version in versions {
                if version == "(system app overlay)" || version == "(system app/env)" {
                    to_delete.push(name.clone());
                } else {
                    to_delete.push(format!("{}/{}", name, version));
                }
            }
        }
        
        // Collect from ref overlays
        to_delete.extend(ref_overlays.clone());
        
        if to_delete.is_empty() {
            utils::print_warning("No matching overlays to delete.");
            return Ok(());
        }

        // Use remove command to delete
        use crate::commands::remove::{RemoveArgs, handle_remove};
        let remove_args = RemoveArgs {
            terms: to_delete,
        };
        handle_remove(&remove_args, config)?;
    }

    Ok(())
}

fn collect_app_overlays(
    images_dir: &Path,
    filters: &[String],
) -> Result<HashMap<String, Vec<String>>> {
    if !images_dir.exists() {
        return Ok(HashMap::new());
    }

    let entries = fs::read_dir(images_dir)
        .with_context(|| format!("Unable to list app overlays in {}", images_dir.display()))?;

    let mut grouped: HashMap<String, HashSet<String>> = HashMap::new();

    for entry in entries {
        let entry = entry?;
        let file_name = entry.file_name();
        let file_name_str = file_name.to_string_lossy();

        // Skip directories
        if entry.file_type()?.is_dir() {
            continue;
        }

        // Check if it's an overlay file
        if !is_overlay(&file_name_str) {
            continue;
        }

        // App overlays have at most one "--" delimiter
        let delim_count = file_name_str.matches("--").count();
        if delim_count > 1 {
            continue;
        }

        // Remove extension
        let name_version = file_name_str
            .trim_end_matches(".sqf")
            .trim_end_matches(".sqsh")
            .trim_end_matches(".squashfs")
            .trim_end_matches(".img")
            .trim_end_matches(".sif");

        let normalized = normalize_name_version(name_version).to_lowercase();
        if !matches_filters(&normalized, filters) {
            continue;
        }

        let (name, version) = if let Some(pos) = name_version.find("--") {
            let (n, v) = name_version.split_at(pos);
            (n.to_string(), v[2..].to_string())
        } else {
            (name_version.to_string(), "(system app/env)".to_string())
        };

        if name.is_empty() {
            continue;
        }

        grouped
            .entry(name)
            .or_insert_with(HashSet::new)
            .insert(version);
    }

    let mut result = HashMap::new();
    for (name, versions) in grouped {
        let mut version_list: Vec<String> = versions.into_iter().collect();
        version_list.sort();
        result.insert(name, version_list);
    }

    Ok(result)
}

fn collect_reference_overlays(images_dir: &Path, filters: &[String]) -> Result<Vec<String>> {
    if !images_dir.exists() {
        return Ok(Vec::new());
    }

    let entries = fs::read_dir(images_dir).with_context(|| {
        format!(
            "Unable to list reference overlays in {}",
            images_dir.display()
        )
    })?;

    let mut names = Vec::new();

    for entry in entries {
        let entry = entry?;
        let file_name = entry.file_name();
        let file_name_str = file_name.to_string_lossy();

        // Skip directories
        if entry.file_type()?.is_dir() {
            continue;
        }

        // Check if it's an overlay file
        if !is_overlay(&file_name_str) {
            continue;
        }

        // Reference overlays have more than one "--" delimiter
        let delim_count = file_name_str.matches("--").count();
        if delim_count <= 1 {
            continue;
        }

        // Remove extension
        let name_version = file_name_str
            .trim_end_matches(".sqf")
            .trim_end_matches(".sqsh")
            .trim_end_matches(".squashfs")
            .trim_end_matches(".img")
            .trim_end_matches(".sif");

        let normalized = normalize_name_version(name_version).to_lowercase();
        if !matches_filters(&normalized, filters) {
            continue;
        }

        let display_name = name_version.replace("--", "/");
        names.push(display_name);
    }

    names.sort();
    Ok(names)
}

fn normalize_filters(filters: &[String]) -> Vec<String> {
    filters
        .iter()
        .filter_map(|f| {
            let trimmed = f.trim();
            if trimmed.is_empty() {
                None
            } else {
                Some(normalize_name_version(trimmed).to_lowercase())
            }
        })
        .collect()
}

fn matches_filters(name: &str, filters: &[String]) -> bool {
    if filters.is_empty() {
        return true;
    }
    filters.iter().all(|filter| name.contains(filter))
}

fn is_overlay(filename: &str) -> bool {
    is_img(filename) || is_sqf(filename) || is_sif(filename)
}

/// Normalize name/version by removing special characters
fn normalize_name_version(s: &str) -> String {
    s.chars()
        .map(|c| if c.is_alphanumeric() || c == '-' { c } else { '-' })
        .collect()
}
