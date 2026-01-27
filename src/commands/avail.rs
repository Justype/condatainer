use anyhow::{Context, Result};
use clap::Args;
use regex::Regex;
use serde::Deserialize;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::Path;

use crate::config::AppConfig;
use crate::utils::{self, is_img, is_sqf, is_sif};

#[derive(Args, Debug)]
pub struct AvailArgs {
    #[arg(help = "Search terms to filter build scripts (AND logic)")]
    pub terms: Vec<String>,

    #[arg(long, help = "Show only local build scripts")]
    pub local: bool,

    #[arg(long, help = "Show only remote build scripts from GitHub")]
    pub remote: bool,

    #[arg(short = 'i', long = "install", visible_alias = "add", visible_short_alias = 'a', help = "Install the selected build scripts (used with terms)")]
    pub install: bool,
}

#[derive(Debug, Clone)]
struct PackageInfo {
    name: String,
    path: String,
    is_container: bool,
    is_installed: bool,
    is_remote: bool,
}

#[derive(Debug, Deserialize)]
struct RemoteScriptEntry {
    relative_path: String,
    #[allow(dead_code)]
    deps: Vec<String>,
    #[allow(dead_code)]
    #[serde(alias = "sbatch")]  // Support old field name from Go version
    requires_scheduler: bool,
    #[allow(dead_code)]
    whatis: String,
}

pub fn handle_avail(args: &AvailArgs, config: &AppConfig) -> Result<()> {
    let filters = normalize_filters(&args.terms);

    // Determine which sources to include
    let include_local = !args.remote; // include local unless --remote only
    let include_remote = !args.local; // include remote unless --local only

    // Get installed overlays for marking
    let installed_overlays = get_installed_overlays(config)?;

    let mut packages = Vec::new();

    // Get local scripts
    if include_local {
        let local_scripts = get_local_build_scripts(config)?;
        for (name, info) in local_scripts {
            packages.push(PackageInfo {
                name: name.clone(),
                path: info.0,
                is_container: info.1,
                is_installed: installed_overlays.contains(&name),
                is_remote: false,
            });
        }
    }

    // Get remote scripts
    if include_remote {
        match get_remote_build_scripts() {
            Ok(remote_scripts) => {
                for (name, info) in remote_scripts {
                    // Skip if already have local version (local takes precedence)
                    if include_local && packages.iter().any(|p| p.name == name) {
                        continue;
                    }
                    packages.push(PackageInfo {
                        name: name.clone(),
                        path: info.0,
                        is_container: info.1,
                        is_installed: installed_overlays.contains(&name),
                        is_remote: true,
                    });
                }
            }
            Err(e) => {
                utils::print_debug(&format!("Failed to fetch remote scripts: {}", e));
            }
        }
    }

    // Filter packages
    let filtered = filter_packages(&packages, &filters);

    if filtered.is_empty() {
        if args.local {
            utils::print_warning("No matching local build scripts found.");
        } else if args.remote {
            utils::print_warning("No matching remote build scripts found.");
        } else {
            utils::print_warning("No matching build scripts found.");
        }
        return Ok(());
    }

    // Sort by name
    let mut sorted = filtered;
    sorted.sort_by(|a, b| a.name.cmp(&b.name));

    // Print results
    for pkg in &sorted {
        println!("{}", format_package_line(pkg, &filters));
    }

    // Print summary when showing both
    if !args.local && !args.remote {
        let local_count = sorted.iter().filter(|p| !p.is_remote).count();
        let remote_count = sorted.iter().filter(|p| p.is_remote).count();
        if remote_count > 0 {
            println!();
            utils::print_message(&format!(
                "Found {} local and {} remote build scripts.",
                utils::style_number(&local_count.to_string()),
                utils::style_number(&remote_count.to_string())
            ));
        }
    }

    // Handle installation if --install flag is set
    if args.install && !args.terms.is_empty() {
        // Filter uninstalled packages
        let uninstalled: Vec<String> = sorted
            .iter()
            .filter(|p| !p.is_installed)
            .map(|p| p.name.clone())
            .collect();

        if uninstalled.is_empty() {
            utils::print_warning("All matching packages are already installed.");
            return Ok(());
        }

        // Check if in interactive shell
        if utils::is_interactive_shell() {
            println!();
            utils::print_message("Packages to be installed:");
            for pkg in &uninstalled {
                println!(" - {}", utils::style_name(pkg));
            }
            
            use std::io::{self, Write};
            print!("{}" , "[CNT] Do you want to install these packages? [y/N]: ");
            io::stdout().flush()?;
            
            let mut input = String::new();
            io::stdin().read_line(&mut input)?;
            let choice = input.trim().to_lowercase();
            
            if choice != "y" {
                println!("Installation cancelled.");
                return Ok(());
            }

            // Install packages using create command
            use crate::commands::create::{CreateArgs, handle_create};
            let create_args = CreateArgs {
                packages: uninstalled,
                name: None,
                prefix: None,
                file: None,
                base_image: None,
                source: None,
                temp_size: "20G".to_string(),
                zstd_fast: false,
                zstd_medium: false,
                zstd: false,
                zstd_high: false,
                gzip: false,
                lz4: false,
                submit: false,
            };
            handle_create(create_args, config)?;
        } else {
            utils::print_warning("Cannot install in non-interactive mode. Use 'condatainer create' manually.");
        }
    }

    Ok(())
}

fn get_installed_overlays(config: &AppConfig) -> Result<HashSet<String>> {
    let mut installed = HashSet::new();

    if !config.images_dir.exists() {
        return Ok(installed);
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
        installed.insert(normalized);
    }

    Ok(installed)
}

fn get_local_build_scripts(config: &AppConfig) -> Result<HashMap<String, (String, bool)>> {
    let mut scripts = HashMap::new();

    let build_scripts_dir = &config.build_scripts_dir;
    if !build_scripts_dir.exists() {
        return Ok(scripts);
    }

    walk_dir(build_scripts_dir, build_scripts_dir, &mut scripts)?;

    Ok(scripts)
}

fn walk_dir(
    dir: &Path,
    base_dir: &Path,
    scripts: &mut HashMap<String, (String, bool)>,
) -> Result<()> {
    let entries = fs::read_dir(dir)?;

    for entry in entries {
        let entry = entry?;
        let path = entry.path();

        if path.is_dir() {
            walk_dir(&path, base_dir, scripts)?;
            continue;
        }

        let file_name = entry.file_name();
        let file_name_str = file_name.to_string_lossy();

        // Skip .py and .sh files (helper scripts)
        if file_name_str.ends_with(".py") || file_name_str.ends_with(".sh") {
            continue;
        }

        // Skip template files
        if path.to_string_lossy().contains("template") {
            continue;
        }

        // Generate key (relative path)
        let rel_path = path
            .strip_prefix(base_dir)
            .ok()
            .and_then(|p| p.to_str())
            .unwrap_or("");

        // Handle .def files
        let is_container = rel_path.ends_with(".def");
        let mut key = rel_path.to_string();

        if is_container {
            // Skip base image/overlay scripts
            if rel_path.starts_with("base_image") || rel_path.starts_with("base-overlay") {
                continue;
            }
            key = key.trim_end_matches(".def").to_string();
        }

        scripts.insert(key, (path.to_string_lossy().to_string(), is_container));
    }

    Ok(())
}

fn get_remote_build_scripts() -> Result<HashMap<String, (String, bool)>> {
    let url = format!(
        "https://raw.githubusercontent.com/{}/main/metadata/build-scripts.json.gz",
        crate::config::GITHUB_REPO
    );

    // Fetch and decompress
    let response = ureq::get(&url)
        .timeout(std::time::Duration::from_secs(30))
        .call()
        .context("Failed to fetch remote metadata")?;

    let mut decoder = flate2::read::GzDecoder::new(response.into_reader());
    let mut data = String::new();
    std::io::Read::read_to_string(&mut decoder, &mut data)
        .context("Failed to decompress metadata")?;

    // Parse JSON
    let metadata: HashMap<String, RemoteScriptEntry> =
        serde_json::from_str(&data).context("Failed to parse metadata")?;

    // Convert to our format
    let mut scripts = HashMap::new();
    for (name, entry) in metadata {
        let is_container = entry.relative_path.ends_with(".def");
        scripts.insert(name, (entry.relative_path, is_container));
    }

    Ok(scripts)
}

fn normalize_filters(filters: &[String]) -> Vec<String> {
    filters
        .iter()
        .filter_map(|f| {
            let trimmed = f.trim();
            if trimmed.is_empty() {
                None
            } else {
                Some(trimmed.to_lowercase())
            }
        })
        .collect()
}

fn filter_packages(packages: &[PackageInfo], filters: &[String]) -> Vec<PackageInfo> {
    if filters.is_empty() {
        return packages.to_vec();
    }

    packages
        .iter()
        .filter(|pkg| {
            let name_lower = pkg.name.to_lowercase();
            filters.iter().all(|filter| name_lower.contains(filter))
        })
        .cloned()
        .collect()
}

fn format_package_line(pkg: &PackageInfo, filters: &[String]) -> String {
    let mut line = pkg.name.clone();

    // Build suffix
    let mut suffixes = Vec::new();
    if pkg.is_installed {
        suffixes.push(utils::style_success("installed"));
    }
    if pkg.is_container {
        suffixes.push("container".to_string());
    }
    if pkg.is_remote {
        suffixes.push(utils::style_hint("remote"));
    }

    if !suffixes.is_empty() {
        line = format!("{} ({})", line, suffixes.join(", "));
    }

    // Highlight search terms if present
    if !filters.is_empty() {
        for term in filters {
            let re = Regex::new(&format!("(?i){}", regex::escape(term))).unwrap();
            line = re
                .replace_all(&line, |caps: &regex::Captures| {
                    utils::style_highlight(&caps[0])
                })
                .to_string();
        }
    }

    line
}

fn is_overlay(filename: &str) -> bool {
    is_img(filename) || is_sqf(filename) || is_sif(filename)
}
