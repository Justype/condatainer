//! Create command - build new SquashFS overlays

use anyhow::{anyhow, Context, Result};
use clap::Parser;
use std::path::PathBuf;

use crate::build::BuildObject;
use crate::config::AppConfig;
use crate::utils;

#[derive(Debug, Parser)]
pub struct CreateArgs {
    /// Packages to install (e.g., samtools/1.16 bcftools)
    #[arg(value_name = "PACKAGES")]
    pub packages: Vec<String>,

    /// Custom name for the resulting overlay file
    #[arg(short = 'n', long)]
    pub name: Option<String>,

    /// Custom prefix path for the overlay file
    #[arg(short = 'p', long)]
    pub prefix: Option<String>,

    /// Path to definition file (.yaml, .sh, .def)
    #[arg(short = 'f', long)]
    pub file: Option<PathBuf>,

    /// Base image to use instead of default
    #[arg(short = 'b', long)]
    pub base_image: Option<PathBuf>,

    /// Remote source URI (e.g., docker://ubuntu:22.04)
    #[arg(short = 's', long)]
    pub source: Option<String>,

    /// Size of temporary overlay
    #[arg(long, default_value = "20G")]
    pub temp_size: String,

    // Compression flags
    /// Use zstd compression level 3
    #[arg(long)]
    pub zstd_fast: bool,

    /// Use zstd compression level 8
    #[arg(long)]
    pub zstd_medium: bool,

    /// Use zstd compression level 14
    #[arg(long)]
    pub zstd: bool,

    /// Use zstd compression level 19
    #[arg(long)]
    pub zstd_high: bool,

    /// Use gzip compression
    #[arg(long)]
    pub gzip: bool,

    /// Use LZ4 compression (default)
    #[arg(long)]
    pub lz4: bool,

    /// Submit build as scheduler job (if scheduler detected)
    #[arg(long)]
    pub submit: bool,
}

/// Main handler for the create command
pub fn handle_create(args: CreateArgs, config: &AppConfig) -> Result<()> {
    utils::print_debug(&format!("Create args: {:?}", args));

    // 1. Validation Logic
    validate_arguments(&args)?;

    // 2. Ensure base image exists
    ensure_base_image(config)?;

    // 3. Handle Compression Config
    let compress_args = determine_compression(&args, config)?;
    utils::print_debug(&format!("Compression args: {}", compress_args));

    // 4. Handle temp size
    let temp_size_mb = parse_temp_size(&args.temp_size)?;
    utils::print_debug(&format!("Temp size: {} MB", temp_size_mb));

    // 5. Normalize package names
    let normalized_packages: Vec<String> = args
        .packages
        .iter()
        .map(|pkg| normalize_name_version(pkg))
        .collect();

    // 6. Execute create based on mode
    if args.source.is_some() {
        // Mode: --source (external image like docker://ubuntu)
        run_create_from_source(&args, config, &compress_args)?;
    } else if args.prefix.is_some() {
        // Mode: --prefix (with --file for YAML/def/sh)
        run_create_with_prefix(&args, config, &compress_args)?;
    } else if args.name.is_some() {
        // Mode: --name (multiple packages or YAML into one sqf)
        run_create_with_name(&args, config, &normalized_packages, &compress_args)?;
    } else {
        // Mode: Default (each package gets its own sqf via BuildObject)
        run_create_packages(&args, config, &normalized_packages, &compress_args)?;
    }

    Ok(())
}

/// Validate command arguments
fn validate_arguments(args: &CreateArgs) -> Result<()> {
    if args.packages.is_empty() && args.file.is_none() && args.source.is_none() {
        return Err(anyhow!(
            "At least one of [packages], --file, or --source must be provided"
        ));
    }

    if !args.packages.is_empty() && args.prefix.is_some() {
        return Err(anyhow!("Packages cannot be used with --prefix"));
    }

    if args.prefix.is_some() && args.name.is_some() {
        return Err(anyhow!(
            "Cannot use both --prefix and --name at the same time"
        ));
    }

    if args.source.is_some() && args.name.is_none() && args.prefix.is_none() {
        return Err(anyhow!(
            "When using --source, either --name or --prefix must be provided"
        ));
    }

    if args.file.is_some() && args.prefix.is_none() {
        return Err(anyhow!("When using --file, --prefix must be provided"));
    }

    Ok(())
}

/// Ensure base image exists
fn ensure_base_image(config: &AppConfig) -> Result<()> {
    crate::apptainer::ensure_base_image(config)
}

/// Determine compression arguments based on flags
fn determine_compression(args: &CreateArgs, config: &AppConfig) -> Result<String> {
    // Priority order for compression flags
    if args.lz4 {
        return Ok("-comp lz4".to_string());
    } else if args.zstd {
        return Ok("-comp zstd -Xcompression-level 14".to_string());
    } else if args.zstd_fast {
        return Ok("-comp zstd -Xcompression-level 3".to_string());
    } else if args.zstd_medium {
        return Ok("-comp zstd -Xcompression-level 8".to_string());
    } else if args.zstd_high {
        return Ok("-comp zstd -Xcompression-level 19".to_string());
    } else if args.gzip {
        return Ok("-comp gzip".to_string());
    }

    // Auto-detect: use zstd-medium if supported
    let apptainer = crate::apptainer::Apptainer::from_config(config);
    if let Ok(version) = apptainer.get_version() {
        if check_zstd_support(&version) {
            utils::print_debug(&format!(
                "Auto-selected zstd-medium compression (Apptainer {} supports zstd)",
                utils::style_number(&version)
            ));
            return Ok("-comp zstd -Xcompression-level 8".to_string());
        }
    }

    // Default to lz4
    Ok("-comp lz4".to_string())
}

/// Check if Apptainer version supports zstd compression
fn check_zstd_support(version: &str) -> bool {
    // zstd support added in Apptainer 1.1.0
    if let Some((major, minor, _)) = parse_version(version) {
        return major > 1 || (major == 1 && minor >= 1);
    }
    false
}

/// Parse version string like "1.2.3" into (major, minor, patch)
fn parse_version(version: &str) -> Option<(u32, u32, u32)> {
    let parts: Vec<&str> = version.split('.').collect();
    if parts.len() < 2 {
        return None;
    }
    let major = parts[0].parse::<u32>().ok()?;
    let minor = parts[1].parse::<u32>().ok()?;
    let patch = parts.get(2).and_then(|p| p.parse::<u32>().ok()).unwrap_or(0);
    Some((major, minor, patch))
}

/// Parse temp size string (e.g., "20G", "1024M") to MB
fn parse_temp_size(size_str: &str) -> Result<u64> {
    let size_str = size_str.trim().to_uppercase();
    
    // Extract number and unit
    let (num_str, unit) = if let Some(idx) = size_str.find(|c: char| c.is_alphabetic()) {
        size_str.split_at(idx)
    } else {
        // No unit, assume MB
        return size_str.parse::<u64>()
            .context("Invalid size format");
    };

    let num: f64 = num_str.parse()
        .context("Invalid number in size")?;

    let multiplier = match unit {
        "K" | "KB" => 1.0 / 1024.0,
        "M" | "MB" => 1.0,
        "G" | "GB" => 1024.0,
        "T" | "TB" => 1024.0 * 1024.0,
        _ => return Err(anyhow!("Unknown size unit: {}", unit)),
    };

    Ok((num * multiplier) as u64)
}

/// Normalize package name/version (replace @ with /)
fn normalize_name_version(name: &str) -> String {
    name.replace('@', "/")
}

/// Run create for individual packages (default mode)
fn run_create_packages(
    _args: &CreateArgs,
    config: &AppConfig,
    packages: &[String],
    _compress_args: &str,
) -> Result<()> {
    use crate::build::{BuildGraph, CondaBuildObject};
    
    utils::print_message("Creating separate overlays for each package...");
    
    // Create build objects for each package
    let mut build_objects: Vec<Box<dyn crate::build::BuildObject>> = Vec::new();
    
    for pkg in packages {
        let bo = CondaBuildObject::new(
            pkg.clone(),
            pkg.clone(), // Single package mode
            &config.images_dir,
            &config.tmp_dir,
        );
        
        utils::print_debug(&format!("[CREATE] BuildObject created for: {}", pkg));
        build_objects.push(Box::new(bo));
    }
    
    // Create build graph
    let mut graph = BuildGraph::new(
        build_objects,
        config.images_dir.clone(),
        config.tmp_dir.clone(),
        false, // Not for scheduler
    )?;
    
    // Run the build
    graph.run()?;
    
    Ok(())
}

/// Run create with a custom name (single overlay for multiple packages or YAML)
fn run_create_with_name(
    args: &CreateArgs,
    config: &AppConfig,
    packages: &[String],
    _compress_args: &str,
) -> Result<()> {
    use crate::build::CondaBuildObject;
    
    let name = args.name.as_ref().unwrap();
    let normalized_name = normalize_name_version(name);
    
    if normalized_name.matches('/').count() > 1 {
        return Err(anyhow!("--name cannot contain more than one '/'"));
    }

    // Build absolute path for the sqf file
    let target_path = config.images_dir.join(format!("{}.sqf", name));

    // Check if already exists
    if target_path.exists() {
        utils::print_message(&format!(
            "Overlay {} already exists at {}. Skipping creation.",
            utils::style_name(&target_path.file_name().unwrap().to_string_lossy()),
            utils::style_path(&target_path.to_string_lossy())
        ));
        return Ok(());
    }

    utils::print_debug(&format!("[CREATE] Creating overlay with name: {}", name));

    // Determine build source
    let build_source = if let Some(file) = &args.file {
        // YAML file mode
        if !file.to_string_lossy().ends_with(".yml") && !file.to_string_lossy().ends_with(".yaml") {
            return Err(anyhow!("File must be .yml or .yaml for conda environments with --name"));
        }
        file.to_string_lossy().to_string()
    } else if !packages.is_empty() {
        // Multiple packages mode - join with commas
        packages.join(",")
    } else {
        return Err(anyhow!("--name requires either packages or --file"));
    };

    // Create CondaBuildObject
    let bo = CondaBuildObject::new(
        normalized_name,
        build_source,
        &config.images_dir,
        &config.tmp_dir,
    );

    // Build it
    bo.build(false)?;
    
    utils::print_success(&format!("Created overlay: {}", utils::style_name(name)));
    
    Ok(())
}

/// Run create with a custom prefix (from external file)
fn run_create_with_prefix(
    args: &CreateArgs,
    config: &AppConfig,
    _compress_args: &str,
) -> Result<()> {
    use crate::build::{CondaBuildObject, DefBuildObject, ShellBuildObject};
    
    let prefix = args.prefix.as_ref().unwrap();
    let file = args.file.as_ref().unwrap();

    if !file.exists() {
        return Err(anyhow!("File {} not found", utils::style_path(&file.to_string_lossy())));
    }

    utils::print_message(&format!(
        "Creating overlay with prefix {} from {}",
        utils::style_name(prefix),
        utils::style_path(&file.to_string_lossy())
    ));

    let file_str = file.to_string_lossy().to_string();
    let normalized_prefix = normalize_name_version(prefix);
    
    // Determine build type based on file extension
    let result = if file_str.ends_with(".yml") || file_str.ends_with(".yaml") {
        // Conda YAML file
        let bo = CondaBuildObject::new(
            normalized_prefix,
            file_str,
            &config.images_dir,
            &config.tmp_dir,
        );
        bo.build(false)
    } else if file_str.ends_with(".def") {
        // Apptainer definition file
        let bo = DefBuildObject::new(
            normalized_prefix,
            file_str,
            &config.images_dir,
            &config.tmp_dir,
        );
        bo.build(false)
    } else if file_str.ends_with(".sh") {
        // Shell script
        let bo = ShellBuildObject::new(
            normalized_prefix,
            file_str,
            &config.images_dir,
            &config.tmp_dir,
            false, // Not a reference overlay
        );
        bo.build(false)
    } else {
        return Err(anyhow!(
            "Unsupported file type. Supported: .yml, .yaml, .sh, .def"
        ));
    };
    
    result?;
    utils::print_success(&format!("Created overlay: {}", utils::style_name(prefix)));
    
    Ok(())
}

/// Run create from external source (docker, def file, etc.)
fn run_create_from_source(
    args: &CreateArgs,
    config: &AppConfig,
    _compress_args: &str,
) -> Result<()> {
    use crate::build::DefBuildObject;
    
    let source = args.source.as_ref().unwrap();
    
    let (target_name, target_path) = if let Some(prefix) = &args.prefix {
        let normalized = normalize_name_version(prefix);
        let path = config.images_dir.join(format!("{}.sqf", normalized));
        (normalized, path)
    } else if let Some(name) = &args.name {
        let normalized = normalize_name_version(name);
        let path = config.images_dir.join(format!("{}.sqf", normalized));
        (normalized, path)
    } else {
        return Err(anyhow!("--source requires either --name or --prefix"));
    };
    
    // Check if already exists
    if target_path.exists() {
        utils::print_message(&format!(
            "Overlay {} already exists at {}. Skipping creation.",
            utils::style_name(&target_name),
            utils::style_path(&target_path.to_string_lossy())
        ));
        return Ok(());
    }

    utils::print_message(&format!(
        "Creating overlay from source: {}",
        utils::style_path(source)
    ));

    // For docker:// sources, we use DefBuildObject which can handle remote sources
    // The source string is passed as the build_source
    let bo = DefBuildObject::new(
        target_name.clone(),
        source.clone(),
        &config.images_dir,
        &config.tmp_dir,
    );

    bo.build(false)?;
    
    utils::print_success(&format!("Created overlay: {}", utils::style_name(&target_name)));
    
    Ok(())
}
