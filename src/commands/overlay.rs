//! Overlay management commands - create, resize, check, info, chown

use anyhow::{anyhow, Context, Result};
use clap::{Args, Subcommand};
use std::path::PathBuf;

use crate::config::AppConfig;
use crate::utils;

/// Shared options for overlay create command
#[derive(Debug, Clone)]
pub struct CreateOptions {
    pub path: PathBuf,
    pub size: String,
    pub r#type: String,
    pub fs: String,
    pub fakeroot: bool,
    pub sparse: bool,
    pub file: Option<String>,
}

#[derive(Debug, Args)]
pub struct OverlayArgs {
    #[command(subcommand)]
    pub command: OverlayCommand,
}

#[derive(Debug, Subcommand)]
pub enum OverlayCommand {
    /// Create a new sparse overlay image
    #[command(about = "Create a new sparse overlay image")]
    Create {
        /// Path to the overlay image (default: env.img)
        #[arg(value_name = "IMAGE_PATH", default_value = "env.img")]
        path: PathBuf,

        /// Set overlay size (e.g., 500M, 10G)
        #[arg(short = 's', long, default_value = "10G")]
        size: String,

        /// Overlay profile: small/balanced/large files
        #[arg(short = 't', long, default_value = "balanced")]
        r#type: String,

        /// Filesystem type: ext3 or ext4
        #[arg(long, default_value = "ext3")]
        fs: String,

        /// Create a fakeroot-compatible overlay (owned by root)
        #[arg(long)]
        fakeroot: bool,

        /// Create a sparse overlay image
        #[arg(long)]
        sparse: bool,

        /// Conda environment file to initialize (.yaml or .yml)
        #[arg(short = 'f', long)]
        file: Option<String>,
    },

    /// Expand or shrink an existing overlay image
    #[command(about = "Expand or shrink an existing overlay image")]
    Resize {
        /// Path to the overlay image
        #[arg(value_name = "IMAGE_PATH")]
        path: PathBuf,

        /// New size (e.g., 20G, 2048M)
        #[arg(short = 's', long, required = true)]
        size: String,
    },

    /// Display disk usage and filesystem stats
    #[command(about = "Display disk usage and filesystem stats")]
    Info {
        /// Path to the overlay image
        #[arg(value_name = "IMAGE_PATH")]
        path: PathBuf,
    },

    /// Verify filesystem integrity (e2fsck)
    #[command(about = "Verify filesystem integrity (e2fsck)")]
    Check {
        /// Path to the overlay image
        #[arg(value_name = "IMAGE_PATH")]
        path: PathBuf,

        /// Force check even if filesystem appears clean
        #[arg(short = 'f', long)]
        force: bool,
    },

    /// Recursively set internal files to specific UID/GID
    #[command(about = "Recursively set internal files to specific UID/GID")]
    Chown {
        /// Path to the overlay image
        #[arg(value_name = "IMAGE_PATH")]
        path: PathBuf,

        /// User ID to set
        #[arg(short = 'u', long)]
        uid: Option<u32>,

        /// Group ID to set
        #[arg(short = 'g', long)]
        gid: Option<u32>,

        /// Set UID and GID to 0 (root); overrides -u and -g
        #[arg(long)]
        root: bool,

        /// Path inside the overlay to change
        #[arg(short = 'p', long, default_value = "/ext3")]
        internal_path: String,
    },
}

pub fn handle_overlay(args: &OverlayArgs, _config: &AppConfig) -> Result<()> {
    match &args.command {
        OverlayCommand::Create {
            path,
            size,
            r#type,
            fs,
            fakeroot,
            sparse,
            file,
        } => handle_create(path, size, r#type, fs, *fakeroot, *sparse, file.as_ref()),

        OverlayCommand::Resize { path, size } => handle_resize(path, size),

        OverlayCommand::Info { path } => handle_info(path),

        OverlayCommand::Check { path, force } => handle_check(path, *force),

        OverlayCommand::Chown {
            path,
            uid,
            gid,
            root,
            internal_path,
        } => handle_chown(path, *uid, *gid, *root, internal_path),
    }
}

/// Handle overlay create command (public API for shortcut commands)
pub fn handle_overlay_create(options: &CreateOptions, _config: &AppConfig) -> Result<()> {
    handle_create(
        &options.path,
        &options.size,
        &options.r#type,
        &options.fs,
        options.fakeroot,
        options.sparse,
        options.file.as_ref(),
    )
}

/// Handle overlay create command
fn handle_create(
    path: &PathBuf,
    size_str: &str,
    type_str: &str,
    fs_type: &str,
    fakeroot: bool,
    sparse: bool,
    file: Option<&String>,
) -> Result<()> {
    // Check if file already exists
    if path.exists() {
        utils::print_error(&format!("File already exists: {}", path.display()));
        return Err(anyhow!("File already exists"));
    }

    // Parse size
    let size_mb = parse_size_to_mb(size_str)
        .context(format!("Invalid size format '{}'", size_str))?;

    // Parse profile type
    let profile = match type_str.to_lowercase().as_str() {
        "small" | "conda" | "python" => crate::overlay::Profile::Small,
        "large" | "data" | "genome" => crate::overlay::Profile::Large,
        "balanced" | "default" => crate::overlay::Profile::Default,
        _ => {
            utils::print_warning(&format!("Unknown profile type '{}', using 'balanced'", type_str));
            crate::overlay::Profile::Default
        }
    };

    // Create overlay
    if fakeroot {
        crate::overlay::create_for_root(path, size_mb, profile, sparse)?;
    } else {
        crate::overlay::create_for_current_user(path, size_mb, profile, sparse)?;
    }

    // Initialize conda environment if this is a .img overlay
    if utils::is_img(&path.to_string_lossy()) {
        initialize_conda_env(path, file, fakeroot)?;
    }

    Ok(())
}

/// Initialize conda environment inside the overlay from an environment file or with zlib
fn initialize_conda_env(path: &PathBuf, file: Option<&String>, fakeroot: bool) -> Result<()> {
    let config = AppConfig::new(false, true)?;
    let apptainer = crate::apptainer::Apptainer::from_config(&config);
    
    let overlay_path_str = path.to_string_lossy().to_string();
    
    if let Some(env_file) = file {
        // Initialize with environment file
        let env_file_path = std::path::Path::new(env_file);
        if !env_file_path.exists() {
            utils::print_error(&format!("Conda environment file not found: {}", env_file));
            return Err(anyhow!("Environment file not found"));
        }
        
        let abs_file_path = env_file_path.canonicalize()?;
        utils::print_message(&format!(
            "Initializing conda environment using {}...",
            utils::style_path(&abs_file_path.to_string_lossy())
        ));
        
        // Run mm-create with the file
        let mut exec_opts = crate::apptainer::ExecOptions::new()
            .overlay(&overlay_path_str)
            .env("MAMBA_ROOT_PREFIX=/ext3/env")
            .env("CONDA_PREFIX=/ext3/env")
            .env("CONDA_DEFINE_ENV=env")
            .env("RETICULATE_PYTHON=/ext3/env/bin/python")
            .env("CNT_CONDA_PREFIX=/ext3/env")
            .capture_output();
        
        if fakeroot {
            exec_opts = exec_opts.fakeroot();
        }
        
        let result = apptainer.exec(
            &config.base_image,
            &["mm-create", "-f", &abs_file_path.to_string_lossy(), "-y"],
            Some(exec_opts),
        );
        
        match result {
            Ok(_) => {
                // Clean up
                let mut clean_opts = crate::apptainer::ExecOptions::new()
                    .overlay(&overlay_path_str)
                    .env("MAMBA_ROOT_PREFIX=/ext3/env")
                    .env("CONDA_PREFIX=/ext3/env")
                    .env("CONDA_DEFINE_ENV=env")
                    .env("RETICULATE_PYTHON=/ext3/env/bin/python")
                    .env("CNT_CONDA_PREFIX=/ext3/env")
                    .capture_output();
                
                if fakeroot {
                    clean_opts = clean_opts.fakeroot();
                }
                
                let _ = apptainer.exec(
                    &config.base_image,
                    &["mm-clean", "-a", "-y"],
                    Some(clean_opts),
                );
                
                utils::print_success(&format!(
                    "Conda env is created inside {}",
                    utils::style_path(&path.to_string_lossy())
                ));
                Ok(())
            }
            Err(e) => {
                utils::print_error(&format!(
                    "Failed to create conda env inside {}",
                    utils::style_path(&path.to_string_lossy())
                ));
                Err(anyhow!("Failed to initialize conda environment: {}", e))
            }
        }
    } else {
        // Initialize with minimal zlib package
        utils::print_message("Initializing minimal conda environment with small package (zlib)...");
        
        let mut exec_opts = crate::apptainer::ExecOptions::new()
            .overlay(&overlay_path_str)
            .env("MAMBA_ROOT_PREFIX=/ext3/env")
            .env("CONDA_PREFIX=/ext3/env")
            .env("CONDA_DEFINE_ENV=env")
            .env("RETICULATE_PYTHON=/ext3/env/bin/python")
            .env("CNT_CONDA_PREFIX=/ext3/env")
            .capture_output();
        
        if fakeroot {
            exec_opts = exec_opts.fakeroot();
        }
        
        let result = apptainer.exec(
            &config.base_image,
            &["mm-create", "zlib", "-y"],
            Some(exec_opts),
        );
        
        match result {
            Ok(_) => {
                // Clean up
                let mut clean_opts = crate::apptainer::ExecOptions::new()
                    .overlay(&overlay_path_str)
                    .env("MAMBA_ROOT_PREFIX=/ext3/env")
                    .env("CONDA_PREFIX=/ext3/env")
                    .env("CONDA_DEFINE_ENV=env")
                    .env("RETICULATE_PYTHON=/ext3/env/bin/python")
                    .env("CNT_CONDA_PREFIX=/ext3/env")
                    .capture_output();
                
                if fakeroot {
                    clean_opts = clean_opts.fakeroot();
                }
                
                let _ = apptainer.exec(
                    &config.base_image,
                    &["mm-clean", "-a", "-y"],
                    Some(clean_opts),
                );
                
                utils::print_success(&format!(
                    "Conda env is created inside {}",
                    utils::style_path(&path.to_string_lossy())
                ));
                Ok(())
            }
            Err(e) => {
                utils::print_error(&format!(
                    "Failed to create conda env inside {}",
                    utils::style_path(&path.to_string_lossy())
                ));
                Err(anyhow!("Failed to initialize conda environment: {}", e))
            }
        }
    }
}

/// Handle overlay resize command
fn handle_resize(path: &PathBuf, size_str: &str) -> Result<()> {
    let size_mb = parse_size_to_mb(size_str)
        .context(format!("Invalid size format '{}'", size_str))?;

    crate::overlay::resize(path, size_mb)
}

/// Handle overlay info command
fn handle_info(path: &PathBuf) -> Result<()> {
    let stats = crate::overlay::get_stats(path)?;

    // Calculate usage
    let (used_bytes, disk_pct) = stats.usage();
    let total_bytes = stats.total_blocks * stats.block_size;
    let inode_pct = stats.inode_usage();
    let used_inodes = stats.total_inodes - stats.free_inodes;

    // Display information
    println!("{}", utils::style_title("Overlay Statistics"));
    println!("  {:15} {}", "File:", utils::style_path(&path.to_string_lossy()));
    println!("  {:15} {}", "State:", stats.filesystem_state);
    println!("  {:15} {}", "Last Mounted:", stats.last_mounted);
    println!("{}", utils::style_title("Usage"));
    println!(
        "  {:15} {} / {} ({:.1}%)",
        "Disk Space:",
        format_bytes(used_bytes),
        format_bytes(total_bytes),
        disk_pct
    );
    println!(
        "  {:15} {} / {} ({:.1}%)",
        "Inodes:",
        used_inodes,
        stats.total_inodes,
        inode_pct
    );

    Ok(())
}

/// Handle overlay check command
fn handle_check(path: &PathBuf, force: bool) -> Result<()> {
    crate::overlay::check_integrity(path, force)?;

    utils::print_success(&format!(
        "Filesystem check completed for {}",
        utils::style_path(&path.to_string_lossy())
    ));

    Ok(())
}

/// Handle overlay chown command
fn handle_chown(
    path: &PathBuf,
    uid: Option<u32>,
    gid: Option<u32>,
    root: bool,
    internal_path: &str,
) -> Result<()> {
    let (target_uid, target_gid) = if root {
        (0, 0)
    } else {
        let current_uid = unsafe { libc::getuid() };
        let current_gid = unsafe { libc::getgid() };
        (
            uid.unwrap_or(current_uid),
            gid.unwrap_or(current_gid),
        )
    };

    crate::overlay::change_ownership(path, target_uid, target_gid, internal_path)
}

/// Parse size string (e.g., "20G", "1024M") to megabytes
fn parse_size_to_mb(size_str: &str) -> Result<u32> {
    let size_str = size_str.trim().to_uppercase();
    
    // Extract number and unit
    let (num_str, unit) = if let Some(idx) = size_str.find(|c: char| c.is_alphabetic()) {
        size_str.split_at(idx)
    } else {
        // No unit, assume MB
        return size_str.parse::<u32>()
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

    Ok((num * multiplier) as u32)
}

/// Format bytes into human-readable format
fn format_bytes(bytes: i64) -> String {
    const UNITS: &[&str] = &["B", "KB", "MB", "GB", "TB"];
    let mut size = bytes as f64;
    let mut unit_idx = 0;

    while size >= 1024.0 && unit_idx < UNITS.len() - 1 {
        size /= 1024.0;
        unit_idx += 1;
    }

    format!("{:.2} {}", size, UNITS[unit_idx])
}
