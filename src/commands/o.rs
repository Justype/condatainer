//! Shortcut command for overlay create

use anyhow::Result;
use clap::Parser;
use std::path::PathBuf;

use crate::commands::overlay::{handle_overlay_create, CreateOptions};
use crate::config::AppConfig;

#[derive(Parser, Debug)]
pub struct OArgs {
    /// Path to the overlay image (default: env.img)
    #[arg(value_name = "IMAGE_PATH", default_value = "env.img")]
    pub path: PathBuf,

    /// Set overlay size (e.g., 500M, 10G)
    #[arg(short = 's', long, default_value = "10G")]
    pub size: String,

    /// Overlay profile: small/balanced/large files
    #[arg(short = 't', long, default_value = "balanced")]
    pub r#type: String,

    /// Filesystem type: ext3 or ext4
    #[arg(long, default_value = "ext3")]
    pub fs: String,

    /// Create a fakeroot-compatible overlay (owned by root)
    #[arg(long)]
    pub fakeroot: bool,

    /// Create a sparse overlay image
    #[arg(long)]
    pub sparse: bool,

    /// Initialize with Conda environment file (.yaml/.yml)
    #[arg(short = 'f', long)]
    pub file: Option<String>,
}

pub fn handle_o(args: &OArgs, config: &AppConfig) -> Result<()> {
    let options = CreateOptions {
        path: args.path.clone(),
        size: args.size.clone(),
        r#type: args.r#type.clone(),
        fs: args.fs.clone(),
        fakeroot: args.fakeroot,
        sparse: args.sparse,
        file: args.file.clone(),
    };
    
    handle_overlay_create(&options, config)
}
