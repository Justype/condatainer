//! Build object trait and common implementations

use anyhow::{Context, Result};
use std::path::{Path, PathBuf};

/// Type of build object
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BuildType {
    Conda,
    Def,
    Shell,
    Reference,
}

/// Script specifications for scheduler jobs
#[derive(Debug, Clone, Default)]
pub struct ScriptSpecs {
    /// Raw scheduler flags (e.g., SLURM flags)
    pub raw_flags: Vec<String>,
    /// Number of CPUs required
    pub ncpus: Option<usize>,
    /// Memory required in MB
    pub memory_mb: Option<u64>,
    /// Time limit
    pub time_limit: Option<String>,
}

impl ScriptSpecs {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_cpus(mut self, ncpus: usize) -> Self {
        self.ncpus = Some(ncpus);
        self
    }

    pub fn with_memory(mut self, memory_mb: u64) -> Self {
        self.memory_mb = Some(memory_mb);
        self
    }

    pub fn with_time(mut self, time_limit: String) -> Self {
        self.time_limit = Some(time_limit);
        self
    }

    pub fn add_flag<S: Into<String>>(mut self, flag: S) -> Self {
        self.raw_flags.push(flag.into());
        self
    }
}

/// Trait representing a build target with its metadata and dependencies
pub trait BuildObject: std::fmt::Debug {
    /// Get name and version (e.g., "samtools/1.16")
    fn name_version(&self) -> &str;
    
    /// Get build source (path to script, YAML file, etc.)
    fn build_source(&self) -> &str;
    
    /// Get list of dependencies
    fn dependencies(&self) -> &[String];
    
    /// Check if already installed
    fn is_installed(&self) -> bool;
    
    /// Get build type
    fn build_type(&self) -> BuildType;
    
    /// Path to temporary overlay during build
    fn tmp_overlay_path(&self) -> &Path;
    
    /// Path to final target overlay
    fn target_overlay_path(&self) -> &Path;
    
    /// Path to container directory
    fn cnt_dir_path(&self) -> &Path;
    
    /// Get scheduler specifications
    fn script_specs(&self) -> Option<&ScriptSpecs>;
    
    /// Check if requires scheduler submission
    fn requires_scheduler(&self) -> bool {
        self.script_specs()
            .map(|specs| !specs.raw_flags.is_empty())
            .unwrap_or(false)
    }
    
    /// Build the overlay
    fn build(&self, build_deps: bool) -> Result<()>;
    
    /// Get list of missing dependencies
    fn get_missing_dependencies(&self) -> Result<Vec<String>>;
    
    /// Create temporary overlay
    fn create_tmp_overlay(&self, force: bool) -> Result<()>;
    
    /// Cleanup after build (successful or failed)
    fn cleanup(&self, failed: bool) -> Result<()>;
}

/// Base implementation providing common functionality
#[derive(Debug, Clone)]
pub struct BaseBuildObject {
    pub name_version: String,
    pub build_source: String,
    pub dependencies: Vec<String>,
    pub tmp_overlay_path: PathBuf,
    pub target_overlay_path: PathBuf,
    pub cnt_dir_path: PathBuf,
    pub ncpus: usize,
    pub is_remote: bool,
    pub script_specs: Option<ScriptSpecs>,
}

impl BaseBuildObject {
    pub fn new(
        name_version: String,
        build_source: String,
        images_dir: &Path,
        tmp_dir: &Path,
    ) -> Self {
        let tmp_overlay_path = tmp_dir.join(format!("{}.tmp.ext3", name_version.replace('/', "--")));
        let target_overlay_path = images_dir.join(format!("{}.sqf", name_version.replace('/', "--")));
        let cnt_dir_path = images_dir.join(name_version.replace('/', "--"));

        Self {
            name_version,
            build_source,
            dependencies: Vec::new(),
            tmp_overlay_path,
            target_overlay_path,
            cnt_dir_path,
            ncpus: super::default_build_ncpus(),
            is_remote: false,
            script_specs: None,
        }
    }

    pub fn is_installed(&self) -> bool {
        self.target_overlay_path.exists()
    }

    pub fn get_missing_dependencies(&self) -> Result<Vec<String>> {
        // TODO: Implement actual dependency checking
        // For now, assume all dependencies are satisfied
        Ok(Vec::new())
    }

    pub fn create_tmp_overlay(&self, force: bool) -> Result<()> {
        if self.tmp_overlay_path.exists() {
            if !force {
                return Err(super::BuildError::TmpOverlayExists(
                    self.tmp_overlay_path.display().to_string()
                ).into());
            }
            std::fs::remove_file(&self.tmp_overlay_path)?;
        }

        // Use overlay module to create the temporary overlay
        crate::overlay::create_tmp_overlay(&self.tmp_overlay_path, 20480)?;
        
        Ok(())
    }


    pub fn cleanup(&self, _failed: bool) -> Result<()> {
        // Remove remote build source if downloaded
        if self.is_remote && !self.build_source.is_empty() {
            if let Err(e) = std::fs::remove_file(&self.build_source) {
                crate::utils::print_debug(&format!(
                    "Failed to remove remote build script {}: {}",
                    crate::utils::style_path(&self.build_source),
                    e
                ));
            } else {
                crate::utils::print_debug(&format!(
                    "Removed remote build script {}",
                    crate::utils::style_path(&self.build_source)
                ));
            }
        }

        // Remove tmp overlay
        if !self.tmp_overlay_path.as_os_str().is_empty() {
            if let Err(e) = std::fs::remove_file(&self.tmp_overlay_path) {
                crate::utils::print_debug(&format!(
                    "Failed to remove temporary overlay {}: {}",
                    crate::utils::style_path(&self.tmp_overlay_path.to_string_lossy()),
                    e
                ));
            } else {
                crate::utils::print_debug(&format!(
                    "Removed temporary overlay {}",
                    crate::utils::style_path(&self.tmp_overlay_path.to_string_lossy())
                ));
            }
        }

        Ok(())
    }
}
