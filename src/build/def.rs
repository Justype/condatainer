//! Apptainer definition file-based build object implementation

use anyhow::Result;
use std::path::Path;

use super::object::{BaseBuildObject, BuildObject, BuildType, ScriptSpecs};

/// Apptainer definition file-based build object
#[derive(Debug, Clone)]
pub struct DefBuildObject {
    pub base: BaseBuildObject,
}

impl DefBuildObject {
    pub fn new(name_version: String, build_source: String, images_dir: &Path, tmp_dir: &Path) -> Self {
        let base = BaseBuildObject::new(name_version, build_source, images_dir, tmp_dir);
        Self { base }
    }
}

impl BuildObject for DefBuildObject {
    fn name_version(&self) -> &str {
        &self.base.name_version
    }

    fn build_source(&self) -> &str {
        &self.base.build_source
    }

    fn dependencies(&self) -> &[String] {
        &self.base.dependencies
    }

    fn is_installed(&self) -> bool {
        self.base.is_installed()
    }

    fn build_type(&self) -> BuildType {
        BuildType::Def
    }

    fn tmp_overlay_path(&self) -> &Path {
        &self.base.tmp_overlay_path
    }

    fn target_overlay_path(&self) -> &Path {
        &self.base.target_overlay_path
    }

    fn cnt_dir_path(&self) -> &Path {
        &self.base.cnt_dir_path
    }

    fn script_specs(&self) -> Option<&ScriptSpecs> {
        self.base.script_specs.as_ref()
    }

    fn build(&self, _build_deps: bool) -> Result<()> {
        let target_path = self.base.target_overlay_path.canonicalize()
            .unwrap_or_else(|_| self.base.target_overlay_path.clone());
        let styled_overlay = crate::utils::style_name(self.name_version());

        // Check if already exists
        if target_path.exists() {
            crate::utils::print_message(&format!(
                "Overlay {} already exists at {}. Skipping creation.",
                styled_overlay,
                crate::utils::style_path(&target_path.to_string_lossy())
            ));
            return Ok(());
        }

        let build_mode = if self.requires_scheduler() {
            "scheduler"
        } else {
            "local"
        };

        crate::utils::print_message(&format!(
            "Building overlay {} ({} build) from {}",
            styled_overlay,
            crate::utils::style_action(build_mode),
            crate::utils::style_path(&self.base.build_source)
        ));

        crate::utils::print_message(&format!(
            "Running apptainer build from {}",
            crate::utils::style_path(&self.base.build_source)
        ));

        // Get apptainer instance from config
        let config = crate::config::AppConfig::new(false, false)?;
        let apptainer = crate::apptainer::Apptainer::from_config(&config);

        // Build SIF using apptainer build --fakeroot
        // Note: tmp_overlay_path is a SIF file at this stage
        let build_opts = crate::apptainer::build::BuildOptions::default();

        if let Err(e) = apptainer.build(&self.base.build_source, &self.base.tmp_overlay_path, Some(build_opts)) {
            crate::utils::print_error(&format!(
                "Apptainer build failed for {}: {}",
                styled_overlay,
                e
            ));
            self.cleanup(true)?;
            anyhow::bail!("Failed to build SIF from {}: {}", self.base.build_source, e);
        }

        crate::utils::print_message(&format!(
            "Extracting SquashFS to {}",
            crate::utils::style_path(&target_path.to_string_lossy())
        ));

        // Extract SquashFS from SIF
        if let Err(e) = apptainer.dump_sif_to_squashfs(&self.base.tmp_overlay_path, &target_path) {
            crate::utils::print_error(&format!(
                "Failed to export SquashFS for {}: {}",
                styled_overlay,
                e
            ));
            self.cleanup(true)?;
            anyhow::bail!("Failed to dump SquashFS from SIF: {}", e);
        }

        // Set permissions on target overlay
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            if let Err(e) = std::fs::set_permissions(&target_path, std::fs::Permissions::from_mode(0o664)) {
                crate::utils::print_debug(&format!(
                    "Failed to set permissions on {}: {}",
                    target_path.display(),
                    e
                ));
            }
        }

        crate::utils::print_message(&format!("Finished overlay {}", styled_overlay));
        crate::utils::print_debug(&format!(
            "Overlay {} created at {}. Removing temporary overlay...",
            styled_overlay,
            crate::utils::style_path(&target_path.to_string_lossy())
        ));

        // Cleanup temp files
        self.cleanup(false)?;

        Ok(())
    }

    fn get_missing_dependencies(&self) -> Result<Vec<String>> {
        self.base.get_missing_dependencies()
    }

    fn create_tmp_overlay(&self, force: bool) -> Result<()> {
        self.base.create_tmp_overlay(force)
    }

    fn cleanup(&self, failed: bool) -> Result<()> {
        self.base.cleanup(failed)
    }
}
