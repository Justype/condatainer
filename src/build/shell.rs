//! Shell script-based build object implementation

use anyhow::{Context, Result};
use std::path::{Path, PathBuf};

use super::object::{BaseBuildObject, BuildObject, BuildType, ScriptSpecs};

/// Shell script-based build object
#[derive(Debug, Clone)]
pub struct ShellBuildObject {
    pub base: BaseBuildObject,
    pub is_ref: bool, // Whether this is a reference overlay (large datasets)
}

impl ShellBuildObject {
    pub fn new(
        name_version: String,
        build_source: String,
        images_dir: &Path,
        tmp_dir: &Path,
        is_ref: bool,
    ) -> Self {
        let base = BaseBuildObject::new(name_version, build_source, images_dir, tmp_dir);
        Self { base, is_ref }
    }

    /// Create a SquashFS file from the given source directory
    fn create_squashfs(
        &self,
        source_dir: &str,
        target_path: &Path,
        apptainer: &crate::apptainer::Apptainer,
        base_image: &str,
        compress_args: &str,
    ) -> Result<()> {
        let target_abs = target_path.canonicalize()
            .unwrap_or_else(|_| target_path.to_path_buf());

        // Build mksquashfs command
        let bash_script = if source_dir == "/cnt" {
            // For app overlays: set permissions first
            format!(
                r#"
set -euo pipefail

# Trap signals for clean exit
trap 'exit 130' INT TERM

echo "Setting permissions..."
find /cnt -type f -exec chmod ug+rw,o+r {{}} \;
find /cnt -type d -exec chmod ug+rwx,o+rx {{}} \;
mksquashfs /cnt {} -processors {} {} -b 1M
"#,
                target_abs.to_string_lossy(),
                self.base.ncpus,
                compress_args
            )
        } else {
            // For ref overlays: just pack the directory
            format!(
                r#"
set -euo pipefail

# Trap signals for clean exit
trap 'exit 130' INT TERM

mksquashfs {} {} -processors {} {} -b 1M
"#,
                source_dir,
                target_abs.to_string_lossy(),
                self.base.ncpus,
                compress_args
            )
        };

        let config = crate::config::AppConfig::new(false, false)?;
        
        let tmp_overlay_str = self.base.tmp_overlay_path.to_string_lossy().to_string();
        let opts = crate::apptainer::exec::ExecOptions::new()
            .overlay(&tmp_overlay_str);

        crate::utils::print_message(&format!(
            "Packing {} into {}",
            crate::utils::style_path(source_dir),
            crate::utils::style_path(&target_abs.to_string_lossy())
        ));

        if let Err(e) = apptainer.exec(base_image, &["/bin/bash", "-c", &bash_script], Some(opts)) {
            anyhow::bail!("Failed to create SquashFS: {}", e);
        }

        Ok(())
    }
}

impl BuildObject for ShellBuildObject {
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
        BuildType::Shell
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

    fn build(&self, build_deps: bool) -> Result<()> {
        let target_path = self.base.target_overlay_path.canonicalize()
            .unwrap_or_else(|_| self.base.target_overlay_path.clone());
        let styled_overlay = crate::utils::style_name(&self.name_version());

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
            "Building overlay {} ({} build)",
            styled_overlay,
            crate::utils::style_action(build_mode)
        ));

        // Create temporary overlay
        crate::utils::print_message(&format!(
            "Creating temporary overlay {}",
            crate::utils::style_path(&self.base.tmp_overlay_path.to_string_lossy())
        ));
        self.create_tmp_overlay(false)?;

        crate::utils::print_message(&format!(
            "Populating overlay {} via {}",
            styled_overlay,
            crate::utils::style_command(&self.base.build_source)
        ));

        // Check dependencies
        let missing_deps = self.get_missing_dependencies()?;
        if !missing_deps.is_empty() {
            let dep_list = missing_deps.join(", ");
            crate::utils::print_error(&format!(
                "Missing dependencies for {}: {}",
                styled_overlay,
                crate::utils::style_number(&dep_list)
            ));
            
            if build_deps {
                anyhow::bail!(
                    "Missing dependencies: {} (recursive build not yet implemented)",
                    dep_list
                );
            }
            
            anyhow::bail!(
                "Missing dependencies for {}: {}. Please install them first",
                self.name_version(),
                dep_list
            );
        }

        // Parse name and version
        let parts: Vec<&str> = self.name_version().split('/').collect();
        let name = parts[0];
        let version = if parts.len() >= 2 { parts[1] } else { "env" };

        // Get config values
        let config = crate::config::AppConfig::new(false, false)?;
        let base_image = &config.base_image;
        let apptainer = crate::apptainer::Apptainer::from_config(&config);
        let condatainer_dir = &config.base_dir;
        let compress_args = &config.build.compress_args;

        // Build environment settings
        let tmp_overlay_str = self.base.tmp_overlay_path.to_string_lossy().to_string();
        let mut opts = crate::apptainer::exec::ExecOptions::new()
            .env(&format!("app_name={}", name))
            .env(&format!("version={}", version))
            .env(&format!("app_name_version={}", self.name_version()))
            .env("tmp_dir=/ext3/tmp")
            .env("TMPDIR=/ext3/tmp")
            .env(&format!("SLURM_CPUS_PER_TASK={}", self.base.ncpus))
            .env("IN_CONDATAINER=1")
            .overlay(&tmp_overlay_str);

        // Set target_dir based on ref type
        let (source_dir, target_dir) = if self.is_ref {
            // Create cnt_dir for ref overlays (large reference datasets)
            std::fs::create_dir_all(&self.base.cnt_dir_path)
                .with_context(|| format!("Failed to create cnt_dir: {}", self.base.cnt_dir_path.display()))?;
            
            let target = self.base.cnt_dir_path.join(self.name_version());
            opts = opts.env(&format!("target_dir={}", target.to_string_lossy()));
            
            (self.base.cnt_dir_path.to_string_lossy().to_string(), target)
        } else {
            let target = format!("/cnt/{}", self.name_version());
            opts = opts.env(&format!("target_dir={}", target));
            ("/cnt".to_string(), PathBuf::from(target))
        };

        // Build bash script
        let bash_script = format!(
            r#"
set -euo pipefail

# Trap signals for clean exit
trap 'exit 130' INT TERM

mkdir -p $TMPDIR
bash {}
if [ $? -ne 0 ]; then
    echo "Build script {} failed."
    exit 1
fi
"#,
            self.base.build_source,
            self.base.build_source
        );

        crate::utils::print_message(&format!(
            "Executing build script {} for overlay {}",
            crate::utils::style_command(&self.base.build_source),
            styled_overlay
        ));

        // Set up signal handler to ensure cleanup on interrupt
        let interrupted = std::sync::Arc::new(std::sync::atomic::AtomicBool::new(false));
        let interrupted_clone = interrupted.clone();
        let _ = ctrlc::set_handler(move || {
            interrupted_clone.store(true, std::sync::atomic::Ordering::SeqCst);
        });

        // Execute the build script
        let exec_result = apptainer.exec(base_image, &["/bin/bash", "-c", &bash_script], Some(opts));
        
        // Check if interrupted during execution
        if interrupted.load(std::sync::atomic::Ordering::SeqCst) {
            crate::utils::print_warning("Build interrupted by user.");
            crate::utils::print_warning(&format!("Cleaning up temporary files for {}...", styled_overlay));
            self.cleanup(true)?;
            anyhow::bail!("Build interrupted");
        }
        
        // Check for other errors
        if let Err(e) = exec_result {
            crate::utils::print_error(&format!(
                "Build script {} failed for overlay {}: {}",
                crate::utils::style_command(&self.base.build_source),
                styled_overlay,
                e
            ));
            crate::utils::print_warning(&format!("Cleaning up temporary files for {}...", styled_overlay));
            self.cleanup(true)?;
            anyhow::bail!("Build script {} failed: {}", self.base.build_source, e);
        }

        // Create SquashFS
        if self.is_ref {
            // For ref overlays: check files exist, then pack from cnt_dir
            let entries: Vec<_> = std::fs::read_dir(&target_dir)
                .with_context(|| format!("Failed to read target directory: {}", target_dir.display()))?
                .collect::<std::io::Result<_>>()?;
            
            if entries.is_empty() {
                crate::utils::print_error(&format!(
                    "No files produced for overlay {}",
                    styled_overlay
                ));
                self.cleanup(true)?;
                anyhow::bail!(
                    "Overlay build script did not create any files in {}",
                    target_dir.display()
                );
            }

            crate::utils::print_message(&format!(
                "Creating SquashFS from {} for overlay {}",
                crate::utils::style_path(&self.base.cnt_dir_path.to_string_lossy()),
                styled_overlay
            ));
            
            let base_image_str = base_image.to_string_lossy();
            self.create_squashfs(&source_dir, &target_path, &apptainer, &base_image_str, compress_args)?;
        } else {
            // For app overlays: pack from /cnt
            crate::utils::print_message(&format!(
                "Preparing SquashFS from /cnt for overlay {}",
                styled_overlay
            ));
            
            let base_image_str = base_image.to_string_lossy();
            self.create_squashfs(&source_dir, &target_path, &apptainer, &base_image_str, compress_args)?;
        }

        // Set permissions
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            std::fs::set_permissions(&target_path, std::fs::Permissions::from_mode(0o664))?;
        }

        crate::utils::print_message(&format!("Finished overlay {}", styled_overlay));
        crate::utils::print_debug(&format!(
            "Overlay created at {}. Cleaning up...",
            crate::utils::style_path(&target_path.to_string_lossy())
        ));

        // Cleanup
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
