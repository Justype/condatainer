//! Conda-based build object implementation

use anyhow::{Context, Result};
use std::path::Path;

use super::object::{BaseBuildObject, BuildObject, BuildType, ScriptSpecs};

/// Conda-based build object
#[derive(Debug, Clone)]
pub struct CondaBuildObject {
    pub base: BaseBuildObject,
}

impl CondaBuildObject {
    pub fn new(name_version: String, build_source: String, images_dir: &Path, tmp_dir: &Path) -> Self {
        let base = BaseBuildObject::new(name_version, build_source, images_dir, tmp_dir);
        Self { base }
    }
}

impl BuildObject for CondaBuildObject {
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
        BuildType::Conda
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
        let styled_overlay = crate::utils::style_name(&self.name_version());

        // Check if already exists
        if target_path.exists() {
            crate::utils::print_debug(&format!(
                "Overlay {} already exists at {}. Skipping creation.",
                styled_overlay,
                crate::utils::style_path(&target_path.to_string_lossy())
            ));
            return Ok(());
        }

        // Create temporary overlay
        self.create_tmp_overlay(false)?;

        // Determine build mode and create micromamba command
        let (install_cmd, bind_paths) = if self.base.build_source.ends_with(".yml") 
            || self.base.build_source.ends_with(".yaml") {
            // Mode 3: YAML file
            let abs_file = Path::new(&self.base.build_source)
                .canonicalize()
                .context("Failed to resolve YAML file path")?;
            let parent_dir = abs_file.parent()
                .ok_or_else(|| anyhow::anyhow!("YAML file has no parent directory"))?;
            
            let cmd = format!(
                "micromamba create -r /ext3/tmp -c conda-forge -c bioconda -q -y -p /cnt/{} -f {}",
                self.name_version(),
                abs_file.to_string_lossy()
            );
            
            crate::utils::print_debug(&format!(
                "Building SquashFS overlay at {} from YAML file {}",
                crate::utils::style_path(&target_path.to_string_lossy()),
                crate::utils::style_path(&abs_file.to_string_lossy())
            ));
            
            (cmd, vec![parent_dir.to_path_buf()])
        } else if !self.base.build_source.is_empty() {
            // Mode 2: Multiple packages
            let packages: Vec<String> = self.base.build_source
                .split(',')
                .map(|s| s.trim().replace('/', "="))
                .collect();
            
            let cmd = format!(
                "micromamba create -r /ext3/tmp -c conda-forge -c bioconda -q -y -p /cnt/{} {}",
                self.name_version(),
                packages.join(" ")
            );
            
            crate::utils::print_debug(&format!(
                "Building SquashFS overlay at {} with packages: {}",
                crate::utils::style_path(&target_path.to_string_lossy()),
                packages.join(", ")
            ));
            
            (cmd, Vec::new())
        } else {
            // Mode 1: Single package (name/version)
            let parts: Vec<&str> = self.name_version().split('/').collect();
            if parts.len() != 2 {
                anyhow::bail!("Invalid conda package format: {}", self.name_version());
            }
            
            let cmd = format!(
                "micromamba create -r /ext3/tmp -c conda-forge -c bioconda -q -y -p /cnt/{} {}={}",
                self.name_version(),
                parts[0],
                parts[1]
            );
            
            crate::utils::print_debug(&format!(
                "Building SquashFS overlay at {} with package: {}={}",
                crate::utils::style_path(&target_path.to_string_lossy()),
                parts[0],
                parts[1]
            ));
            
            (cmd, Vec::new())
        };

        // Get config values
        let config = crate::config::AppConfig::new(false, false)?;
        let base_image = &config.base_image;
        let apptainer = crate::apptainer::Apptainer::from_config(&config);
        let compress_args = &config.build.compress_args;

        // Build bash script
        let bash_script = format!(
            r#"
set -euo pipefail

# Trap signals for clean exit
trap 'exit 130' INT TERM

mkdir -p $TMPDIR
echo "Creating conda environment in overlay..."
{}

echo "Setting permissions..."
find /cnt -type f -exec chmod ug+rw,o+r {{}} \;
find /cnt -type d -exec chmod ug+rwx,o+rx {{}} \;

echo "Packing overlay to SquashFS..."
mksquashfs /cnt {} -processors {} {} -b 1M
"#,
            install_cmd,
            target_path.to_string_lossy(),
            self.base.ncpus,
            compress_args
        );

        // Build exec command
        let tmp_overlay_str = self.base.tmp_overlay_path.to_string_lossy().to_string();
        let mut opts = crate::apptainer::exec::ExecOptions::new()
            .env("TMPDIR=/ext3/tmp")
            .overlay(&tmp_overlay_str)
            .fakeroot();

        // Add bind paths for YAML files
        for path in &bind_paths {
            let path_str = path.to_string_lossy().to_string();
            opts = opts.bind(&path_str);
        }

        crate::utils::print_message(&format!(
            "Building Conda overlay: {}",
            styled_overlay
        ));

        // Set up signal handler to ensure cleanup on interrupt
        let interrupted = std::sync::Arc::new(std::sync::atomic::AtomicBool::new(false));
        let interrupted_clone = interrupted.clone();
        let _ = ctrlc::set_handler(move || {
            interrupted_clone.store(true, std::sync::atomic::Ordering::SeqCst);
        });

        // Execute the build
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
            crate::utils::print_warning(&format!("Build failed for {}, cleaning up temporary files...", styled_overlay));
            self.cleanup(true)?;
            anyhow::bail!("Failed to build conda package {}: {}", self.name_version(), e);
        }

        // Set permissions on target overlay
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            std::fs::set_permissions(&target_path, std::fs::Permissions::from_mode(0o664))?;
        }

        crate::utils::print_debug(&format!(
            "Overlay {} created at {}. Removing temporary overlay...",
            styled_overlay,
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
