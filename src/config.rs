use anyhow::{Context, Result};
use config as config_crate;
use config_crate::{Config as ConfigLoader, Environment, File};
use serde::{Deserialize, Serialize};
use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Duration;

use crate::utils;

pub const VERSION: &str = "1.0.5";
pub const GITHUB_REPO: &str = "Justype/condatainer";
pub const REMOTE_METADATA_SUFFIX: &str = "metadata/build-scripts.json.gz";
pub const REMOTE_HELPER_SUFFIX: &str = "metadata/helper-scripts.json.gz";
pub const MANUAL_URL: &str =
    "https://github.com/Justype/condatainer/blob/main/docs/manuals/condatainer.md";
pub const PREBUILT_OVERLAYS: [&str; 1] = ["build-essential"];
pub const DEFAULT_COMPRESS_ARGS: &str = "-comp lz4";
pub const MODERN_COMPRESS_ARGS: &str = "-comp zstd -Xcompression-level 8";
pub const CONFIG_FILENAME: &str = "config";
pub const CONFIG_TYPE: &str = "yaml";

#[allow(dead_code)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BuildConfig {
    pub default_cpus: u32,
    pub default_memory_mb: u64,
    #[serde(skip)]
    pub default_time: Duration,
    pub tmp_size_mb: u32,
    pub compress_args: String,
    pub overlay_type: String,
}

#[allow(dead_code)]
#[derive(Debug)]
pub struct AppConfig {
    pub debug: bool,
    pub submit_job: bool,
    pub version: &'static str,
    pub program_dir: PathBuf,
    pub program_base: PathBuf,  // Installation directory (may be read-only)
    pub base_dir: PathBuf,       // User data directory (writable)
    pub build_scripts_dir: PathBuf,
    pub images_dir: PathBuf,
    pub helper_scripts_dir: PathBuf,
    pub tmp_dir: PathBuf,
    pub logs_dir: PathBuf,
    pub apptainer_bin: PathBuf,
    pub base_image: PathBuf,
    pub scheduler_bin: Option<PathBuf>,
    pub build: BuildConfig,
}

/// Serializable config structure for saving to file
#[derive(Debug, Serialize, Deserialize)]
struct SaveableConfig {
    base_dir: Option<String>,
    apptainer_bin: Option<String>,
    scheduler_bin: Option<String>,
    scheduler_type: Option<String>,
    logs_dir: Option<String>,
    submit_job: bool,
    build: SaveableBuildConfig,
}

#[derive(Debug, Serialize, Deserialize)]
struct SaveableBuildConfig {
    default_cpus: u32,
    default_mem_mb: u64,
    default_time: String,
    tmp_size_mb: u32,
    compress_args: String,
    overlay_type: String,
}

impl AppConfig {
    pub fn new(debug: bool, submit_job: bool) -> Result<Self> {
        let mut builder = ConfigLoader::builder()
            .set_default("debug", false)?
            .set_default("submit_job", true)?
            .set_default("tmp_size_mb", 20480)?
            .set_default("default_cpus", 4)?
            .set_default("default_memory_mb", 8192)?
            .set_default("default_time_seconds", 2 * 60 * 60)?
            .set_default("overlay_type", "ext3")?;

        // Add system config file
        builder = builder.add_source(File::with_name("/etc/condatainer/config").required(false));

        // Add user config files (with priority)
        if let Some(config_dir) = dirs::config_dir() {
            let user_conf = config_dir.join("condatainer").join("config");
            if let Some(path_str) = user_conf.to_str() {
                builder = builder.add_source(File::with_name(path_str).required(false));
            }
        }

        if let Some(home) = env::var_os("HOME") {
            let user_conf = PathBuf::from(&home)
                .join(".condatainer")
                .join("config");
            if let Some(path_str) = user_conf.to_str() {
                builder = builder.add_source(File::with_name(path_str).required(false));
            }
        }

        // Add environment variables (highest priority)
        builder = builder
            .add_source(Environment::with_prefix("CONDATAINER").separator("_"))
            .set_override("debug", debug)?
            .set_override("submit_job", submit_job)?;

        if let Ok(base_dir) = env::var("CONDATAINER_BASE_DIR") {
            builder = builder.set_override("base_dir", base_dir)?;
        }

        let settings: Settings = builder.build()?.try_deserialize()?;

        let exe_path = env::current_exe().context("Failed to locate the CondaTainer executable")?;
        let program_dir = exe_path
            .parent()
            .map(PathBuf::from)
            .unwrap_or_else(|| PathBuf::from("."));

        utils::print_debug(&format!("Debug mode: {:#?}", settings));

        // Detect program installation base
        let program_base = detect_program_base();
        
        let base_dir = settings
            .base_dir
            .map(|p| PathBuf::from(shellexpand::full(&p)
                .unwrap_or_else(|_| std::borrow::Cow::Borrowed(&p))
                .to_string()))
            .unwrap_or_else(|| detect_base_dir(&program_base));
        fs::create_dir_all(&base_dir).context("Failed to create base directory")?;

        utils::print_debug(&format!("Program base: {}", program_base.display()));
        utils::print_debug(&format!("Using base directory: {}", base_dir.display()));

        let build_scripts_dir = base_dir.join("build-scripts");
        let images_dir = base_dir.join("images");
        let helper_scripts_dir = base_dir.join("helper-scripts");
        let tmp_dir = base_dir.join("tmp");
        
        let logs_dir = settings.logs_dir
            .map(|p| PathBuf::from(shellexpand::full(&p)
                .unwrap_or_else(|_| std::borrow::Cow::Borrowed(&p))
                .to_string()))
            .or_else(|| env::var("HOME").ok().map(|h| PathBuf::from(h).join("logs")))
            .unwrap_or_else(|| PathBuf::from("logs"));

        let apptainer_bin = settings
            .apptainer_bin
            .map(PathBuf::from)
            .unwrap_or_else(|| detect_apptainer_bin());
        
        let scheduler_bin = settings.scheduler_bin.map(PathBuf::from);
        
        // Determine if we should submit jobs
        let submit_job = if settings.submit_job {
            // Auto-disable if no scheduler binary is available
            if let Some(ref bin) = scheduler_bin {
                validate_binary(&bin.to_string_lossy())
            } else {
                detect_scheduler_bin().is_some()
            }
        } else {
            false
        };

        let compress_args = auto_detect_compression(&apptainer_bin, settings.compress_args);

        Ok(Self {
            debug: settings.debug,
            submit_job,
            version: VERSION,
            program_dir,
            program_base: program_base.clone(),
            base_dir: base_dir.clone(),
            build_scripts_dir,
            images_dir,
            helper_scripts_dir,
            tmp_dir,
            logs_dir,
            apptainer_bin,
            base_image: base_dir.join("images").join("base_image.sif"),
            scheduler_bin,
            build: BuildConfig {
                default_cpus: settings.default_cpus,
                default_memory_mb: settings.default_memory_mb,
                default_time: Duration::from_secs(settings.default_time_seconds),
                tmp_size_mb: settings.tmp_size_mb,
                compress_args,
                overlay_type: settings.overlay_type,
            },
        })
    }

    pub fn metadata_url(&self) -> String {
        format!(
            "https://raw.githubusercontent.com/{}/{}",
            GITHUB_REPO, REMOTE_METADATA_SUFFIX
        )
    }

    pub fn helper_metadata_url(&self) -> String {
        format!(
            "https://raw.githubusercontent.com/{}/{}",
            GITHUB_REPO, REMOTE_HELPER_SUFFIX
        )
    }

    pub fn manual_url(&self) -> &'static str {
        MANUAL_URL
    }

    pub fn base_image_def(&self) -> PathBuf {
        self.build_scripts_dir.join("base_image.def")
    }

    #[allow(dead_code)]
    pub fn base_image_sif(&self) -> PathBuf {
        self.base_image.clone()
    }

    /// Get user config path (~/.config/condatainer/config.yaml or ~/.condatainer/config.yaml)
    pub fn get_user_config_path() -> Result<PathBuf> {
        if let Some(config_dir) = dirs::config_dir() {
            Ok(config_dir.join("condatainer").join(format!("{}.{}", CONFIG_FILENAME, CONFIG_TYPE)))
        } else if let Some(home) = dirs::home_dir() {
            Ok(home.join(".condatainer").join(format!("{}.{}", CONFIG_FILENAME, CONFIG_TYPE)))
        } else {
            anyhow::bail!("Cannot determine user config directory")
        }
    }

    /// Save current config to user config file
    pub fn save_config(&self) -> Result<()> {
        let config_path = Self::get_user_config_path()?;
        let config_dir = config_path.parent()
            .context("Invalid config path")?;
        
        fs::create_dir_all(config_dir)
            .context("Failed to create config directory")?;

        // Convert to saveable format
        let saveable = SaveableConfig {
            base_dir: Some(self.base_dir.to_string_lossy().to_string()),
            apptainer_bin: Some(self.apptainer_bin.to_string_lossy().to_string()),
            scheduler_bin: self.scheduler_bin.as_ref()
                .map(|p| p.to_string_lossy().to_string()),
            scheduler_type: None, // Could be enhanced to store scheduler type
            logs_dir: Some(self.logs_dir.to_string_lossy().to_string()),
            submit_job: self.submit_job,
            build: SaveableBuildConfig {
                default_cpus: self.build.default_cpus,
                default_mem_mb: self.build.default_memory_mb,
                default_time: format!("{}s", self.build.default_time.as_secs()),
                tmp_size_mb: self.build.tmp_size_mb,
                compress_args: self.build.compress_args.clone(),
                overlay_type: self.build.overlay_type.clone(),
            },
        };

        let config_content = serde_yaml::to_string(&saveable)
            .context("Failed to serialize config")?;
        
        fs::write(&config_path, config_content)
            .context("Failed to write config file")?;

        Ok(())
    }
}

#[derive(Debug, Deserialize)]
#[serde(default)]
struct Settings {
    debug: bool,
    submit_job: bool,
    base_dir: Option<String>,
    apptainer_bin: Option<String>,
    scheduler_bin: Option<String>,
    scheduler_type: Option<String>,
    logs_dir: Option<String>,
    tmp_size_mb: u32,
    default_cpus: u32,
    default_memory_mb: u64,
    default_time_seconds: u64,
    compress_args: Option<String>,
    overlay_type: String,
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            debug: false,
            submit_job: true,
            base_dir: None,
            apptainer_bin: None,
            scheduler_bin: None,
            scheduler_type: None,
            logs_dir: None,
            tmp_size_mb: 20480,
            default_cpus: 4,
            default_memory_mb: 8192,
            default_time_seconds: 2 * 60 * 60,
            compress_args: None,
            overlay_type: "ext3".to_string(),
        }
    }
}

/// Detect the program installation base directory
fn detect_program_base() -> PathBuf {
    if let Ok(exe_path) = env::current_exe() {
        if let Some(exe_dir) = exe_path.parent() {
            // Check if executable is in a 'bin' directory
            if exe_dir.file_name().and_then(|n| n.to_str()) == Some("bin") {
                // Return the parent of the bin directory
                if let Some(parent) = exe_dir.parent() {
                    return parent.to_path_buf();
                }
            }
            // If not in bin, use the executable's directory
            return exe_dir.to_path_buf();
        }
    }
    // Fallback to current directory
    env::current_dir().unwrap_or_else(|_| PathBuf::from("."))
}

/// Check if a directory is writable
pub fn is_writable(path: &Path) -> bool {
    // Try to create the directory if it doesn't exist
    if !path.exists() {
        if let Ok(()) = fs::create_dir_all(path) {
            // Successfully created, it's writable
            return true;
        }
        return false;
    }
    
    // Directory exists, try to write a test file
    let test_file = path.join(".condatainer_write_test");
    let writable = fs::write(&test_file, b"test").is_ok();
    if writable {
        let _ = fs::remove_file(&test_file);
    }
    writable
}

/// Detect base directory for user data (must be writable)
fn detect_base_dir(program_base: &Path) -> PathBuf {
    // First, check if program_base is writable and not in $HOME/bin or $HOME/.local/bin
    if let Ok(home) = env::var("HOME") {
        let home_path = Path::new(&home);
        let home_bin = home_path.join("bin");
        let home_local_bin = home_path.join(".local").join("bin");
        let program_bin = program_base.join("bin");
        
        // If program is NOT in $HOME/bin or $HOME/.local/bin, try using program_base
        if program_bin != home_bin && program_bin != home_local_bin {
            if is_writable(program_base) {
                return program_base.to_path_buf();
            }
        }
    }
    
    // Fall back to normal user directory detection
    detect_base_dir_normal()
}

fn detect_base_dir_normal() -> PathBuf {
    if let Ok(scratch) = env::var("SCRATCH") {
        return Path::new(&scratch).join("condatainer");
    }
    if let Ok(home) = env::var("HOME") {
        return Path::new(&home).join("condatainer");
    }
    env::current_dir()
        .map(|cwd| cwd.join("condatainer"))
        .unwrap_or_else(|_| PathBuf::from("condatainer"))
}

pub fn is_inside_container() -> bool {
    if env::var("IN_CONDATAINER").is_ok() {
        return true;
    }
    if env::var("APPTAINER_NAME").is_ok() || env::var("SINGULARITY_NAME").is_ok() {
        return true;
    }
    if Path::new("/.singularity.d").exists() || Path::new("/.apptainer.d").exists() {
        return true;
    }
    false
}

fn detect_apptainer_bin() -> PathBuf {
    if is_inside_container() {
        const CONTAINER_PATHS: [&str; 3] = [
            "/usr/bin/apptainer",
            "/usr/local/bin/apptainer",
            "/opt/apptainer/bin/apptainer",
        ];
        for path in CONTAINER_PATHS {
            let candidate = Path::new(path);
            if candidate.exists() {
                return candidate.to_path_buf();
            }
        }
    }
    if let Some(path) = which_in_path("apptainer") {
        return path;
    }
    PathBuf::from("apptainer")
}

fn which_in_path(bin: &str) -> Option<PathBuf> {
    env::var_os("PATH").and_then(|paths| {
        env::split_paths(&paths).find_map(|dir| {
            let candidate = dir.join(bin);
            if candidate.exists() {
                Some(candidate)
            } else {
                None
            }
        })
    })
}

fn detect_apptainer_version(apptainer_bin: &Path) -> Option<(u64, u64, u64)> {
    let output = Command::new(apptainer_bin).arg("--version").output().ok()?;
    if !output.status.success() {
        return None;
    }
    let stdout = String::from_utf8(output.stdout).ok()?;
    parse_version(&stdout)
}

fn parse_version(text: &str) -> Option<(u64, u64, u64)> {
    for chunk in text.split(|c: char| !c.is_digit(10) && c != '.') {
        if chunk.is_empty() {
            continue;
        }
        let parts: Vec<&str> = chunk.split('.').collect();
        if parts.len() < 2 {
            continue;
        }
        let mut nums = Vec::new();
        for part in parts.iter().take(3) {
            if part.is_empty() {
                break;
            }
            if let Ok(n) = part.parse::<u64>() {
                nums.push(n);
            } else {
                nums.clear();
                break;
            }
        }
        if nums.len() < 2 {
            continue;
        }
        while nums.len() < 3 {
            nums.push(0);
        }
        return Some((nums[0], nums[1], nums[2]));
    }
    None
}
/// Validate if a binary exists and is executable
pub fn validate_binary(bin_path: &str) -> bool {
    if bin_path.is_empty() {
        return false;
    }

    let path = Path::new(bin_path);
    
    // If it's an absolute path, check directly
    if path.is_absolute() {
        if let Ok(metadata) = fs::metadata(path) {
            #[cfg(unix)]
            {
                use std::os::unix::fs::PermissionsExt;
                return metadata.permissions().mode() & 0o111 != 0;
            }
            #[cfg(not(unix))]
            return metadata.is_file();
        }
        return false;
    }

    // Otherwise, try to find it in PATH
    which_in_path(bin_path).is_some()
}

/// Detect scheduler binary and type
/// Returns (binary_path, scheduler_type) if found
pub fn detect_scheduler_bin() -> Option<(PathBuf, String)> {
    // Try SLURM first (most common in HPC)
    if let Some(path) = which_in_path("sbatch") {
        return Some((path, "SLURM".to_string()));
    }

    // Try PBS/Torque
    if let Some(path) = which_in_path("qsub") {
        // Check if it's SGE by looking for SGE_ROOT
        if env::var("SGE_ROOT").is_ok() {
            return Some((path, "SGE".to_string()));
        }
        return Some((path, "PBS".to_string()));
    }

    // Try LSF
    if let Some(path) = which_in_path("bsub") {
        return Some((path, "LSF".to_string()));
    }

    None
}

/// Parse duration from string (e.g., "2h", "30m", "1h30m", or "02:00:00")
pub fn parse_duration(s: &str) -> Result<Duration> {
    // Try to parse as simple duration first (e.g., "2h", "30m")
    if let Some(stripped) = s.strip_suffix('h') {
        let hours: u64 = stripped.parse()
            .context("Invalid hours value")?;
        return Ok(Duration::from_secs(hours * 3600));
    }
    
    if let Some(stripped) = s.strip_suffix('m') {
        let minutes: u64 = stripped.parse()
            .context("Invalid minutes value")?;
        return Ok(Duration::from_secs(minutes * 60));
    }
    
    if let Some(stripped) = s.strip_suffix('s') {
        let seconds: u64 = stripped.parse()
            .context("Invalid seconds value")?;
        return Ok(Duration::from_secs(seconds));
    }

    // Try to parse as HH:MM:SS format
    if s.contains(':') {
        let parts: Vec<&str> = s.split(':').collect();
        if parts.len() == 3 {
            let hours: u64 = parts[0].parse()
                .context("Invalid hours in HH:MM:SS format")?;
            let minutes: u64 = parts[1].parse()
                .context("Invalid minutes in HH:MM:SS format")?;
            let seconds: u64 = parts[2].parse()
                .context("Invalid seconds in HH:MM:SS format")?;
            return Ok(Duration::from_secs(hours * 3600 + minutes * 60 + seconds));
        }
    }

    // Try combined format (e.g., "1h30m")
    let mut total_seconds = 0u64;
    let mut current_num = String::new();
    
    for ch in s.chars() {
        if ch.is_digit(10) {
            current_num.push(ch);
        } else if !current_num.is_empty() {
            let num: u64 = current_num.parse()
                .context("Invalid number in duration")?;
            match ch {
                'h' => total_seconds += num * 3600,
                'm' => total_seconds += num * 60,
                's' => total_seconds += num,
                _ => anyhow::bail!("Invalid duration unit: {}", ch),
            }
            current_num.clear();
        }
    }

    if total_seconds > 0 {
        Ok(Duration::from_secs(total_seconds))
    } else {
        anyhow::bail!("Invalid duration format: {}", s)
    }
}

/// Auto-detect binaries and update config if needed
pub fn auto_detect_and_save(config: &mut AppConfig) -> Result<bool> {
    let mut updated = false;

    // Auto-detect apptainer if not inside container
    if !is_inside_container() {
        if !validate_binary(&config.apptainer_bin.to_string_lossy()) {
            if let Some(detected) = which_in_path("apptainer") {
                config.apptainer_bin = detected;
                updated = true;
            }
        }
    }

    // Auto-detect scheduler
    if config.scheduler_bin.is_none() 
        || !validate_binary(&config.scheduler_bin.as_ref().unwrap().to_string_lossy()) {
        if let Some((bin, _type)) = detect_scheduler_bin() {
            config.scheduler_bin = Some(bin);
            updated = true;
        }
    }

    if updated {
        config.save_config()?;
    }

    Ok(updated)
}

/// Force re-detection of binaries from current environment and save
pub fn force_detect_and_save(config: &mut AppConfig) -> Result<bool> {
    let mut updated = false;

    // Always re-detect apptainer if not inside container
    if !is_inside_container() {
        if let Some(detected) = which_in_path("apptainer") {
            if config.apptainer_bin != detected {
                config.apptainer_bin = detected;
                updated = true;
            }
        }
    }

    // Always re-detect scheduler
    if let Some((bin, _type)) = detect_scheduler_bin() {
        let needs_update = match &config.scheduler_bin {
            Some(current) => current != &bin,
            None => true,
        };
        if needs_update {
            config.scheduler_bin = Some(bin);
            updated = true;
        }
    }

    // Always save (even if nothing changed, to create the file)
    config.save_config()?;

    Ok(updated)
}

/// Auto-detect compression based on apptainer version
pub fn auto_detect_compression(apptainer_bin: &Path, explicit_setting: Option<String>) -> String {
    // If user explicitly set compress_args, respect it
    if let Some(setting) = explicit_setting {
        if !setting.is_empty() {
            return setting;
        }
    }

    // Auto-detect based on version
    if let Some((major, minor, _)) = detect_apptainer_version(apptainer_bin) {
        if major > 1 || (major == 1 && minor >= 4) {
            return MODERN_COMPRESS_ARGS.to_string();
        }
    }
    
    DEFAULT_COMPRESS_ARGS.to_string()
}