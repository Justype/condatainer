//! Overlay image creation and management

use anyhow::{Context, Result};
use std::path::Path;
use std::process::Command;

use crate::utils;

/// Profile for filesystem tuning
#[derive(Debug, Clone, Copy)]
pub enum Profile {
    /// Small files (Conda packages, Python site-packages) - 1 inode per 4KB
    Small,
    /// Default profile - 1 inode per 16KB
    Default,
    /// Large files (genomes, FASTA, databases) - 1 inode per 1MB
    Large,
}

impl Profile {
    fn block_size(&self) -> u32 {
        4096
    }

    fn inode_ratio(&self) -> u32 {
        match self {
            Profile::Small => 4096,
            Profile::Default => 16384,
            Profile::Large => 1048576,
        }
    }

    fn reserved_percent(&self) -> u8 {
        0
    }
}

/// Options for creating an overlay
#[derive(Debug, Clone)]
pub struct CreateOptions {
    /// Path to the overlay image file
    pub path: String,
    /// Size in megabytes
    pub size_mb: u32,
    /// Owner user ID
    pub uid: u32,
    /// Owner group ID
    pub gid: u32,
    /// Filesystem tuning profile
    pub profile: Profile,
    /// Create sparse file (vs allocated)
    pub sparse: bool,
    /// Filesystem type: "ext3" or "ext4"
    pub fs_type: String,
}

impl Default for CreateOptions {
    fn default() -> Self {
        Self {
            path: String::new(),
            size_mb: 20480, // 20GB
            uid: unsafe { libc::getuid() },
            gid: unsafe { libc::getgid() },
            profile: Profile::Default,
            sparse: true,
            fs_type: "ext3".to_string(),
        }
    }
}

impl CreateOptions {
    pub fn new<P: AsRef<Path>>(path: P, size_mb: u32) -> Self {
        Self {
            path: path.as_ref().to_string_lossy().to_string(),
            size_mb,
            ..Default::default()
        }
    }

    pub fn with_profile(mut self, profile: Profile) -> Self {
        self.profile = profile;
        self
    }

    pub fn with_sparse(mut self, sparse: bool) -> Self {
        self.sparse = sparse;
        self
    }

    pub fn with_fs_type<S: Into<String>>(mut self, fs_type: S) -> Self {
        self.fs_type = fs_type.into();
        self
    }

    pub fn for_root(mut self) -> Self {
        self.uid = 0;
        self.gid = 0;
        self
    }
}

/// Create an ext3/ext4 overlay with specified options
pub fn create_overlay(opts: &CreateOptions) -> Result<()> {
    // Validate filesystem type
    let fs_type = opts.fs_type.to_lowercase();
    if fs_type != "ext3" && fs_type != "ext4" {
        anyhow::bail!("Unsupported filesystem type '{}': must be ext3 or ext4", opts.fs_type);
    }

    // Check for required tools
    check_dependencies(&["dd", "mke2fs"])?;

    let type_str = if opts.sparse {
        if opts.uid == 0 && opts.gid == 0 {
            "Sparse Fakeroot "
        } else {
            "Sparse "
        }
    } else {
        if opts.uid == 0 && opts.gid == 0 {
            "Allocated Fakeroot "
        } else {
            "Allocated "
        }
    };

    utils::print_message(&format!(
        "Creating {}overlay {} ({} MB, {}, block size {}, inode ratio {})",
        type_str,
        utils::style_path(&opts.path),
        utils::style_number(opts.size_mb),
        utils::style_info(&fs_type),
        utils::style_number(opts.profile.block_size()),
        utils::style_number(opts.profile.inode_ratio())
    ));

    // 1. Create file with dd
    create_file_with_dd(&opts.path, opts.size_mb, opts.sparse)?;

    // 2. Create filesystem with mke2fs
    utils::print_message(&format!(
        "Formatting {} as {}",
        utils::style_path(&opts.path),
        utils::style_info(&fs_type)
    ));
    create_filesystem(&opts.path, &fs_type, opts.profile)?;

    // 3. Inject overlay structure (upper/work directories) using debugfs
    utils::print_message(&format!(
        "Injecting overlayFS structure via debugfs for {}",
        utils::style_path(&opts.path)
    ));
    inject_overlay_structure(&opts.path, opts.uid, opts.gid)?;

    // 4. Set file ownership (if running as root or with sudo)
    if opts.uid != 0 || opts.gid != 0 {
        set_ownership(&opts.path, opts.uid, opts.gid)?;
    }

    let type_str = if opts.sparse { "Sparse " } else { "Allocated " }.to_string()
        + if opts.uid == 0 && opts.gid == 0 { "Fakeroot " } else { "" };

    utils::print_success(&format!(
        "Created {}overlay: {}",
        type_str,
        utils::style_path(&opts.path)
    ));

    Ok(())
}

/// Create a file using dd
fn create_file_with_dd(path: &str, size_mb: u32, sparse: bool) -> Result<()> {
    let mut args = vec![
        "if=/dev/zero".to_string(),
        format!("of={}", path),
        "bs=1M".to_string(),
    ];

    if sparse {
        // Sparse: seek to end, write nothing (instant)
        args.push("count=0".to_string());
        args.push(format!("seek={}", size_mb));
    } else {
        // Allocated: write zeros (slow but guarantees space)
        args.push(format!("count={}", size_mb));
        args.push("status=progress".to_string());
    }

    utils::print_debug(&format!("Running: dd {}", args.join(" ")));

    let output = Command::new("dd")
        .args(&args)
        .output()
        .context("Failed to execute dd")?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        anyhow::bail!("dd failed: {}", stderr);
    }

    Ok(())
}

/// Create ext3/ext4 filesystem using mke2fs
fn create_filesystem(path: &str, fs_type: &str, profile: Profile) -> Result<()> {
    let args = vec![
        "-t".to_string(),
        fs_type.to_string(),
        "-b".to_string(),
        profile.block_size().to_string(),
        "-i".to_string(),
        profile.inode_ratio().to_string(),
        "-m".to_string(),
        profile.reserved_percent().to_string(),
        "-F".to_string(), // Force
        path.to_string(),
    ];

    utils::print_debug(&format!("Running: mke2fs {}", args.join(" ")));

    let output = Command::new("mke2fs")
        .args(&args)
        .output()
        .context("Failed to execute mke2fs")?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        anyhow::bail!("mke2fs failed: {}", stderr);
    }

    Ok(())
}

/// Inject overlayFS structure using debugfs
fn inject_overlay_structure(path: &str, uid: u32, gid: u32) -> Result<()> {
    use std::io::Write;
    
    let script = format!(
        "mkdir upper\nmkdir work\nset_inode_field upper uid {}\nset_inode_field upper gid {}\nset_inode_field work uid {}\nset_inode_field work gid {}\nquit\n",
        uid, gid, uid, gid
    );

    utils::print_debug(&format!("Running debugfs script on {}", path));

    let mut child = Command::new("debugfs")
        .arg("-w")
        .arg(path)
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .context("Failed to spawn debugfs")?;

    if let Some(mut stdin) = child.stdin.take() {
        stdin.write_all(script.as_bytes())
            .context("Failed to write to debugfs stdin")?;
    }

    let output = child.wait_with_output()
        .context("Failed to wait for debugfs")?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        let stdout = String::from_utf8_lossy(&output.stdout);
        anyhow::bail!("debugfs failed: {}\n{}", stderr, stdout);
    }

    Ok(())
}

/// Set ownership of the overlay file
fn set_ownership(path: &str, uid: u32, gid: u32) -> Result<()> {
    utils::print_debug(&format!(
        "Setting ownership of {} to {}:{}",
        utils::style_path(path),
        uid,
        gid
    ));

    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        let path_obj = Path::new(path);
        
        // Set permissions to 664 (rw-rw-r--)
        std::fs::set_permissions(path_obj, std::fs::Permissions::from_mode(0o664))?;
        
        // Try to change ownership (may fail if not root)
        unsafe {
            let c_path = std::ffi::CString::new(path)?;
            let result = libc::chown(c_path.as_ptr(), uid, gid);
            if result != 0 {
                utils::print_debug(&format!(
                    "Could not change ownership (may need root): {}",
                    std::io::Error::last_os_error()
                ));
            }
        }
    }

    Ok(())
}

/// Create a temporary ext3 overlay for building
pub fn create_tmp_overlay<P: AsRef<Path>>(path: P, size_mb: u32) -> Result<()> {
    let opts = CreateOptions::new(path, size_mb)
        .with_profile(Profile::Small) // Use small profile for conda packages
        .with_sparse(true)
        .with_fs_type("ext3");

    create_overlay(&opts)
}

/// Create an ext3 overlay for current user
pub fn create_for_current_user<P: AsRef<Path>>(
    path: P,
    size_mb: u32,
    profile: Profile,
    sparse: bool,
) -> Result<()> {
    let opts = CreateOptions::new(path, size_mb)
        .with_profile(profile)
        .with_sparse(sparse)
        .with_fs_type("ext3");

    create_overlay(&opts)
}

/// Create an ext3 overlay for root (fakeroot builds)
pub fn create_for_root<P: AsRef<Path>>(
    path: P,
    size_mb: u32,
    profile: Profile,
    sparse: bool,
) -> Result<()> {
    let opts = CreateOptions::new(path, size_mb)
        .with_profile(profile)
        .with_sparse(sparse)
        .with_fs_type("ext3")
        .for_root();

    create_overlay(&opts)
}

// =============================================================================
// Public API for overlay operations
// =============================================================================

/// Filesystem statistics
#[derive(Debug, Default)]
pub struct FilesystemStats {
    pub block_size: i64,
    pub total_blocks: i64,
    pub free_blocks: i64,
    pub total_inodes: i64,
    pub free_inodes: i64,
    pub filesystem_state: String,
    pub last_mounted: String,
}

impl FilesystemStats {
    /// Calculate used space in bytes and percentage
    pub fn usage(&self) -> (i64, f64) {
        if self.total_blocks == 0 {
            return (0, 0.0);
        }
        let used_blocks = self.total_blocks - self.free_blocks;
        let used_bytes = used_blocks * self.block_size;
        let percent = (used_blocks as f64 / self.total_blocks as f64) * 100.0;
        (used_bytes, percent)
    }

    /// Calculate inode usage percentage
    pub fn inode_usage(&self) -> f64 {
        if self.total_inodes == 0 {
            return 0.0;
        }
        let used_inodes = self.total_inodes - self.free_inodes;
        (used_inodes as f64 / self.total_inodes as f64) * 100.0
    }
}

/// Get filesystem statistics using tune2fs
pub fn get_stats<P: AsRef<Path>>(path: P) -> Result<FilesystemStats> {
    check_dependencies(&["tune2fs"])?;

    let output = Command::new("tune2fs")
        .arg("-l")
        .arg(path.as_ref())
        .output()
        .context("Failed to execute tune2fs")?;

    if !output.status.success() {
        anyhow::bail!("tune2fs failed: {}", String::from_utf8_lossy(&output.stderr));
    }

    let stdout = String::from_utf8_lossy(&output.stdout);
    let mut stats = FilesystemStats::default();

    for line in stdout.lines() {
        if let Some((key, value)) = line.split_once(':') {
            let key = key.trim();
            let value = value.trim();

            match key {
                "Block size" => stats.block_size = value.parse().unwrap_or(0),
                "Block count" => stats.total_blocks = value.parse().unwrap_or(0),
                "Free blocks" => stats.free_blocks = value.parse().unwrap_or(0),
                "Inode count" => stats.total_inodes = value.parse().unwrap_or(0),
                "Free inodes" => stats.free_inodes = value.parse().unwrap_or(0),
                "Filesystem state" => stats.filesystem_state = value.to_string(),
                "Last mounted on" => stats.last_mounted = value.to_string(),
                _ => {}
            }
        }
    }

    Ok(stats)
}

/// Check filesystem integrity using e2fsck
pub fn check_integrity<P: AsRef<Path>>(path: P, force: bool) -> Result<()> {
    check_dependencies(&["e2fsck"])?;

    let path_str = path.as_ref().to_string_lossy().to_string();
    let mut args = vec!["-p"];
    if force {
        args.push("-f");
    }
    args.push(&path_str);

    utils::print_message(&format!(
        "Checking integrity of {} (force={})",
        utils::style_path(&path.as_ref().to_string_lossy()),
        force
    ));

    utils::print_debug(&format!("[CHECK] e2fsck {}", args.join(" ")));

    let output = Command::new("e2fsck")
        .args(&args)
        .output()
        .context("Failed to execute e2fsck")?;

    // e2fsck uses non-standard exit codes:
    // 0 = No errors
    // 1 = File system errors corrected (Success for our purposes)
    // 2+ = Critical errors or reboot required
    if let Some(code) = output.status.code() {
        if code > 1 {
            let stderr = String::from_utf8_lossy(&output.stderr);
            let stdout = String::from_utf8_lossy(&output.stdout);
            anyhow::bail!(
                "e2fsck failed with exit code {}: {}\n{}",
                code,
                stderr,
                stdout
            );
        }
    }

    if !output.stdout.is_empty() || !output.stderr.is_empty() {
        utils::print_debug(&format!(
            "e2fsck output: {}",
            String::from_utf8_lossy(&output.stdout)
        ));
    }

    Ok(())
}

/// Resize an overlay image to a new size
pub fn resize<P: AsRef<Path>>(path: P, size_mb: u32) -> Result<()> {
    check_dependencies(&["resize2fs", "e2fsck"])?;

    let path_ref = path.as_ref();
    if !path_ref.exists() {
        anyhow::bail!("Overlay image not found: {}", path_ref.display());
    }

    // Get current size
    let metadata = std::fs::metadata(path_ref)?;
    let current_size_bytes = metadata.len();
    let new_size_bytes = (size_mb as u64) * 1024 * 1024;

    if current_size_bytes == new_size_bytes {
        utils::print_message(&format!(
            "Size unchanged ({} MiB) for {}",
            size_mb,
            utils::style_path(&path_ref.to_string_lossy())
        ));
        return Ok(());
    }

    utils::print_message(&format!(
        "Resizing {} to {} MiB",
        utils::style_path(&path_ref.to_string_lossy()),
        size_mb
    ));

    // Pre-check integrity
    utils::print_message("Pre-check filesystem integrity...");
    check_integrity(path_ref, true)?;

    // Perform resize
    if new_size_bytes < current_size_bytes {
        // SHRINK
        utils::print_message(&format!("Shrinking overlay image to {} MiB", size_mb));
        
        let size_arg = format!("{}M", size_mb);
        run_command(
            "shrink filesystem",
            path_ref,
            "resize2fs",
            &["-p", &path_ref.to_string_lossy(), &size_arg],
        )?;

        std::fs::File::options()
            .write(true)
            .open(path_ref)?
            .set_len(new_size_bytes)
            .context("Failed to truncate file")?;

    } else {
        // EXPAND
        utils::print_message(&format!("Expanding overlay image to {} MiB", size_mb));
        
        std::fs::File::options()
            .write(true)
            .open(path_ref)?
            .set_len(new_size_bytes)
            .context("Failed to expand file")?;

        run_command(
            "expand filesystem",
            path_ref,
            "resize2fs",
            &["-p", &path_ref.to_string_lossy()],
        )?;
    }

    // Final verification
    utils::print_message("Final filesystem integrity check...");
    check_integrity(path_ref, true)?;

    utils::print_success(&format!(
        "Overlay image resized successfully to {} MiB: {}",
        size_mb,
        utils::style_path(&path_ref.to_string_lossy())
    ));

    Ok(())
}

/// Change ownership of files inside an overlay image
pub fn change_ownership<P: AsRef<Path>>(
    path: P,
    uid: u32,
    gid: u32,
    internal_path: &str,
) -> Result<()> {
    check_dependencies(&["debugfs"])?;

    let path_ref = path.as_ref();
    if !path_ref.exists() {
        anyhow::bail!("Overlay image not found: {}", path_ref.display());
    }

    // Map logical path to internal /upper structure
    let target_path = if internal_path == "/" {
        "/upper".to_string()
    } else {
        format!("/upper/{}", internal_path.trim_start_matches('/'))
    };

    utils::print_message(&format!(
        "Scanning {} inside {} (uid={} gid={})",
        utils::style_path(&target_path),
        utils::style_path(&path_ref.to_string_lossy()),
        uid,
        gid
    ));

    // Scan for inodes to update
    let inodes = scan_inodes(path_ref, &target_path)?;

    if inodes.is_empty() {
        utils::print_note(&format!(
            "No inodes found to modify at {}",
            utils::style_path(&target_path)
        ));
        return Ok(());
    }

    utils::print_message(&format!(
        "Updating {} inodes (uid={} gid={}) in {}",
        inodes.len(),
        uid,
        gid,
        utils::style_path(&target_path)
    ));

    // Build debugfs commands
    let mut commands = Vec::new();

    // Ensure skeleton paths are accessible if chowning root
    if internal_path == "/" || internal_path.is_empty() {
        for skeleton_path in &["/upper", "/work", "/work/work"] {
            commands.push(format!("sif {} uid {}", skeleton_path, uid));
            commands.push(format!("sif {} gid {}", skeleton_path, gid));
        }
    }

    // Fix discovered inodes
    for inode in &inodes {
        commands.push(format!("sif <{}> uid {}", inode, uid));
        commands.push(format!("sif <{}> gid {}", inode, gid));
    }
    commands.push("quit".to_string());

    // Execute batch commands
    let script = commands.join("\n");
    let mut child = Command::new("debugfs")
        .arg("-w")
        .arg(path_ref)
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .context("Failed to spawn debugfs")?;

    if let Some(mut stdin) = child.stdin.take() {
        use std::io::Write;
        stdin.write_all(script.as_bytes())
            .context("Failed to write to debugfs stdin")?;
    }

    let output = child.wait_with_output()
        .context("Failed to wait for debugfs")?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        anyhow::bail!("debugfs failed: {}", stderr);
    }

    utils::print_success(&format!(
        "Permissions updated for {}",
        utils::style_path(&path_ref.to_string_lossy())
    ));

    Ok(())
}

// =============================================================================
// Internal helper functions
// =============================================================================

/// Check that required system tools are available
fn check_dependencies(tools: &[&str]) -> Result<()> {
    let mut missing = Vec::new();

    for tool in tools {
        if which::which(tool).is_err() {
            missing.push(*tool);
        }
    }

    if !missing.is_empty() {
        anyhow::bail!(
            "Missing required system tools: {}",
            missing.join(", ")
        );
    }

    Ok(())
}

/// Run a command and wrap failures with context
fn run_command(op: &str, path: &Path, tool: &str, args: &[&str]) -> Result<()> {
    utils::print_debug(&format!("[{}] Running {} {}", op, tool, args.join(" ")));

    let output = Command::new(tool)
        .args(args)
        .output()
        .context(format!("Failed to execute {}", tool))?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        let stdout = String::from_utf8_lossy(&output.stdout);
        anyhow::bail!(
            "{} failed for {}: {}\n{}",
            tool,
            path.display(),
            stderr,
            stdout
        );
    }

    Ok(())
}

/// Scan inodes in an overlay image using debugfs
fn scan_inodes(img_path: &Path, start_path: &str) -> Result<Vec<String>> {
    use std::collections::{HashSet, VecDeque};

    // Get start inode via 'stat'
    let output = Command::new("debugfs")
        .arg("-R")
        .arg(format!("stat {}", start_path))
        .arg(img_path)
        .output()
        .context("Failed to execute debugfs stat")?;

    if !output.status.success() {
        anyhow::bail!(
            "Path {} not found in overlay layer (it may only exist in base image)",
            start_path
        );
    }

    // Parse inode number from stat output
    let stdout = String::from_utf8_lossy(&output.stdout);
    let inode_re = regex::Regex::new(r"Inode:\s+(\d+)").unwrap();
    let start_inode = inode_re
        .captures(&stdout)
        .and_then(|cap| cap.get(1))
        .ok_or_else(|| anyhow::anyhow!("Could not parse inode for {}", start_path))?
        .as_str()
        .to_string();

    // BFS walk the tree
    let mut visited: HashSet<String> = HashSet::new();
    let mut to_visit: VecDeque<String> = VecDeque::new();
    let mut results = vec![start_inode.clone()];

    to_visit.push_back(start_inode);

    let mut count = 0;
    let mut last_print = std::time::Instant::now();

    while let Some(current) = to_visit.pop_front() {
        if visited.contains(&current) {
            continue;
        }
        visited.insert(current.clone());

        // List contents: ls -l <INODE>
        let output = Command::new("debugfs")
            .arg("-R")
            .arg(format!("ls -l <{}>", current))
            .arg(img_path)
            .output()
            .context("Failed to execute debugfs ls")?;

        if !output.status.success() {
            continue;
        }

        let stdout = String::from_utf8_lossy(&output.stdout);
        for line in stdout.lines() {
            let parts: Vec<&str> = line.split_whitespace().collect();

            // Expected debugfs output format: "  2   40755 (2)      0      0    ..."
            if parts.len() < 6 {
                continue;
            }

            let inode = parts[0];
            let mode = parts[1];
            let name = parts[parts.len() - 1];

            if name == "." || name == ".." {
                continue;
            }

            // Validate inode is a number
            if inode.parse::<u64>().is_err() {
                continue;
            }

            results.push(inode.to_string());

            // Check if directory (Mode starts with 4 or 04)
            if mode.starts_with('4') || mode.starts_with("04") {
                to_visit.push_back(inode.to_string());
            }

            count += 1;
            if count % 1000 == 0 && last_print.elapsed().as_millis() > 100 {
                use std::io::Write;
                print!("\r    Scanning inodes... {} found", count);
                std::io::stdout().flush().ok();
                last_print = std::time::Instant::now();
            }
        }
    }

    // Clear progress line
    use std::io::Write;
    print!("\r{}\r", " ".repeat(50));
    std::io::stdout().flush().ok();
    
    utils::print_message(&format!("Finished scanning {} inodes", results.len()));

    Ok(results)
}
