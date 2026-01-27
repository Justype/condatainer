use colored::*;
use std::fs;
use std::io::{self, IsTerminal};
use std::path::Path;
use std::sync::atomic::{AtomicBool, Ordering};

/// Debug mode flag - controls whether debug output is visible
static DEBUG_MODE: AtomicBool = AtomicBool::new(false);

/// Project prefix for all log messages
const PROJECT_PREFIX: &str = "[CNT]";

/// Set debug mode
pub fn set_debug_mode(enabled: bool) {
    DEBUG_MODE.store(enabled, Ordering::Relaxed);
}

/// Check if debug mode is enabled
pub fn is_debug_mode() -> bool {
    DEBUG_MODE.load(Ordering::Relaxed)
}

// ---------------------------------------------------------
// Semantic Styles - Simple one-line color functions
// ---------------------------------------------------------

pub fn style_error(msg: &str) -> String { msg.red().to_string() }
pub fn style_success(msg: &str) -> String { msg.green().to_string() }
pub fn style_warning(msg: &str) -> String { msg.yellow().to_string() }
pub fn style_hint(msg: &str) -> String { msg.cyan().to_string() }
pub fn style_note(msg: &str) -> String { msg.magenta().to_string() }
pub fn style_info(msg: &str) -> String { msg.magenta().to_string() }
pub fn style_debug(msg: &str) -> String { msg.bright_black().to_string() }
pub fn style_command(cmd: &str) -> String { cmd.bright_black().to_string() }
pub fn style_action(act: &str) -> String { act.yellow().to_string() }
pub fn style_title(title: &str) -> String { title.cyan().bold().to_string() }
pub fn style_name(name: &str) -> String { name.yellow().to_string() }
pub fn style_highlight(text: &str) -> String { text.yellow().bold().to_string() }

pub fn style_number<T: std::fmt::Display>(num: T) -> String { 
    format!("{}", num).magenta().to_string() 
}

pub fn style_path(path: &str) -> String {
    if is_img(path) { path.magenta().bold().to_string() }
    else if is_sif(path) || is_sqf(path) { path.cyan().bold().to_string() }
    else { path.blue().bold().to_string() }
}

// ---------------------------------------------------------
// Log Printers - Simplified implementations
// ---------------------------------------------------------

pub fn print_message(msg: &str) {
    println!("{} {}", PROJECT_PREFIX, msg);
}

pub fn print_success(msg: &str) {
    println!("{}{} {}", PROJECT_PREFIX, style_success("[PASS]"), msg);
}

pub fn print_error(msg: &str) {
    eprintln!("{}{} {}", PROJECT_PREFIX, style_error("[ERR] "), msg);
}

pub fn print_warning(msg: &str) {
    eprintln!("{}{} {}", PROJECT_PREFIX, style_warning("[WARN]"), msg);
}

pub fn print_hint(msg: &str) {
    println!("{}{} {}", PROJECT_PREFIX, style_hint("[HINT]"), msg);
}

pub fn print_note(msg: &str) {
    println!("{}{} {}", PROJECT_PREFIX, style_note("[NOTE]"), msg);
}

pub fn print_debug(msg: &str) {
    if is_debug_mode() {
        eprintln!("{}{} {}", PROJECT_PREFIX, style_debug("[DBG] "), msg);
    }
}

// Formatting macros for convenience
#[macro_export]
macro_rules! print_message {
    ($($arg:tt)*) => { $crate::utils::print_message(&format!($($arg)*)); };
}

#[macro_export]
macro_rules! print_success {
    ($($arg:tt)*) => { $crate::utils::print_success(&format!($($arg)*)); };
}

#[macro_export]
macro_rules! print_error {
    ($($arg:tt)*) => { $crate::utils::print_error(&format!($($arg)*)); };
}

#[macro_export]
macro_rules! print_warning {
    ($($arg:tt)*) => { $crate::utils::print_warning(&format!($($arg)*)); };
}

#[macro_export]
macro_rules! print_hint {
    ($($arg:tt)*) => { $crate::utils::print_hint(&format!($($arg)*)); };
}

#[macro_export]
macro_rules! print_note {
    ($($arg:tt)*) => { $crate::utils::print_note(&format!($($arg)*)); };
}

#[macro_export]
macro_rules! print_debug {
    ($($arg:tt)*) => { $crate::utils::print_debug(&format!($($arg)*)); };
}

// ---------------------------------------------------------
// Terminal Detection
// ---------------------------------------------------------

/// Check if stdout is connected to a TTY (interactive terminal)
pub fn is_interactive_shell() -> bool {
    io::stdout().is_terminal()
}

// ---------------------------------------------------------
// File Extension Checks
// ---------------------------------------------------------

/// Check if the path has an ext3 overlay extension (.img)
/// Note: In Apptainer context, .img usually implies a writable ext3 overlay
pub fn is_img(path: &str) -> bool {
    Path::new(path)
        .extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| ext.eq_ignore_ascii_case("img"))
        .unwrap_or(false)
}

/// Check if the path has a SquashFS extension (.sqf, .sqsh, .squashfs)
/// These are read-only compressed images
pub fn is_sqf(path: &str) -> bool {
    Path::new(path)
        .extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| {
            ext.eq_ignore_ascii_case("sqf")
                || ext.eq_ignore_ascii_case("sqsh")
                || ext.eq_ignore_ascii_case("squashfs")
        })
        .unwrap_or(false)
}

/// Check if the path has a Singularity Image Format extension (.sif)
/// This is the native format for Apptainer/Singularity
pub fn is_sif(path: &str) -> bool {
    Path::new(path)
        .extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| ext.eq_ignore_ascii_case("sif"))
        .unwrap_or(false)
}

/// Check if a filename is an overlay file (.sqf, .sqsh, .squashfs, .img, .sif)
pub fn is_overlay(filename: &str) -> bool {
    filename.ends_with(".sqf")
        || filename.ends_with(".sqsh")
        || filename.ends_with(".squashfs")
        || filename.ends_with(".img")
        || filename.ends_with(".sif")
}

// ---------------------------------------------------------
// Permission Management
// ---------------------------------------------------------

/// Standard default file permissions: u=rw, g=rw, o=r
#[cfg(unix)]
pub const PERM_FILE: u32 = 0o664;

/// Standard default directory permissions: u=rwx, g=rwx, o=rx
#[cfg(unix)]
pub const PERM_DIR: u32 = 0o775;

/// Fix permissions automatically using standard defaults (0664/0775)
#[cfg(unix)]
pub fn fix_permissions_default<P: AsRef<Path>>(path: P) -> io::Result<()> {
    fix_permissions(path, PERM_FILE, PERM_DIR)
}

/// Fix permissions automatically for a path
/// If path is a file, sets it to filePerm
/// If path is a directory, recursively sets files to filePerm and subdirs to dirPerm
#[cfg(unix)]
pub fn fix_permissions<P: AsRef<Path>>(path: P, file_perm: u32, dir_perm: u32) -> io::Result<()> {
    use std::os::unix::fs::PermissionsExt;

    let path = path.as_ref();
    let metadata = fs::metadata(path)?;

    // Handle single file
    if metadata.is_file() {
        print_debug(&format!(
            "Fixing permissions for file: {} [0o{:o}]",
            style_path(&path.display().to_string()),
            file_perm
        ));
        let permissions = fs::Permissions::from_mode(file_perm);
        fs::set_permissions(path, permissions)?;
        return Ok(());
    }

    // Handle directory (recursive walk)
    print_debug(&format!(
        "Recursively fixing permissions for directory: {} [Dirs:0o{:o}, Files:0o{:o}]",
        style_path(&path.display().to_string()),
        dir_perm,
        file_perm
    ));

    // First, fix the root directory itself
    let permissions = fs::Permissions::from_mode(dir_perm);
    fs::set_permissions(path, permissions)?;

    // Then walk the contents
    walk_and_fix_permissions(path, file_perm, dir_perm)?;

    Ok(())
}

#[cfg(unix)]
fn walk_and_fix_permissions<P: AsRef<Path>>(
    path: P,
    file_perm: u32,
    dir_perm: u32,
) -> io::Result<()> {
    use std::os::unix::fs::PermissionsExt;

    for entry in fs::read_dir(path.as_ref())? {
        let entry = match entry {
            Ok(e) => e,
            Err(err) => {
                print_warning(&format!(
                    "Skipping inaccessible path during permission fix: {}",
                    err
                ));
                continue;
            }
        };

        let path = entry.path();
        let metadata = match entry.metadata() {
            Ok(m) => m,
            Err(err) => {
                print_warning(&format!(
                    "Could not get metadata for {}: {}",
                    style_path(&path.display().to_string()),
                    err
                ));
                continue;
            }
        };

        let target_mode = if metadata.is_dir() {
            dir_perm
        } else {
            file_perm
        };

        // Only run chmod if permissions actually differ
        if metadata.permissions().mode() & 0o777 != target_mode {
            let permissions = fs::Permissions::from_mode(target_mode);
            if let Err(err) = fs::set_permissions(&path, permissions) {
                print_warning(&format!(
                    "Could not chmod {}: {}",
                    style_path(&path.display().to_string()),
                    err
                ));
            }
        }

        // Recurse into directories
        if metadata.is_dir() {
            walk_and_fix_permissions(&path, file_perm, dir_perm)?;
        }
    }

    Ok(())
}

// Non-Unix platforms - stub implementations
#[cfg(not(unix))]
pub fn fix_permissions_default<P: AsRef<Path>>(_path: P) -> io::Result<()> {
    Ok(())
}

#[cfg(not(unix))]
pub fn fix_permissions<P: AsRef<Path>>(_path: P, _file_perm: u32, _dir_perm: u32) -> io::Result<()> {
    Ok(())
}
