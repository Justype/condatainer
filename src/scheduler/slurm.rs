use std::env;
use std::process::Command;

use super::SchedulerInfo;

/// Detect SLURM scheduler information
pub fn detect(binary_path: &str) -> Option<SchedulerInfo> {
    let version = get_version(binary_path);
    let in_job = is_inside_job();
    let available = !in_job;

    Some(SchedulerInfo {
        scheduler_type: "SLURM".to_string(),
        binary: binary_path.to_string(),
        version,
        in_job,
        available,
    })
}

/// Check if we're inside a SLURM job
pub fn is_inside_job() -> bool {
    // SLURM sets SLURM_JOB_ID when inside a job
    env::var("SLURM_JOB_ID").is_ok() || env::var("SLURM_JOBID").is_ok()
}

/// Get SLURM version
fn get_version(binary_path: &str) -> Option<String> {
    // Try sbatch --version
    let output = if binary_path.contains("sbatch") {
        Command::new(binary_path).arg("--version").output().ok()?
    } else {
        // Try to find sbatch
        let sbatch = which::which("sbatch").ok()?;
        Command::new(sbatch).arg("--version").output().ok()?
    };

    if !output.status.success() {
        return None;
    }

    let stdout = String::from_utf8(output.stdout).ok()?;
    
    // Output like "slurm 23.02.7"
    for line in stdout.lines() {
        if let Some(version_str) = line.split_whitespace().nth(1) {
            return Some(version_str.to_string());
        }
    }

    None
}
