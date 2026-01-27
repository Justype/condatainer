use std::env;

use crate::config::{self, validate_binary, AppConfig};

pub mod slurm;

#[derive(Debug)]
pub struct SchedulerInfo {
    pub scheduler_type: String,
    pub binary: String,
    pub version: Option<String>,
    pub in_job: bool,
    pub available: bool,
}

pub fn detect_scheduler(config: &AppConfig) -> Option<SchedulerInfo> {
    // Check if we have a configured scheduler binary
    if let Some(ref scheduler_bin) = config.scheduler_bin {
        if validate_binary(&scheduler_bin.to_string_lossy()) {
            return detect_scheduler_from_binary(scheduler_bin.to_string_lossy().as_ref());
        }
    }

    // Try to auto-detect
    if let Some((bin, sched_type)) = config::detect_scheduler_bin() {
        return detect_scheduler_from_binary(&bin.to_string_lossy());
    }

    None
}

fn detect_scheduler_from_binary(binary_path: &str) -> Option<SchedulerInfo> {
    // Determine scheduler type from binary name
    let binary_name = std::path::Path::new(binary_path)
        .file_name()?
        .to_str()?;

    let scheduler_type = match binary_name {
        "sbatch" | "srun" | "scancel" | "squeue" => "SLURM",
        "qsub" | "qstat" | "qdel" => {
            // Could be PBS or SGE, check environment
            if env::var("SGE_ROOT").is_ok() {
                "SGE"
            } else {
                "PBS"
            }
        }
        "bsub" | "bjobs" | "bkill" => "LSF",
        _ => return None,
    };

    // Delegate to scheduler-specific implementation
    match scheduler_type {
        "SLURM" => slurm::detect(binary_path),
        _ => {
            // For other schedulers, use basic detection
            Some(SchedulerInfo {
                scheduler_type: scheduler_type.to_string(),
                binary: binary_path.to_string(),
                version: None,
                in_job: is_inside_job(scheduler_type),
                available: !is_inside_job(scheduler_type),
            })
        }
    }
}

pub fn is_inside_job(scheduler_type: &str) -> bool {
    
    match scheduler_type {
        "SLURM" => slurm::is_inside_job(),
        "PBS" => {
            // PBS sets PBS_JOBID when inside a job
            env::var("PBS_JOBID").is_ok()
        }
        "SGE" => {
            // SGE sets JOB_ID when inside a job
            env::var("JOB_ID").is_ok()
        }
        "LSF" => {
            // LSF sets LSB_JOBID when inside a job
            env::var("LSB_JOBID").is_ok()
        }
        _ => false,
    }
}
