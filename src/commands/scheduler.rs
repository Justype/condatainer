use anyhow::Result;
use clap::Args;

use crate::config::AppConfig;
use crate::scheduler;
use crate::utils;

#[derive(Args, Debug)]
pub struct SchedulerArgs {}

pub fn handle_scheduler(_args: &SchedulerArgs, config: &AppConfig) -> Result<()> {
    // Try to detect scheduler
    let scheduler_info = scheduler::detect_scheduler(config);

    if let Some(info) = scheduler_info {
        display_scheduler_info(&info);
    } else {
        // No scheduler found
        utils::print_message(&format!(
            "Scheduler Status: {}",
            utils::style_error("Not Found")
        ));
        println!();
        utils::print_message("No job scheduler detected on this system.");
        utils::print_message("Supported schedulers: SLURM (more coming soon)");
    }

    Ok(())
}

fn display_scheduler_info(info: &scheduler::SchedulerInfo) {
    println!("Scheduler Information:");
    println!("  Type:      {}", utils::style_info(&info.scheduler_type));
    println!("  Binary:    {}", utils::style_path(&info.binary));

    if let Some(ref version) = info.version {
        println!("  Version:   {}", utils::style_number(version));
    }

    if info.in_job {
        println!(
            "  Status:    {} (inside job)",
            utils::style_error("Unavailable")
        );
        println!();
        println!("You are currently inside a scheduled job (detected via environment).");
        println!("Job submission is disabled to prevent nested job submissions.");
    } else if info.available {
        println!("  Status:    {}", utils::style_success("Available"));
        println!();
        println!("The scheduler is available and ready for job submission.");
    } else {
        println!("  Status:    {}", utils::style_error("Unavailable"));
        println!();
        println!("Scheduler detected but not available for job submission.");
    }

    // Note: Advanced features like GPU info and cluster limits could be added later
    // by parsing sinfo/scontrol output for SLURM, etc.
}
