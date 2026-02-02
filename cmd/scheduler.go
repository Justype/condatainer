package cmd

import (
	"fmt"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

var schedulerCmd = &cobra.Command{
	Use:   "scheduler",
	Short: "Display scheduler information",
	Long: `Display information about the detected job scheduler.

Shows scheduler type (SLURM, PBS, etc.), binary path, version, and availability status.`,
	Example: `  condatainer scheduler           # Show scheduler information
  condatainer sched              # Short alias`,
	Run: runScheduler,
}

func init() {
	rootCmd.AddCommand(schedulerCmd)
}

func runScheduler(cmd *cobra.Command, args []string) {
	// Try to detect scheduler
	sched, err := scheduler.DetectSchedulerWithBinary(config.Global.SchedulerBin)

	if err != nil {
		// If we're inside a scheduled job, show a concise message and exit
		if scheduler.IsInsideJob() {
			utils.PrintMessage("Scheduler Status: %s", utils.StyleWarning("Unavailable (inside job)"))
			utils.PrintMessage("")
			utils.PrintMessage("You are currently inside a scheduled job; job submission is disabled to prevent nested submissions.")
			return
		}

		// No scheduler found
		utils.PrintMessage("Scheduler Status: %s", utils.StyleError("Not Found"))
		utils.PrintMessage("")
		utils.PrintMessage("No job scheduler detected on this system.")
		utils.PrintMessage("Supported schedulers: SLURM (more coming soon)")
		return
	}

	// Get scheduler info
	info := sched.GetInfo()

	// Display scheduler information (no [CNT] prefix for structured output)
	fmt.Println("Scheduler Information:")
	fmt.Printf("  Type:      %s\n", utils.StyleInfo(info.Type))
	fmt.Printf("  Binary:    %s\n", utils.StylePath(info.Binary))

	if info.Version != "" {
		fmt.Printf("  Version:   %s\n", utils.StyleNumber(info.Version))
	}

	if info.InJob {
		fmt.Printf("  Status:    %s (inside job)\n", utils.StyleError("Unavailable"))
		fmt.Println()
		fmt.Println("You are currently inside a scheduled job (detected via environment).")
		fmt.Println("Job submission is disabled to prevent nested job submissions.")
		// Do not query or print cluster information when inside a job
		return
	} else if info.Available {
		fmt.Printf("  Status:    %s\n", utils.StyleSuccess("Available"))
		fmt.Println()
		fmt.Println("The scheduler is available and ready for job submission.")
	} else {
		fmt.Printf("  Status:    %s\n", utils.StyleError("Unavailable"))
		fmt.Println()
		fmt.Println("Scheduler detected but not available for job submission.")
	}

	// Try to get cluster info
	clusterInfo, err := sched.GetClusterInfo()
	if err == nil && clusterInfo != nil {
		// Display GPU information if available
		if len(clusterInfo.AvailableGpus) > 0 {
			fmt.Println()
			fmt.Println("Available GPUs:")
			for _, gpu := range clusterInfo.AvailableGpus {
				if gpu.Partition != "" {
					fmt.Printf("  %s: %d total, %d available (partition: %s)\n",
						utils.StyleName(gpu.Type), gpu.Total, gpu.Available, gpu.Partition)
				} else {
					fmt.Printf("  %s: %d total, %d available\n",
						utils.StyleName(gpu.Type), gpu.Total, gpu.Available)
				}
			}
		}

		// Display resource limits if available
		if len(clusterInfo.Limits) > 0 {
			fmt.Println()
			fmt.Println("Resource Limits:")
			for _, limit := range clusterInfo.Limits {
				partition := limit.Partition
				if partition == "" {
					partition = "default"
				}
				fmt.Printf("  Partition: %s\n", utils.StyleName(partition))
				if limit.MaxCpus > 0 {
					fmt.Printf("    Max CPUs:   %s\n", utils.StyleNumber(fmt.Sprintf("%d", limit.MaxCpus)))
				}
				if limit.MaxMemMB > 0 {
					fmt.Printf("    Max Memory: %s\n", utils.StyleNumber(fmt.Sprintf("%d MB", limit.MaxMemMB)))
				}
				if limit.MaxGpus > 0 {
					fmt.Printf("    Max GPUs:   %s\n", utils.StyleNumber(fmt.Sprintf("%d", limit.MaxGpus)))
				}
				if limit.MaxTime > 0 {
					fmt.Printf("    Max Time:   %s\n", utils.StyleNumber(limit.MaxTime.String()))
				}
				if limit.DefaultTime > 0 {
					fmt.Printf("    Default Time: %s\n", utils.StyleNumber(limit.DefaultTime.String()))
				}
			}
		}
	}
}
