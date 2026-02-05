package cmd

import (
	"fmt"
	"sort"
	"time"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

// formatDuration formats a duration with days for readability
func formatDuration(d time.Duration) string {
	days := int(d.Hours()) / 24
	hours := int(d.Hours()) % 24
	minutes := int(d.Minutes()) % 60

	if days > 0 {
		if hours > 0 {
			return fmt.Sprintf("%dd %dh", days, hours)
		}
		return fmt.Sprintf("%dd", days)
	}
	if hours > 0 {
		if minutes > 0 {
			return fmt.Sprintf("%dh %dm", hours, minutes)
		}
		return fmt.Sprintf("%dh", hours)
	}
	return fmt.Sprintf("%dm", minutes)
}

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
		// Display max resource limits
		maxLimits, _ := scheduler.ReadMaxSpecs(sched)

		// Determine max CPUs and memory (use partition limits if set, otherwise node info)
		maxCpus := 0
		var maxMemMB int64 = 0
		var maxTime = ""

		if maxLimits != nil {
			maxCpus = maxLimits.MaxCpus
			maxMemMB = maxLimits.MaxMemMB
			if maxLimits.MaxTime > 0 {
				maxTime = formatDuration(maxLimits.MaxTime)
			}
		}

		// Fall back to node resources if partition limits don't specify max
		if maxCpus == 0 && clusterInfo.MaxCpusPerNode > 0 {
			maxCpus = clusterInfo.MaxCpusPerNode
		}
		if maxMemMB == 0 && clusterInfo.MaxMemMBPerNode > 0 {
			maxMemMB = clusterInfo.MaxMemMBPerNode
		}

		if maxCpus > 0 || maxMemMB > 0 || maxTime != "" {
			fmt.Println()
			fmt.Println("Max Resource Limits:")
			if maxCpus > 0 {
				fmt.Printf("  Max CPUs:   %s\n", utils.StyleNumber(fmt.Sprintf("%d", maxCpus)))
			}
			if maxMemMB > 0 {
				var memStr string
				if maxMemMB >= 1024 {
					memStr = fmt.Sprintf("%d MB (%.0f GB)", maxMemMB, float64(maxMemMB)/1024)
				} else {
					memStr = fmt.Sprintf("%d MB", maxMemMB)
				}
				fmt.Printf("  Max Memory: %s\n", utils.StyleNumber(memStr))
			}
			if maxTime != "" {
				fmt.Printf("  Max Time:   %s\n", utils.StyleNumber(maxTime))
			}
		}

		// Display GPU information if available (merged by type)
		if len(clusterInfo.AvailableGpus) > 0 {
			// Merge GPUs with the same type name
			gpuMap := make(map[string]*scheduler.GpuInfo)
			for _, gpu := range clusterInfo.AvailableGpus {
				if existing, ok := gpuMap[gpu.Type]; ok {
					existing.Total += gpu.Total
					existing.Available += gpu.Available
				} else {
					gpuMap[gpu.Type] = &scheduler.GpuInfo{
						Type:      gpu.Type,
						Total:     gpu.Total,
						Available: gpu.Available,
					}
				}
			}

			// Sort GPU names
			gpuNames := make([]string, 0, len(gpuMap))
			for name := range gpuMap {
				gpuNames = append(gpuNames, name)
			}
			sort.Strings(gpuNames)

			fmt.Println()
			fmt.Println("Available GPUs:")
			for _, name := range gpuNames {
				gpu := gpuMap[name]
				fmt.Printf("  %s: %d total, %d available\n",
					utils.StyleName(gpu.Type), gpu.Total, gpu.Available)
			}
		}
	}
}
