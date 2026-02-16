package cmd

import (
	"fmt"
	"os"
	"sort"
	"strings"
	"text/tabwriter"
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

var schedulerShowPartitions bool
var schedulerShowQueues bool
var schedulerShowCpuOnly bool
var schedulerShowGpuOnly bool

var schedulerCmd = &cobra.Command{
	Use:   "scheduler",
	Short: "Display scheduler information",
	Long: `Display information about the detected job scheduler.

Shows scheduler type (SLURM, PBS, etc.), binary path, version, and availability status.
Use -p to show per-partition resource limits.
Use --cpu or --gpu to filter by node type.`,
	Example: `  condatainer scheduler           # Show scheduler information
  condatainer scheduler -p        # Show per-partition limits (or use -q for queues)
  condatainer scheduler --gpu     # Show only GPU partitions
  condatainer scheduler -p --cpu  # Show per-partition limits for CPU-only partitions`,
	Run: runScheduler,
}

func init() {
	rootCmd.AddCommand(schedulerCmd)
	schedulerCmd.Flags().BoolVarP(&schedulerShowPartitions, "partitions", "p", false, "Show per-partition resource limits")
	schedulerCmd.Flags().BoolVarP(&schedulerShowQueues, "queue", "q", false, "Show per-queue resource limits (alias for -p)")
	schedulerCmd.Flags().BoolVar(&schedulerShowCpuOnly, "cpu", false, "Show only CPU-only partitions (no GPUs)")
	schedulerCmd.Flags().BoolVar(&schedulerShowGpuOnly, "gpu", false, "Show only GPU partitions")
}

func runScheduler(cmd *cobra.Command, args []string) {
	// Auto-enable -p if --cpu or --gpu is set, or if -q is used
	if schedulerShowCpuOnly || schedulerShowGpuOnly || schedulerShowQueues {
		schedulerShowPartitions = true
	}

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
		// Build GPU map by partition
		gpusByPartition := make(map[string][]scheduler.GpuInfo)
		for _, gpu := range clusterInfo.AvailableGpus {
			partition := gpu.Partition
			if partition == "" {
				partition = "default"
			}
			gpusByPartition[partition] = append(gpusByPartition[partition], gpu)
		}

		// Display per-partition limits if requested
		if schedulerShowPartitions {
			// Filter limits based on --cpu or --gpu flags (only for per-partition display)
			filteredLimits := make([]scheduler.ResourceLimits, 0)
			for _, limit := range clusterInfo.Limits {
				partitionName := limit.Partition
				if partitionName == "" {
					partitionName = "default"
				}

				hasGpu := len(gpusByPartition[partitionName]) > 0

				// Apply filters only if --cpu or --gpu is set
				if schedulerShowCpuOnly && hasGpu {
					continue // Skip GPU partitions when --cpu is set
				}
				if schedulerShowGpuOnly && !hasGpu {
					continue // Skip CPU-only partitions when --gpu is set
				}

				filteredLimits = append(filteredLimits, limit)
			}

			if len(filteredLimits) > 0 {
				fmt.Println()
				if schedulerShowCpuOnly {
					fmt.Println("Partition Resource Limits (CPU-only):")
				} else if schedulerShowGpuOnly {
					fmt.Println("Partition Resource Limits (GPU):")
				} else {
					fmt.Println("Partition Resource Limits:")
				}
				fmt.Println()

				// Create table writer
				w := tabwriter.NewWriter(os.Stdout, 0, 0, 2, ' ', 0)
				fmt.Fprintln(w, "PARTITION\tCPUs\tMEMORY\tTIME\tNODES\tGPU TYPES")
				fmt.Fprintln(w, "---------\t----\t------\t----\t-----\t---------")
				for _, limit := range filteredLimits {
					partitionName := limit.Partition
					if partitionName == "" {
						partitionName = "default"
					}

					// Format values
					cpuStr := "-"
					if limit.MaxCpus > 0 {
						cpuStr = fmt.Sprintf("%d", limit.MaxCpus)
					}

					memStr := "-"
					if limit.MaxMemMB > 0 {
						if limit.MaxMemMB >= 1024 {
							memStr = fmt.Sprintf("%.0f GB", float64(limit.MaxMemMB)/1024)
						} else {
							memStr = fmt.Sprintf("%d MB", limit.MaxMemMB)
						}
					}

					timeStr := "-"
					if limit.MaxTime > 0 {
						timeStr = formatDuration(limit.MaxTime)
					}

					nodesStr := "-"
					if limit.MaxNodes > 0 {
						nodesStr = fmt.Sprintf("%d", limit.MaxNodes)
					}

					// Collect GPU types for this partition
					gpuTypesStr := "-"
					if gpus, ok := gpusByPartition[partitionName]; ok && len(gpus) > 0 {
						// Merge same GPU types
						gpuMap := make(map[string]*scheduler.GpuInfo)
						for _, gpu := range gpus {
							if existing, ok := gpuMap[gpu.Type]; ok {
								existing.Total += gpu.Total
								existing.Available += gpu.Available
							} else {
								gpuCopy := gpu
								gpuMap[gpu.Type] = &gpuCopy
							}
						}

						// Sort GPU names
						gpuNames := make([]string, 0, len(gpuMap))
						for name := range gpuMap {
							gpuNames = append(gpuNames, name)
						}
						sort.Strings(gpuNames)

						// Format GPU types
						gpuLines := make([]string, 0)
						for _, gpuType := range gpuNames {
							gpu := gpuMap[gpuType]
							gpuLines = append(gpuLines, fmt.Sprintf("%s (%d/%d)", gpu.Type, gpu.Available, gpu.Total))
						}
						gpuTypesStr = strings.Join(gpuLines, "\n")
					}

					// Print row (handle multi-line GPU types)
					if strings.Contains(gpuTypesStr, "\n") {
						// Multiple GPU types - print first line with all data
						gpuTypeLines := strings.Split(gpuTypesStr, "\n")
						fmt.Fprintf(w, "%s\t%s\t%s\t%s\t%s\t%s\n",
							partitionName, cpuStr, memStr, timeStr, nodesStr, gpuTypeLines[0])
						// Print remaining GPU types on separate rows
						for i := 1; i < len(gpuTypeLines); i++ {
							fmt.Fprintf(w, "\t\t\t\t\t%s\n", gpuTypeLines[i])
						}
					} else {
						// Single or no GPU type - print normally
						fmt.Fprintf(w, "%s\t%s\t%s\t%s\t%s\t%s\n",
							partitionName, cpuStr, memStr, timeStr, nodesStr, gpuTypesStr)
					}
				}

				w.Flush()

				return
			}

			// No partitions match the filter
			fmt.Println()
			if schedulerShowGpuOnly {
				utils.PrintWarning("No GPU partitions available")
			} else if schedulerShowCpuOnly {
				utils.PrintWarning("No CPU-only partitions available")
			}
			return
		}

		// Display max resource limits (aggregated from all partitions)
		maxCpus := 0
		var maxMemMB int64 = 0
		var maxTime time.Duration = 0

		// Find max values across all partitions (not filtered)
		for _, limit := range clusterInfo.Limits {
			if limit.MaxCpus > maxCpus {
				maxCpus = limit.MaxCpus
			}
			if limit.MaxMemMB > maxMemMB {
				maxMemMB = limit.MaxMemMB
			}
			if limit.MaxTime > maxTime {
				maxTime = limit.MaxTime
			}
		}

		if maxCpus > 0 || maxMemMB > 0 || maxTime > 0 {
			fmt.Println()
			fmt.Println("Max Resource Limits (across all partitions):")
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
			if maxTime > 0 {
				fmt.Printf("  Max Time:   %s\n", utils.StyleNumber(formatDuration(maxTime)))
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
