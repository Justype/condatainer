package scheduler

import (
	"fmt"
	"sort"
	"strings"
)

// GpuTier represents the performance tier of GPUs for compatibility matching
type GpuTier int

const (
	TierUnknown GpuTier = iota
	TierEntry           // Consumer/Entry-level GPUs
	TierMid             // Mid-range GPUs (P100, V100, T4, etc.)
	TierHigh            // High-end GPUs (A100, H100, etc.)
)

// GpuCompatibility defines GPU compatibility rules and performance tiers
type GpuCompatibility struct {
	Type        string   // GPU type (e.g., "h100", "a100")
	Tier        GpuTier  // Performance tier
	Compute     float64  // Compute capability (e.g., 8.0, 9.0)
	MemoryGB    int      // Typical memory in GB
	Aliases     []string // Alternative names for this GPU
	UpgradeTo   []string // Better GPUs that can substitute this one
	DowngradeTo []string // Fallback GPUs if this is unavailable
}

// gpuDatabase is a knowledge base of GPU types and their compatibility
var gpuDatabase = map[string]GpuCompatibility{
	// NVIDIA H-series (Hopper architecture)
	"h100": {
		Type:        "h100",
		Tier:        TierHigh,
		Compute:     9.0,
		MemoryGB:    80,
		Aliases:     []string{"h100", "hopper"},
		UpgradeTo:   []string{},
		DowngradeTo: []string{"a100", "a40"},
	},

	// Note: MIG profiles are dynamically discovered from SLURM via sinfo
	// and automatically added to the database at runtime by ensureMigProfileInDatabase()

	// NVIDIA A-series (Ampere architecture)
	"a100": {
		Type:        "a100",
		Tier:        TierHigh,
		Compute:     8.0,
		MemoryGB:    40,
		Aliases:     []string{"a100", "ampere"},
		UpgradeTo:   []string{"h100"},
		DowngradeTo: []string{"a40", "a30", "v100"},
	},
	"a40": {
		Type:        "a40",
		Tier:        TierHigh,
		Compute:     8.6,
		MemoryGB:    48,
		Aliases:     []string{"a40"},
		UpgradeTo:   []string{"a100", "h100"},
		DowngradeTo: []string{"a30", "v100"},
	},
	"a30": {
		Type:        "a30",
		Tier:        TierMid,
		Compute:     8.0,
		MemoryGB:    24,
		Aliases:     []string{"a30"},
		UpgradeTo:   []string{"a40", "a100"},
		DowngradeTo: []string{"v100", "t4"},
	},

	// NVIDIA V-series (Volta architecture)
	"v100": {
		Type:        "v100",
		Tier:        TierMid,
		Compute:     7.0,
		MemoryGB:    32,
		Aliases:     []string{"v100", "volta"},
		UpgradeTo:   []string{"a30", "a100"},
		DowngradeTo: []string{"p100", "t4"},
	},

	// NVIDIA P-series (Pascal architecture)
	"p100": {
		Type:        "p100",
		Tier:        TierMid,
		Compute:     6.0,
		MemoryGB:    16,
		Aliases:     []string{"p100", "pascal"},
		UpgradeTo:   []string{"v100", "a30"},
		DowngradeTo: []string{"p40", "k80"},
	},
	"p40": {
		Type:        "p40",
		Tier:        TierMid,
		Compute:     6.1,
		MemoryGB:    24,
		Aliases:     []string{"p40"},
		UpgradeTo:   []string{"v100", "a30"},
		DowngradeTo: []string{"p100"},
	},

	// NVIDIA T-series (Turing architecture)
	"t4": {
		Type:        "t4",
		Tier:        TierMid,
		Compute:     7.5,
		MemoryGB:    16,
		Aliases:     []string{"t4", "turing"},
		UpgradeTo:   []string{"a30", "v100"},
		DowngradeTo: []string{"p100"},
	},

	// NVIDIA K-series (Kepler architecture)
	"k80": {
		Type:        "k80",
		Tier:        TierEntry,
		Compute:     3.7,
		MemoryGB:    24,
		Aliases:     []string{"k80", "kepler"},
		UpgradeTo:   []string{"p100", "v100"},
		DowngradeTo: []string{},
	},

	// Generic fallback
	"gpu": {
		Type:        "gpu",
		Tier:        TierUnknown,
		Compute:     0.0,
		MemoryGB:    0,
		Aliases:     []string{"gpu", "any"},
		UpgradeTo:   []string{},
		DowngradeTo: []string{},
	},
}

// NormalizeGpuType normalizes GPU type strings to canonical form
func NormalizeGpuType(gpuType string) string {
	normalized := strings.ToLower(strings.TrimSpace(gpuType))

	// Check if this is a MIG profile (contains dots - e.g., "1g.10gb")
	// MIG profiles should be kept mostly intact
	if strings.Contains(normalized, ".") && (strings.Contains(normalized, "g.") || strings.Contains(normalized, "gb")) {
		// This looks like a MIG profile
		// First check if it's already in the database
		if _, found := gpuDatabase[normalized]; found {
			return normalized
		}

		// Check aliases for MIG profiles
		for canonical, compat := range gpuDatabase {
			for _, alias := range compat.Aliases {
				if normalized == alias {
					return canonical
				}
			}
		}

		// Return as-is if not found (it might be a valid MIG profile we don't know about)
		return normalized
	}

	// Remove common prefixes/suffixes (before removing dashes)
	normalized = strings.TrimPrefix(normalized, "nvidia-")
	normalized = strings.TrimPrefix(normalized, "nvidia_")
	normalized = strings.TrimPrefix(normalized, "tesla-")
	normalized = strings.TrimPrefix(normalized, "tesla_")
	normalized = strings.TrimSuffix(normalized, "-sxm")
	normalized = strings.TrimSuffix(normalized, "-pcie")

	// Handle common variations
	normalized = strings.ReplaceAll(normalized, "_", "")
	normalized = strings.ReplaceAll(normalized, "-", "")

	// Map aliases to canonical names
	for canonical, compat := range gpuDatabase {
		for _, alias := range compat.Aliases {
			if normalized == alias {
				return canonical
			}
		}
	}

	return normalized
}

// GpuConversionOption represents a possible GPU conversion/fallback
type GpuConversionOption struct {
	OriginalSpec  *GpuSpec // Original GPU request
	SuggestedSpec *GpuSpec // Suggested GPU to use
	Available     int      // Number of this GPU type available
	Total         int      // Total number of this GPU type
	Partition     string   // Partition where this GPU is available
	IsUpgrade     bool     // True if suggested GPU is better than original
	IsDowngrade   bool     // True if suggested GPU is worse than original
	Reason        string   // Explanation for this suggestion
}

// FindAvailableGpuTypes returns all GPU types available on the cluster
func FindAvailableGpuTypes(info *ClusterInfo) []GpuInfo {
	if info == nil {
		return []GpuInfo{}
	}

	available := make([]GpuInfo, 0)
	for _, gpu := range info.AvailableGpus {
		if gpu.Available > 0 {
			available = append(available, gpu)
		}
	}

	return available
}

// FindCompatibleGpu finds compatible GPU alternatives when requested type is unavailable
func FindCompatibleGpu(requestedSpec *GpuSpec, info *ClusterInfo) ([]*GpuConversionOption, error) {
	if requestedSpec == nil {
		return nil, nil
	}

	if info == nil {
		return nil, ErrClusterInfoUnavailable
	}

	options := make([]*GpuConversionOption, 0)
	requestedType := NormalizeGpuType(requestedSpec.Type)
	requestedCompat, knownType := gpuDatabase[requestedType]

	// First, check if the exact requested type is available
	for _, gpu := range info.AvailableGpus {
		normalizedClusterType := NormalizeGpuType(gpu.Type)
		if normalizedClusterType == requestedType && gpu.Available >= requestedSpec.Count {
			options = append(options, &GpuConversionOption{
				OriginalSpec: requestedSpec,
				SuggestedSpec: &GpuSpec{
					Type:  gpu.Type,
					Count: requestedSpec.Count,
					Raw:   fmt.Sprintf("%s:%d", gpu.Type, requestedSpec.Count),
				},
				Available:   gpu.Available,
				Total:       gpu.Total,
				Partition:   gpu.Partition,
				IsUpgrade:   false,
				IsDowngrade: false,
				Reason:      "Exact match available",
			})
		}
	}

	// If exact match found and available, return it first
	if len(options) > 0 {
		return options, nil
	}

	// If requested type is unknown, try to find any available GPU
	if !knownType {
		for _, gpu := range info.AvailableGpus {
			if gpu.Available >= requestedSpec.Count {
				options = append(options, &GpuConversionOption{
					OriginalSpec: requestedSpec,
					SuggestedSpec: &GpuSpec{
						Type:  gpu.Type,
						Count: requestedSpec.Count,
						Raw:   fmt.Sprintf("%s:%d", gpu.Type, requestedSpec.Count),
					},
					Available:   gpu.Available,
					Total:       gpu.Total,
					Partition:   gpu.Partition,
					IsUpgrade:   false,
					IsDowngrade: false,
					Reason:      fmt.Sprintf("Unknown GPU type '%s', suggesting available GPU", requestedSpec.Type),
				})
			}
		}
		return options, nil
	}

	// Try to find compatible upgrades (better GPUs)
	for _, upgradeType := range requestedCompat.UpgradeTo {
		for _, gpu := range info.AvailableGpus {
			normalizedClusterType := NormalizeGpuType(gpu.Type)
			if normalizedClusterType == upgradeType && gpu.Available >= requestedSpec.Count {
				upgradeCompat := gpuDatabase[upgradeType]
				options = append(options, &GpuConversionOption{
					OriginalSpec: requestedSpec,
					SuggestedSpec: &GpuSpec{
						Type:  gpu.Type,
						Count: requestedSpec.Count,
						Raw:   fmt.Sprintf("%s:%d", gpu.Type, requestedSpec.Count),
					},
					Available:   gpu.Available,
					Total:       gpu.Total,
					Partition:   gpu.Partition,
					IsUpgrade:   true,
					IsDowngrade: false,
					Reason: fmt.Sprintf("Upgrade from %s (compute %.1f) to %s (compute %.1f)",
						requestedType, requestedCompat.Compute, upgradeType, upgradeCompat.Compute),
				})
			}
		}
	}

	// Try to find compatible downgrades (fallback GPUs)
	for _, downgradeType := range requestedCompat.DowngradeTo {
		for _, gpu := range info.AvailableGpus {
			normalizedClusterType := NormalizeGpuType(gpu.Type)
			if normalizedClusterType == downgradeType && gpu.Available >= requestedSpec.Count {
				downgradeCompat := gpuDatabase[downgradeType]
				options = append(options, &GpuConversionOption{
					OriginalSpec: requestedSpec,
					SuggestedSpec: &GpuSpec{
						Type:  gpu.Type,
						Count: requestedSpec.Count,
						Raw:   fmt.Sprintf("%s:%d", gpu.Type, requestedSpec.Count),
					},
					Available:   gpu.Available,
					Total:       gpu.Total,
					Partition:   gpu.Partition,
					IsUpgrade:   false,
					IsDowngrade: true,
					Reason: fmt.Sprintf("Downgrade from %s (compute %.1f) to %s (compute %.1f)",
						requestedType, requestedCompat.Compute, downgradeType, downgradeCompat.Compute),
				})
			}
		}
	}

	// For MIG profiles or when no matches found, try memory-based matching
	if len(options) == 0 && knownType && requestedCompat.MemoryGB > 0 {
		memoryOptions := findGpusByMemory(requestedCompat.MemoryGB, requestedSpec.Count, info)
		options = append(options, memoryOptions...)
	}

	// If no compatible options found, suggest any available GPU
	if len(options) == 0 {
		for _, gpu := range info.AvailableGpus {
			if gpu.Available >= requestedSpec.Count {
				options = append(options, &GpuConversionOption{
					OriginalSpec: requestedSpec,
					SuggestedSpec: &GpuSpec{
						Type:  gpu.Type,
						Count: requestedSpec.Count,
						Raw:   fmt.Sprintf("%s:%d", gpu.Type, requestedSpec.Count),
					},
					Available:   gpu.Available,
					Total:       gpu.Total,
					Partition:   gpu.Partition,
					IsUpgrade:   false,
					IsDowngrade: false,
					Reason:      "No compatible GPU found, suggesting available alternative",
				})
			}
		}
	}

	// Sort options: exact match > upgrades > same tier > downgrades
	sort.Slice(options, func(i, j int) bool {
		// Prefer upgrades over downgrades
		if options[i].IsUpgrade && !options[j].IsUpgrade {
			return true
		}
		if !options[i].IsUpgrade && options[j].IsUpgrade {
			return false
		}

		// Prefer non-downgrades
		if !options[i].IsDowngrade && options[j].IsDowngrade {
			return true
		}
		if options[i].IsDowngrade && !options[j].IsDowngrade {
			return false
		}

		// Prefer higher availability
		return options[i].Available > options[j].Available
	})

	return options, nil
}

// findGpusByMemory finds GPUs that have sufficient memory (primarily for MIG profiles)
// Returns options sorted by memory size (smallest sufficient memory first)
func findGpusByMemory(requestedMemoryGB int, count int, info *ClusterInfo) []*GpuConversionOption {
	options := make([]*GpuConversionOption, 0)

	for _, gpu := range info.AvailableGpus {
		if gpu.Available < count {
			continue
		}

		normalizedType := NormalizeGpuType(gpu.Type)
		gpuCompat, found := gpuDatabase[normalizedType]

		if !found || gpuCompat.MemoryGB == 0 {
			continue
		}

		// Only suggest GPUs with enough memory
		if gpuCompat.MemoryGB >= requestedMemoryGB {
			isUpgrade := gpuCompat.MemoryGB > requestedMemoryGB
			reason := fmt.Sprintf("Memory match: %dGB requested, %dGB available",
				requestedMemoryGB, gpuCompat.MemoryGB)

			options = append(options, &GpuConversionOption{
				OriginalSpec: &GpuSpec{
					Type:  fmt.Sprintf("gpu_with_%dgb", requestedMemoryGB),
					Count: count,
				},
				SuggestedSpec: &GpuSpec{
					Type:  gpu.Type,
					Count: count,
					Raw:   fmt.Sprintf("%s:%d", gpu.Type, count),
				},
				Available:   gpu.Available,
				Total:       gpu.Total,
				Partition:   gpu.Partition,
				IsUpgrade:   isUpgrade,
				IsDowngrade: false,
				Reason:      reason,
			})
		}
	}

	// Sort by memory size (prefer smallest sufficient memory to avoid waste)
	sort.Slice(options, func(i, j int) bool {
		iType := NormalizeGpuType(options[i].SuggestedSpec.Type)
		jType := NormalizeGpuType(options[j].SuggestedSpec.Type)
		iMem := gpuDatabase[iType].MemoryGB
		jMem := gpuDatabase[jType].MemoryGB
		return iMem < jMem
	})

	return options
}

// ConvertGpuSpec converts a GPU spec to use an available GPU type
// Returns the best available option or an error if no GPUs are available
func ConvertGpuSpec(requestedSpec *GpuSpec, info *ClusterInfo) (*GpuSpec, error) {
	options, err := FindCompatibleGpu(requestedSpec, info)
	if err != nil {
		return nil, err
	}

	if len(options) == 0 {
		return nil, fmt.Errorf("%w: no compatible GPU available for %s:%d",
			ErrInvalidGpuSpec, requestedSpec.Type, requestedSpec.Count)
	}

	// Return the best option (first in sorted list)
	return options[0].SuggestedSpec, nil
}

// ValidateGpuAvailability checks if the requested GPU is available
// Returns nil if available, or an error with suggestions if not
func ValidateGpuAvailability(requestedSpec *GpuSpec, info *ClusterInfo) error {
	if requestedSpec == nil {
		return nil
	}

	options, err := FindCompatibleGpu(requestedSpec, info)
	if err != nil {
		return err
	}

	if len(options) == 0 {
		return &GpuValidationError{
			RequestedType:  requestedSpec.Type,
			RequestedCount: requestedSpec.Count,
			Suggestions:    []string{},
			AvailableGpus:  info.AvailableGpus,
		}
	}

	// Check if exact match is available
	if options[0].IsUpgrade == false && options[0].IsDowngrade == false {
		// Exact match found
		return nil
	}

	// Exact match not available, but alternatives exist
	suggestions := make([]string, 0)
	for _, opt := range options {
		suggestions = append(suggestions, fmt.Sprintf("%s (available: %d, %s)",
			opt.SuggestedSpec.Type, opt.Available, opt.Reason))
	}

	return &GpuValidationError{
		RequestedType:  requestedSpec.Type,
		RequestedCount: requestedSpec.Count,
		Suggestions:    suggestions,
		AvailableGpus:  info.AvailableGpus,
	}
}

// GetGpuInfo returns information about a specific GPU type
func GetGpuInfo(gpuType string) (GpuCompatibility, bool) {
	normalized := NormalizeGpuType(gpuType)
	info, found := gpuDatabase[normalized]
	return info, found
}

// AddCustomGpuCompat allows adding custom GPU compatibility rules
// Useful for site-specific GPU types not in the default database
func AddCustomGpuCompat(compat GpuCompatibility) {
	normalized := NormalizeGpuType(compat.Type)
	gpuDatabase[normalized] = compat
}

// IsMigProfile checks if a GPU type string is a MIG profile
// MIG profiles contain dots and underscores (e.g., nvidia_h100_80gb_hbm3_1g.10gb)
func IsMigProfile(gpuType string) bool {
	return strings.Contains(gpuType, ".") && strings.Contains(gpuType, "_")
}

// EnsureMigProfileInDatabase adds a MIG profile to the GPU database if not already present
// This is called automatically when SLURM discovers new MIG profiles
func EnsureMigProfileInDatabase(migType string) {
	normalized := NormalizeGpuType(migType)

	// Check if already in database
	if _, found := gpuDatabase[normalized]; found {
		return
	}

	// Extract memory from the MIG profile name
	// Format: nvidia_h100_80gb_hbm3_1g.10gb -> 10GB
	//         nvidia_h100_80gb_hbm3_2g.20gb -> 20GB
	//         a100_1g.5gb -> 5GB
	memoryGB := ExtractMemoryFromMigProfile(migType)

	// Determine the base GPU type (e.g., h100 from nvidia_h100_80gb_hbm3_1g.10gb)
	baseType := ExtractBaseGpuFromMigProfile(migType)
	baseCompat, baseKnown := gpuDatabase[baseType]

	tier := TierUnknown
	compute := 0.0
	upgradeTo := []string{}
	downgradeTo := []string{}

	if baseKnown {
		tier = baseCompat.Tier
		compute = baseCompat.Compute
		// Can always upgrade to the full GPU
		upgradeTo = append(upgradeTo, baseType)
	}

	// Add the MIG profile to the database
	AddCustomGpuCompat(GpuCompatibility{
		Type:        normalized,
		Tier:        tier,
		Compute:     compute,
		MemoryGB:    memoryGB,
		Aliases:     []string{},
		UpgradeTo:   upgradeTo,
		DowngradeTo: downgradeTo,
	})
}

// ExtractMemoryFromMigProfile extracts memory size in GB from a MIG profile name
// Examples:
//
//	nvidia_h100_80gb_hbm3_1g.10gb -> 10
//	a100_2g.10gb -> 10
//	nvidia_a100_40gb_1g.5gb -> 5
func ExtractMemoryFromMigProfile(migType string) int {
	// Look for pattern like "1g.10gb" or "2g.20gb"
	parts := strings.Split(migType, ".")
	if len(parts) < 2 {
		return 0
	}

	// Get the last part (e.g., "10gb")
	lastPart := strings.ToLower(parts[len(parts)-1])
	lastPart = strings.TrimSuffix(lastPart, "gb")

	memoryGB := 0
	fmt.Sscanf(lastPart, "%d", &memoryGB)
	return memoryGB
}

// ExtractBaseGpuFromMigProfile extracts the base GPU type from a MIG profile
// Examples:
//
//	nvidia_h100_80gb_hbm3_1g.10gb -> h100
//	a100_1g.5gb -> a100
//	nvidia_a30_24gb_4g.24gb -> a30
func ExtractBaseGpuFromMigProfile(migType string) string {
	lower := strings.ToLower(migType)

	// Common patterns
	if strings.Contains(lower, "h100") {
		return "h100"
	}
	if strings.Contains(lower, "a100") {
		return "a100"
	}
	if strings.Contains(lower, "a30") {
		return "a30"
	}

	return "gpu"
}
