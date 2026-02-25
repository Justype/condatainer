package scheduler

import (
	"fmt"
	"testing"
)

func TestNormalizeGpuType(t *testing.T) {
	// Dynamically add MIG profiles for testing (simulates SLURM discovery)
	migProfiles := []string{
		"nvidia_h100_80gb_hbm3_1g.10gb",
		"nvidia_h100_80gb_hbm3_2g.20gb",
		"nvidia_h100_80gb_hbm3_3g.40gb",
	}
	for _, profile := range migProfiles {
		if IsMigProfile(profile) {
			EnsureMigProfileInDatabase(profile)
		}
	}

	tests := []struct {
		input    string
		expected string
	}{
		{"a100", "a100"},
		{"A100", "a100"},
		{"NVIDIA-A100", "a100"},
		{"nvidia_a100", "a100"},
		{"a100-sxm", "a100"},
		{"a100-pcie", "a100"},
		{"V100", "v100"},
		{"tesla-v100", "v100"},
		{"h100", "h100"},
		{"Hopper", "h100"},
		// MIG profiles (dynamically discovered)
		{"nvidia_h100_80gb_hbm3_1g.10gb", "nvidia_h100_80gb_hbm3_1g.10gb"},
		{"NVIDIA_H100_80GB_HBM3_1G.10GB", "nvidia_h100_80gb_hbm3_1g.10gb"},
		{"nvidia_h100_80gb_hbm3_2g.20gb", "nvidia_h100_80gb_hbm3_2g.20gb"},
		{"nvidia_h100_80gb_hbm3_3g.40gb", "nvidia_h100_80gb_hbm3_3g.40gb"},
	}

	for _, tt := range tests {
		t.Run(tt.input, func(t *testing.T) {
			result := NormalizeGpuType(tt.input)
			if result != tt.expected {
				t.Errorf("NormalizeGpuType(%q) = %q; want %q", tt.input, result, tt.expected)
			}
		})
	}
}

func TestGetGpuInfo(t *testing.T) {
	tests := []struct {
		gpuType string
		wantOk  bool
	}{
		{"a100", true},
		{"h100", true},
		{"v100", true},
		{"p100", true},
		{"t4", true},
		{"k80", true},
		{"unknown_gpu", false},
	}

	for _, tt := range tests {
		t.Run(tt.gpuType, func(t *testing.T) {
			info, ok := GetGpuInfo(tt.gpuType)
			if ok != tt.wantOk {
				t.Errorf("GetGpuInfo(%q) ok = %v; want %v", tt.gpuType, ok, tt.wantOk)
			}
			if ok && info.Type != tt.gpuType {
				t.Errorf("GetGpuInfo(%q).Type = %q; want %q", tt.gpuType, info.Type, tt.gpuType)
			}
		})
	}
}

func TestFindCompatibleGpu(t *testing.T) {
	// Mock cluster info with various GPUs
	clusterInfo := &ClusterInfo{
		AvailableGpus: []GpuInfo{
			{Type: "h100", Total: 8, Available: 4, Partition: "gpu-high"},
			{Type: "a100", Total: 16, Available: 8, Partition: "gpu-standard"},
			{Type: "v100", Total: 32, Available: 16, Partition: "gpu-standard"},
			{Type: "p100", Total: 24, Available: 0, Partition: "gpu-legacy"},
		},
	}

	tests := []struct {
		name         string
		requestedGpu *GpuSpec
		wantMinOpts  int    // Minimum number of options expected
		wantFirst    string // Expected first (best) option
		wantUpgrade  bool
	}{
		{
			name:         "Exact match available (a100)",
			requestedGpu: &GpuSpec{Type: "a100", Count: 2},
			wantMinOpts:  1,
			wantFirst:    "a100",
			wantUpgrade:  false,
		},
		{
			name:         "Exact match and upgrades available (v100)",
			requestedGpu: &GpuSpec{Type: "v100", Count: 2},
			wantMinOpts:  1, // At least v100 exact match
			wantFirst:    "v100",
			wantUpgrade:  false,
		},
		{
			name:         "No exact match, upgrade available",
			requestedGpu: &GpuSpec{Type: "p100", Count: 2},
			wantMinOpts:  1, // v100 upgrade available
			wantFirst:    "v100",
			wantUpgrade:  true,
		},
		{
			name:         "Request more than available",
			requestedGpu: &GpuSpec{Type: "h100", Count: 8},
			wantMinOpts:  0, // Need 8 but only 4 available
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			options, err := FindCompatibleGpu(tt.requestedGpu, 1, clusterInfo)
			if err != nil {
				t.Fatalf("FindCompatibleGpu() error = %v", err)
			}

			if len(options) < tt.wantMinOpts {
				t.Errorf("FindCompatibleGpu() returned %d options; want at least %d", len(options), tt.wantMinOpts)
			}

			if len(options) > 0 && tt.wantFirst != "" {
				firstType := NormalizeGpuType(options[0].SuggestedSpec.Type)
				wantType := NormalizeGpuType(tt.wantFirst)
				if firstType != wantType {
					t.Errorf("First option type = %q; want %q", firstType, wantType)
				}

				if options[0].IsUpgrade != tt.wantUpgrade {
					t.Errorf("First option IsUpgrade = %v; want %v", options[0].IsUpgrade, tt.wantUpgrade)
				}
			}
		})
	}
}

func TestConvertGpuSpec(t *testing.T) {
	clusterInfo := &ClusterInfo{
		AvailableGpus: []GpuInfo{
			{Type: "a100", Total: 16, Available: 8, Partition: "gpu"},
			{Type: "v100", Total: 32, Available: 16, Partition: "gpu"},
		},
	}

	tests := []struct {
		name         string
		requestedGpu *GpuSpec
		wantType     string
		wantError    bool
	}{
		{
			name:         "Exact match",
			requestedGpu: &GpuSpec{Type: "a100", Count: 2},
			wantType:     "a100",
			wantError:    false,
		},
		{
			name:         "Convert p100 to v100 (upgrade)",
			requestedGpu: &GpuSpec{Type: "p100", Count: 4},
			wantType:     "v100",
			wantError:    false,
		},
		{
			name:         "No GPU available for count",
			requestedGpu: &GpuSpec{Type: "a100", Count: 20},
			wantError:    true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result, err := ConvertGpuSpec(tt.requestedGpu, 1, clusterInfo)

			if tt.wantError {
				if err == nil {
					t.Errorf("ConvertGpuSpec() expected error but got none")
				}
				return
			}

			if err != nil {
				t.Fatalf("ConvertGpuSpec() unexpected error = %v", err)
			}

			resultType := NormalizeGpuType(result.Type)
			wantType := NormalizeGpuType(tt.wantType)
			if resultType != wantType {
				t.Errorf("ConvertGpuSpec() type = %q; want %q", resultType, wantType)
			}

			if result.Count != tt.requestedGpu.Count {
				t.Errorf("ConvertGpuSpec() count = %d; want %d", result.Count, tt.requestedGpu.Count)
			}
		})
	}
}

// TestCrossGpuTypeConversion tests conversion between different GPU types
func TestCrossGpuTypeConversion(t *testing.T) {
	tests := []struct {
		name         string
		clusterInfo  *ClusterInfo
		requestedGpu *GpuSpec
		wantType     string
		wantError    bool
		wantMinOpts  int
	}{
		{
			name: "Request a100 but only h100 available - should upgrade",
			clusterInfo: &ClusterInfo{
				AvailableGpus: []GpuInfo{
					{Type: "h100", Total: 8, Available: 4, Partition: "gpu"},
				},
			},
			requestedGpu: &GpuSpec{Type: "a100", Count: 2},
			wantType:     "h100",
			wantError:    false,
			wantMinOpts:  1,
		},
		{
			name: "Request h100 but only a100 available - should downgrade",
			clusterInfo: &ClusterInfo{
				AvailableGpus: []GpuInfo{
					{Type: "a100", Total: 16, Available: 8, Partition: "gpu"},
				},
			},
			requestedGpu: &GpuSpec{Type: "h100", Count: 2},
			wantType:     "a100",
			wantError:    false,
			wantMinOpts:  1,
		},
		{
			name: "Request a100 with multiple options available",
			clusterInfo: &ClusterInfo{
				AvailableGpus: []GpuInfo{
					{Type: "h100", Total: 8, Available: 4, Partition: "gpu-high"},
					{Type: "a40", Total: 16, Available: 8, Partition: "gpu-mid"},
					{Type: "v100", Total: 32, Available: 16, Partition: "gpu-legacy"},
				},
			},
			requestedGpu: &GpuSpec{Type: "a100", Count: 2},
			wantType:     "h100", // Should prefer h100 (upgrade) over a40 or v100
			wantError:    false,
			wantMinOpts:  1,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// Test FindCompatibleGpu
			options, err := FindCompatibleGpu(tt.requestedGpu, 1, tt.clusterInfo)
			if err != nil {
				t.Fatalf("FindCompatibleGpu() error = %v", err)
			}

			if len(options) < tt.wantMinOpts {
				t.Errorf("FindCompatibleGpu() returned %d options; want at least %d", len(options), tt.wantMinOpts)
			}

			if len(options) > 0 {
				firstType := NormalizeGpuType(options[0].SuggestedSpec.Type)
				wantType := NormalizeGpuType(tt.wantType)
				if firstType != wantType {
					t.Errorf("First option type = %q; want %q", firstType, wantType)
				}

				t.Logf("Requested: %s x%d -> Suggested: %s x%d (%s)",
					tt.requestedGpu.Type, tt.requestedGpu.Count,
					options[0].SuggestedSpec.Type, options[0].SuggestedSpec.Count,
					options[0].Reason)
			}

			// Test ConvertGpuSpec
			result, err := ConvertGpuSpec(tt.requestedGpu, 1, tt.clusterInfo)

			if tt.wantError {
				if err == nil {
					t.Errorf("ConvertGpuSpec() expected error but got none")
				}
				return
			}

			if err != nil {
				t.Fatalf("ConvertGpuSpec() unexpected error = %v", err)
			}

			resultType := NormalizeGpuType(result.Type)
			wantType := NormalizeGpuType(tt.wantType)
			if resultType != wantType {
				t.Errorf("ConvertGpuSpec() type = %q; want %q", resultType, wantType)
			}

			if result.Count != tt.requestedGpu.Count {
				t.Errorf("ConvertGpuSpec() count = %d; want %d", result.Count, tt.requestedGpu.Count)
			}
		})
	}
}

func TestConvertGpuSpecOriginal(t *testing.T) {
	clusterInfo := &ClusterInfo{
		AvailableGpus: []GpuInfo{
			{Type: "a100", Total: 16, Available: 8, Partition: "gpu"},
			{Type: "v100", Total: 32, Available: 16, Partition: "gpu"},
		},
	}

	tests := []struct {
		name         string
		requestedGpu *GpuSpec
		wantType     string
		wantError    bool
	}{
		{
			name:         "Exact match",
			requestedGpu: &GpuSpec{Type: "a100", Count: 2},
			wantType:     "a100",
			wantError:    false,
		},
		{
			name:         "Convert p100 to v100 (upgrade)",
			requestedGpu: &GpuSpec{Type: "p100", Count: 4},
			wantType:     "v100",
			wantError:    false,
		},
		{
			name:         "No GPU available for count",
			requestedGpu: &GpuSpec{Type: "a100", Count: 20},
			wantError:    true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result, err := ConvertGpuSpec(tt.requestedGpu, 1, clusterInfo)

			if tt.wantError {
				if err == nil {
					t.Errorf("ConvertGpuSpec() expected error but got none")
				}
				return
			}

			if err != nil {
				t.Fatalf("ConvertGpuSpec() unexpected error = %v", err)
			}

			resultType := NormalizeGpuType(result.Type)
			wantType := NormalizeGpuType(tt.wantType)
			if resultType != wantType {
				t.Errorf("ConvertGpuSpec() type = %q; want %q", resultType, wantType)
			}

			if result.Count != tt.requestedGpu.Count {
				t.Errorf("ConvertGpuSpec() count = %d; want %d", result.Count, tt.requestedGpu.Count)
			}
		})
	}
}

func TestValidateGpuAvailability(t *testing.T) {
	clusterInfo := &ClusterInfo{
		AvailableGpus: []GpuInfo{
			{Type: "a100", Total: 16, Available: 8, Partition: "gpu"},
			{Type: "v100", Total: 32, Available: 0, Partition: "gpu"},
		},
	}

	tests := []struct {
		name         string
		requestedGpu *GpuSpec
		wantError    bool
		wantGpuError bool
	}{
		{
			name:         "Available GPU",
			requestedGpu: &GpuSpec{Type: "a100", Count: 4},
			wantError:    false,
		},
		{
			name:         "GPU exists but none available",
			requestedGpu: &GpuSpec{Type: "v100", Count: 1},
			wantError:    true,
			wantGpuError: true,
		},
		{
			name:         "Unknown GPU type with available alternatives",
			requestedGpu: &GpuSpec{Type: "unknown", Count: 1},
			wantError:    false, // Unknown GPUs will suggest any available GPU
			wantGpuError: false,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			err := ValidateGpuAvailability(tt.requestedGpu, 1, clusterInfo)

			if tt.wantError {
				if err == nil {
					t.Errorf("ValidateGpuAvailability() expected error but got none")
				}
				if tt.wantGpuError && !IsGpuValidationError(err) {
					t.Errorf("ValidateGpuAvailability() error is not GpuValidationError: %v", err)
				}
			} else {
				if err != nil {
					t.Errorf("ValidateGpuAvailability() unexpected error = %v", err)
				}
			}
		})
	}
}

func TestAddCustomGpuCompat(t *testing.T) {
	// Add a custom GPU
	customGpu := GpuCompatibility{
		Type:        "mi250x",
		Tier:        TierHigh,
		Compute:     9.0,
		MemoryGB:    128,
		Aliases:     []string{"mi250x", "mi250"},
		UpgradeTo:   []string{},
		DowngradeTo: []string{"v100"},
	}

	AddCustomGpuCompat(customGpu)

	// Verify it was added
	info, found := GetGpuInfo("mi250x")
	if !found {
		t.Errorf("Custom GPU 'mi250x' not found after adding")
	}

	if info.Type != "mi250x" {
		t.Errorf("Custom GPU type = %q; want 'mi250x'", info.Type)
	}

	if info.MemoryGB != 128 {
		t.Errorf("Custom GPU memory = %d; want 128", info.MemoryGB)
	}
}

func TestFindAvailableGpuTypes(t *testing.T) {
	clusterInfo := &ClusterInfo{
		AvailableGpus: []GpuInfo{
			{Type: "a100", Total: 16, Available: 8, Partition: "gpu"},
			{Type: "v100", Total: 32, Available: 0, Partition: "gpu"},
			{Type: "h100", Total: 8, Available: 4, Partition: "gpu-high"},
		},
	}

	available := FindAvailableGpuTypes(clusterInfo)

	// Should only return GPUs with Available > 0
	if len(available) != 2 {
		t.Errorf("FindAvailableGpuTypes() returned %d GPUs; want 2", len(available))
	}

	// Check that v100 (with 0 available) is not in the list
	for _, gpu := range available {
		if gpu.Type == "v100" {
			t.Errorf("FindAvailableGpuTypes() included v100 with 0 available")
		}
	}
}

// Example function demonstrating GPU conversion workflow
func ExampleFindCompatibleGpu() {
	clusterInfo := &ClusterInfo{
		AvailableGpus: []GpuInfo{
			{Type: "h100", Total: 8, Available: 4, Partition: "gpu-high"},
			{Type: "a100", Total: 16, Available: 8, Partition: "gpu-standard"},
			{Type: "v100", Total: 32, Available: 16, Partition: "gpu-standard"},
		},
	}

	requestedSpec := &GpuSpec{
		Type:  "v100",
		Count: 2,
	}

	options, err := FindCompatibleGpu(requestedSpec, 1, clusterInfo)
	if err != nil {
		fmt.Printf("Error: %v\n", err)
		return
	}

	fmt.Printf("Found %d compatible GPU options:\n", len(options))
	for i, opt := range options {
		status := "MATCH"
		if opt.IsUpgrade {
			status = "UPGRADE"
		} else if opt.IsDowngrade {
			status = "DOWNGRADE"
		}

		fmt.Printf("%d. %s x%d [%s] - %d available\n",
			i+1, opt.SuggestedSpec.Type, opt.SuggestedSpec.Count, status, opt.Available)
	}

	// Output:
	// Found 1 compatible GPU options:
	// 1. v100 x2 [MATCH] - 16 available
}

// Example demonstrating automatic GPU conversion
func ExampleConvertGpuSpec() {
	clusterInfo := &ClusterInfo{
		AvailableGpus: []GpuInfo{
			{Type: "a100", Total: 16, Available: 8, Partition: "gpu"},
			{Type: "v100", Total: 32, Available: 16, Partition: "gpu"},
		},
	}

	requestedSpec := &GpuSpec{
		Type:  "p100", // Not available
		Count: 2,
	}

	convertedSpec, err := ConvertGpuSpec(requestedSpec, 1, clusterInfo)
	if err != nil {
		fmt.Printf("Error: %v\n", err)
		return
	}

	fmt.Printf("Requested: %s x%d\n", requestedSpec.Type, requestedSpec.Count)
	fmt.Printf("Converted to: %s x%d\n", convertedSpec.Type, convertedSpec.Count)

	// Output:
	// Requested: p100 x2
	// Converted to: v100 x2
}

// TestParseSlurmGpuMIG tests parsing of MIG GPU specifications
func TestParseSlurmGpuMIG(t *testing.T) {
	tests := []struct {
		name      string
		input     string
		wantType  string
		wantCount int
		wantErr   bool
	}{
		{
			name:      "MIG profile without count",
			input:     "nvidia_h100_80gb_hbm3_1g.10gb",
			wantType:  "nvidia_h100_80gb_hbm3_1g.10gb",
			wantCount: 1,
			wantErr:   false,
		},
		{
			name:      "MIG profile with count",
			input:     "nvidia_h100_80gb_hbm3_1g.10gb:2",
			wantType:  "nvidia_h100_80gb_hbm3_1g.10gb",
			wantCount: 2,
			wantErr:   false,
		},
		{
			name:      "MIG 2g.20gb profile",
			input:     "nvidia_h100_80gb_hbm3_2g.20gb:3",
			wantType:  "nvidia_h100_80gb_hbm3_2g.20gb",
			wantCount: 3,
			wantErr:   false,
		},
		{
			name:      "MIG 3g.40gb profile",
			input:     "nvidia_h100_80gb_hbm3_3g.40gb:1",
			wantType:  "nvidia_h100_80gb_hbm3_3g.40gb",
			wantCount: 1,
			wantErr:   false,
		},
		{
			name:      "Regular GPU a100",
			input:     "a100:2",
			wantType:  "a100",
			wantCount: 2,
			wantErr:   false,
		},
		{
			name:      "GPU gres format",
			input:     "gpu:a100:4",
			wantType:  "a100",
			wantCount: 4,
			wantErr:   false,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result, err := parseSlurmGpu(tt.input)

			if tt.wantErr {
				if err == nil {
					t.Errorf("parseSlurmGpu(%q) expected error but got none", tt.input)
				}
				return
			}

			if err != nil {
				t.Fatalf("parseSlurmGpu(%q) unexpected error = %v", tt.input, err)
			}

			if result.Type != tt.wantType {
				t.Errorf("parseSlurmGpu(%q) type = %q; want %q", tt.input, result.Type, tt.wantType)
			}

			if result.Count != tt.wantCount {
				t.Errorf("parseSlurmGpu(%q) count = %d; want %d", tt.input, result.Count, tt.wantCount)
			}
		})
	}
}

// TestGetGpuInfoMIG tests getting GPU info for MIG profiles
func TestGetGpuInfoMIG(t *testing.T) {
	// Dynamically add MIG profiles for testing (simulates SLURM discovery)
	migProfiles := []string{
		"nvidia_h100_80gb_hbm3_1g.10gb",
		"nvidia_h100_80gb_hbm3_2g.20gb",
		"nvidia_h100_80gb_hbm3_3g.40gb",
	}
	for _, profile := range migProfiles {
		if IsMigProfile(profile) {
			EnsureMigProfileInDatabase(profile)
		}
	}

	tests := []struct {
		gpuType    string
		wantOk     bool
		wantMemory int
	}{
		{"nvidia_h100_80gb_hbm3_1g.10gb", true, 10},
		{"nvidia_h100_80gb_hbm3_2g.20gb", true, 20},
		{"nvidia_h100_80gb_hbm3_3g.40gb", true, 40},
	}

	for _, tt := range tests {
		t.Run(tt.gpuType, func(t *testing.T) {
			info, ok := GetGpuInfo(tt.gpuType)
			if ok != tt.wantOk {
				t.Errorf("GetGpuInfo(%q) ok = %v; want %v", tt.gpuType, ok, tt.wantOk)
			}
			if ok && info.MemoryGB != tt.wantMemory {
				t.Errorf("GetGpuInfo(%q).MemoryGB = %d; want %d", tt.gpuType, info.MemoryGB, tt.wantMemory)
			}
		})
	}
}

// TestA100MigProfiles tests A100 MIG profiles with shorter naming format
func TestA100MigProfiles(t *testing.T) {
	// Test parsing A100 MIG profiles with shorter naming convention
	tests := []struct {
		name         string
		input        string
		wantType     string
		wantCount    int
		wantMemoryGB int
		wantErr      bool
	}{
		{
			name:         "A100 1g.5gb profile",
			input:        "a100_1g.5gb",
			wantType:     "a100_1g.5gb",
			wantCount:    1,
			wantMemoryGB: 5,
			wantErr:      false,
		},
		{
			name:         "A100 1g.5gb with count",
			input:        "a100_1g.5gb:2",
			wantType:     "a100_1g.5gb",
			wantCount:    2,
			wantMemoryGB: 5,
			wantErr:      false,
		},
		{
			name:         "A100 2g.10gb profile",
			input:        "a100_2g.10gb",
			wantType:     "a100_2g.10gb",
			wantCount:    1,
			wantMemoryGB: 10,
			wantErr:      false,
		},
		{
			name:         "A100 2g.10gb with count",
			input:        "a100_2g.10gb:1",
			wantType:     "a100_2g.10gb",
			wantCount:    1,
			wantMemoryGB: 10,
			wantErr:      false,
		},
		{
			name:         "A100 3g.20gb profile",
			input:        "a100_3g.20gb",
			wantType:     "a100_3g.20gb",
			wantCount:    1,
			wantMemoryGB: 20,
			wantErr:      false,
		},
		{
			name:         "A100 3g.20gb with count",
			input:        "a100_3g.20gb:1",
			wantType:     "a100_3g.20gb",
			wantCount:    1,
			wantMemoryGB: 20,
			wantErr:      false,
		},
		{
			name:         "A100 4g.20gb profile",
			input:        "a100_4g.20gb",
			wantType:     "a100_4g.20gb",
			wantCount:    1,
			wantMemoryGB: 20,
			wantErr:      false,
		},
		{
			name:         "A100 4g.20gb with count",
			input:        "a100_4g.20gb:1",
			wantType:     "a100_4g.20gb",
			wantCount:    1,
			wantMemoryGB: 20,
			wantErr:      false,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// Test parsing
			result, err := parseSlurmGpu(tt.input)

			if tt.wantErr {
				if err == nil {
					t.Errorf("parseSlurmGpu(%q) expected error but got none", tt.input)
				}
				return
			}

			if err != nil {
				t.Fatalf("parseSlurmGpu(%q) unexpected error = %v", tt.input, err)
			}

			if result.Type != tt.wantType {
				t.Errorf("parseSlurmGpu(%q) type = %q; want %q", tt.input, result.Type, tt.wantType)
			}

			if result.Count != tt.wantCount {
				t.Errorf("parseSlurmGpu(%q) count = %d; want %d", tt.input, result.Count, tt.wantCount)
			}

			// Test dynamic discovery
			if IsMigProfile(tt.wantType) {
				EnsureMigProfileInDatabase(tt.wantType)
			}

			// Verify it was added to database
			info, found := GetGpuInfo(tt.wantType)
			if !found {
				t.Errorf("MIG profile %q not found in database after discovery", tt.wantType)
			}

			// Verify memory was extracted correctly
			if found && info.MemoryGB != tt.wantMemoryGB {
				t.Errorf("MIG profile %q memory = %dGB; want %dGB", tt.wantType, info.MemoryGB, tt.wantMemoryGB)
			}

			// Verify base GPU type was identified
			if found && info.Tier != TierHigh {
				t.Errorf("MIG profile %q tier = %v; want TierHigh (inherited from a100)", tt.wantType, info.Tier)
			}

			t.Logf("✓ %s: parsed=%s x%d, memory=%dGB, tier=%v",
				tt.name, result.Type, result.Count, info.MemoryGB, info.Tier)
		})
	}
}

// TestDynamicMigProfileDetection tests dynamic MIG profile discovery
func TestDynamicMigProfileDetection(t *testing.T) {
	// Simulate discovering a new MIG profile not in the hardcoded database
	newMigProfile := "nvidia_a100_80gb_pcie_2g.20gb"

	// Before discovery, it shouldn't be in the database
	_, foundBefore := GetGpuInfo(newMigProfile)

	// Simulate SLURM discovery
	if IsMigProfile(newMigProfile) {
		EnsureMigProfileInDatabase(newMigProfile)
	}

	// After discovery, it should be in the database
	info, foundAfter := GetGpuInfo(newMigProfile)
	if !foundAfter {
		t.Errorf("Dynamic MIG profile '%s' was not added to database", newMigProfile)
	}

	// Check that memory was extracted correctly
	if foundAfter && info.MemoryGB != 20 {
		t.Errorf("MIG profile memory = %dGB; want 20GB", info.MemoryGB)
	}

	// Check that base GPU type was identified
	if foundAfter && info.Tier != TierHigh {
		t.Errorf("MIG profile tier = %v; want TierHigh (inherited from a100)", info.Tier)
	}

	t.Logf("Before discovery: found=%v, After discovery: found=%v, memory=%dGB",
		foundBefore, foundAfter, info.MemoryGB)
}

// TestFindCompatibleGpuMIG tests finding compatible MIG GPUs
func TestFindCompatibleGpuMIG(t *testing.T) {
	// Mock cluster info with MIG profiles
	clusterInfo := &ClusterInfo{
		AvailableGpus: []GpuInfo{
			{Type: "nvidia_h100_80gb_hbm3_1g.10gb", Total: 8, Available: 4, Partition: "gpu-mig"},
			{Type: "nvidia_h100_80gb_hbm3_2g.20gb", Total: 4, Available: 2, Partition: "gpu-mig"},
			{Type: "nvidia_h100_80gb_hbm3_3g.40gb", Total: 2, Available: 1, Partition: "gpu-mig"},
			{Type: "h100", Total: 1, Available: 1, Partition: "gpu-full"},
		},
	}

	tests := []struct {
		name         string
		requestedGpu *GpuSpec
		wantMinOpts  int
		wantUpgrade  bool
	}{
		{
			name:         "Request 1g.10gb MIG",
			requestedGpu: &GpuSpec{Type: "nvidia_h100_80gb_hbm3_1g.10gb", Count: 2},
			wantMinOpts:  1,
			wantUpgrade:  false,
		},
		{
			name:         "Request 2g.20gb MIG - can upgrade to 3g or full",
			requestedGpu: &GpuSpec{Type: "nvidia_h100_80gb_hbm3_2g.20gb", Count: 1},
			wantMinOpts:  1,
			wantUpgrade:  false,
		},
		{
			name:         "Request 1g.10gb - can upgrade to larger MIG or full H100",
			requestedGpu: &GpuSpec{Type: "nvidia_h100_80gb_hbm3_1g.10gb", Count: 1},
			wantMinOpts:  1,
			wantUpgrade:  false,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			options, err := FindCompatibleGpu(tt.requestedGpu, 1, clusterInfo)
			if err != nil {
				t.Fatalf("FindCompatibleGpu() error = %v", err)
			}

			if len(options) < tt.wantMinOpts {
				t.Errorf("FindCompatibleGpu() returned %d options; want at least %d", len(options), tt.wantMinOpts)
			}
		})
	}
}
