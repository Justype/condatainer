package scheduler

import (
	"reflect"
	"testing"
)

// --- SLURM ------------------------------------------------------------------

func TestBuildSlurmDepFlag(t *testing.T) {
	tests := []struct {
		name string
		deps []Dependency
		want string
	}{
		{
			name: "no deps",
			deps: nil,
			want: "",
		},
		{
			name: "empty deps slice",
			deps: []Dependency{},
			want: "",
		},
		{
			name: "dep with no job IDs",
			deps: []Dependency{{Type: DependencyAfterOK, JobIDs: nil}},
			want: "",
		},
		{
			name: "single afterok",
			deps: []Dependency{{Type: DependencyAfterOK, JobIDs: []string{"123"}}},
			want: "--dependency=afterok:123",
		},
		{
			name: "single afternotok",
			deps: []Dependency{{Type: DependencyAfterNotOK, JobIDs: []string{"456"}}},
			want: "--dependency=afternotok:456",
		},
		{
			name: "single afterany",
			deps: []Dependency{{Type: DependencyAfterAny, JobIDs: []string{"789"}}},
			want: "--dependency=afterany:789",
		},
		{
			name: "afterok multiple IDs",
			deps: []Dependency{{Type: DependencyAfterOK, JobIDs: []string{"1", "2", "3"}}},
			want: "--dependency=afterok:1:2:3",
		},
		{
			name: "mixed types",
			deps: []Dependency{
				{Type: DependencyAfterOK, JobIDs: []string{"100", "101"}},
				{Type: DependencyAfterNotOK, JobIDs: []string{"200"}},
				{Type: DependencyAfterAny, JobIDs: []string{"300"}},
			},
			want: "--dependency=afterok:100:101,afternotok:200,afterany:300",
		},
		{
			name: "afternotok and afterany only",
			deps: []Dependency{
				{Type: DependencyAfterNotOK, JobIDs: []string{"10"}},
				{Type: DependencyAfterAny, JobIDs: []string{"20"}},
			},
			want: "--dependency=afternotok:10,afterany:20",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := buildSlurmDepFlag(tt.deps)
			if got != tt.want {
				t.Errorf("buildSlurmDepFlag() = %q; want %q", got, tt.want)
			}
		})
	}
}

// --- PBS --------------------------------------------------------------------

func TestBuildPbsDepFlag(t *testing.T) {
	tests := []struct {
		name string
		deps []Dependency
		want string
	}{
		{
			name: "no deps",
			deps: nil,
			want: "",
		},
		{
			name: "empty deps slice",
			deps: []Dependency{},
			want: "",
		},
		{
			name: "dep with no job IDs",
			deps: []Dependency{{Type: DependencyAfterOK, JobIDs: nil}},
			want: "",
		},
		{
			name: "single afterok",
			deps: []Dependency{{Type: DependencyAfterOK, JobIDs: []string{"123"}}},
			want: "-W depend=afterok:123",
		},
		{
			name: "single afternotok",
			deps: []Dependency{{Type: DependencyAfterNotOK, JobIDs: []string{"456"}}},
			want: "-W depend=afternotok:456",
		},
		{
			name: "single afterany",
			deps: []Dependency{{Type: DependencyAfterAny, JobIDs: []string{"789"}}},
			want: "-W depend=afterany:789",
		},
		{
			name: "afterok multiple IDs (colon-separated)",
			deps: []Dependency{{Type: DependencyAfterOK, JobIDs: []string{"1", "2", "3"}}},
			want: "-W depend=afterok:1:2:3",
		},
		{
			name: "mixed types",
			deps: []Dependency{
				{Type: DependencyAfterOK, JobIDs: []string{"100", "101"}},
				{Type: DependencyAfterNotOK, JobIDs: []string{"200"}},
				{Type: DependencyAfterAny, JobIDs: []string{"300"}},
			},
			want: "-W depend=afterok:100:101,afternotok:200,afterany:300",
		},
		{
			name: "afternotok and afterany only",
			deps: []Dependency{
				{Type: DependencyAfterNotOK, JobIDs: []string{"10"}},
				{Type: DependencyAfterAny, JobIDs: []string{"20"}},
			},
			want: "-W depend=afternotok:10,afterany:20",
		},
		{
			name: "PBS job IDs with host suffix",
			deps: []Dependency{{Type: DependencyAfterOK, JobIDs: []string{"12345.host.example.com"}}},
			want: "-W depend=afterok:12345.host.example.com",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := buildPbsDepFlag(tt.deps)
			if got != tt.want {
				t.Errorf("buildPbsDepFlag() = %q; want %q", got, tt.want)
			}
		})
	}
}

// --- LSF --------------------------------------------------------------------

func TestBuildLsfDepCondition(t *testing.T) {
	tests := []struct {
		name string
		deps []Dependency
		want string
	}{
		{
			name: "no deps",
			deps: nil,
			want: "",
		},
		{
			name: "empty deps slice",
			deps: []Dependency{},
			want: "",
		},
		{
			name: "dep with no job IDs",
			deps: []Dependency{{Type: DependencyAfterOK, JobIDs: nil}},
			want: "",
		},
		{
			name: "afterok → done()",
			deps: []Dependency{{Type: DependencyAfterOK, JobIDs: []string{"123"}}},
			want: "done(123)",
		},
		{
			name: "afternotok → exit()",
			deps: []Dependency{{Type: DependencyAfterNotOK, JobIDs: []string{"456"}}},
			want: "exit(456)",
		},
		{
			name: "afterany → ended()",
			deps: []Dependency{{Type: DependencyAfterAny, JobIDs: []string{"789"}}},
			want: "ended(789)",
		},
		{
			name: "afterok multiple IDs (separate done() per ID)",
			deps: []Dependency{{Type: DependencyAfterOK, JobIDs: []string{"1", "2", "3"}}},
			want: "done(1) && done(2) && done(3)",
		},
		{
			name: "mixed types",
			deps: []Dependency{
				{Type: DependencyAfterOK, JobIDs: []string{"100", "101"}},
				{Type: DependencyAfterNotOK, JobIDs: []string{"200"}},
				{Type: DependencyAfterAny, JobIDs: []string{"300"}},
			},
			want: "done(100) && done(101) && exit(200) && ended(300)",
		},
		{
			name: "afternotok and afterany only",
			deps: []Dependency{
				{Type: DependencyAfterNotOK, JobIDs: []string{"10"}},
				{Type: DependencyAfterAny, JobIDs: []string{"20"}},
			},
			want: "exit(10) && ended(20)",
		},
		{
			name: "unknown type defaults to done()",
			deps: []Dependency{{Type: "unknown", JobIDs: []string{"99"}}},
			want: "done(99)",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := buildLsfDepCondition(tt.deps)
			if got != tt.want {
				t.Errorf("buildLsfDepCondition() = %q; want %q", got, tt.want)
			}
		})
	}
}

// --- SLURM submit args (--kill-on-invalid-dep) -----------------------------------

func TestBuildSlurmSubmitArgs(t *testing.T) {
	const script = "job.sh"
	tests := []struct {
		name string
		deps []Dependency
		want []string
	}{
		{
			name: "no deps — no kill-on-invalid-dep",
			deps: nil,
			want: []string{script},
		},
		{
			name: "empty deps — no kill-on-invalid-dep",
			deps: []Dependency{},
			want: []string{script},
		},
		{
			name: "dep with no IDs — no kill-on-invalid-dep",
			deps: []Dependency{{Type: DependencyAfterOK, JobIDs: nil}},
			want: []string{script},
		},
		{
			name: "afterok dep — includes kill-on-invalid-dep",
			deps: []Dependency{{Type: DependencyAfterOK, JobIDs: []string{"123"}}},
			want: []string{"--dependency=afterok:123", "--kill-on-invalid-dep=yes", script},
		},
		{
			name: "mixed deps — includes kill-on-invalid-dep",
			deps: []Dependency{
				{Type: DependencyAfterOK, JobIDs: []string{"1", "2"}},
				{Type: DependencyAfterNotOK, JobIDs: []string{"3"}},
			},
			want: []string{"--dependency=afterok:1:2,afternotok:3", "--kill-on-invalid-dep=yes", script},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := buildSlurmSubmitArgs(tt.deps, script)
			if !reflect.DeepEqual(got, tt.want) {
				t.Errorf("buildSlurmSubmitArgs() = %v; want %v", got, tt.want)
			}
		})
	}
}

// --- LSF submit args (-ti) -------------------------------------------------------

func TestBuildLsfArgs(t *testing.T) {
	tests := []struct {
		name string
		deps []Dependency
		want []string
	}{
		{
			name: "no deps — no -ti",
			deps: nil,
			want: nil,
		},
		{
			name: "empty deps — no -ti",
			deps: []Dependency{},
			want: nil,
		},
		{
			name: "dep with no IDs — no -ti",
			deps: []Dependency{{Type: DependencyAfterOK, JobIDs: nil}},
			want: nil,
		},
		{
			name: "afterok dep — includes -ti",
			deps: []Dependency{{Type: DependencyAfterOK, JobIDs: []string{"123"}}},
			want: []string{"-w", "done(123)", "-ti"},
		},
		{
			name: "mixed deps — includes -ti",
			deps: []Dependency{
				{Type: DependencyAfterOK, JobIDs: []string{"100"}},
				{Type: DependencyAfterNotOK, JobIDs: []string{"200"}},
				{Type: DependencyAfterAny, JobIDs: []string{"300"}},
			},
			want: []string{"-w", "done(100) && exit(200) && ended(300)", "-ti"},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := buildLsfArgs(tt.deps)
			if !reflect.DeepEqual(got, tt.want) {
				t.Errorf("buildLsfArgs() = %v; want %v", got, tt.want)
			}
		})
	}
}
