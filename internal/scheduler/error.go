package scheduler

import (
	"errors"
	"fmt"
	"strings"
)

// Common errors
var (
	// ErrSchedulerNotAvailable indicates the scheduler is not available
	ErrSchedulerNotAvailable = errors.New("scheduler is not available")

	// ErrSchedulerNotFound indicates the scheduler binary was not found
	ErrSchedulerNotFound = errors.New("scheduler binary not found in PATH")

	// ErrAlreadyInJob indicates we're already inside a scheduler job
	ErrAlreadyInJob = errors.New("already inside a scheduler job")

	// ErrScriptNotFound indicates the script file was not found
	ErrScriptNotFound = errors.New("script file not found")

	// ErrInvalidScriptFormat indicates the script has invalid format
	ErrInvalidScriptFormat = errors.New("invalid script format")

	// ErrJobSubmissionFailed indicates job submission failed
	ErrJobSubmissionFailed = errors.New("job submission failed")

	// ErrJobIDParseFailed indicates parsing job ID from output failed
	ErrJobIDParseFailed = errors.New("failed to parse job ID from scheduler output")

	// ErrClusterInfoUnavailable indicates cluster information is not available
	ErrClusterInfoUnavailable = errors.New("cluster information unavailable")

	// ErrInvalidGpuSpec indicates GPU specification is invalid
	ErrInvalidGpuSpec = errors.New("invalid GPU specification")

	// ErrInvalidTimeFormat indicates time format is invalid
	ErrInvalidTimeFormat = errors.New("invalid time format")

	// ErrInvalidMemoryFormat indicates memory format is invalid
	ErrInvalidMemoryFormat = errors.New("invalid memory format")
)

// ValidationError represents a validation error when specs exceed limits
type ValidationError struct {
	Field     string // Field that failed validation
	Requested int    // Requested value
	Limit     int    // Maximum allowed value
	Partition string // Partition where limit applies
}

func (e *ValidationError) Error() string {
	if e.Partition != "" {
		return fmt.Sprintf("%s: requested %d exceeds limit %d for partition %s",
			e.Field, e.Requested, e.Limit, e.Partition)
	}
	return fmt.Sprintf("%s: requested %d exceeds limit %d",
		e.Field, e.Requested, e.Limit)
}

// Is allows errors.Is to match ValidationError
func (e *ValidationError) Is(target error) bool {
	_, ok := target.(*ValidationError)
	return ok
}

// ParseError represents an error parsing scheduler directives
type ParseError struct {
	Scheduler string // Scheduler name (e.g., "SLURM", "PBS")
	Line      int    // Line number where error occurred
	Content   string // Line content
	Reason    string // Reason for parse failure
}

func (e *ParseError) Error() string {
	if e.Line > 0 {
		return fmt.Sprintf("%s parse error at line %d (%s): %s",
			e.Scheduler, e.Line, e.Content, e.Reason)
	}
	return fmt.Sprintf("%s parse error: %s", e.Scheduler, e.Reason)
}

// SubmissionError represents an error during job submission
type SubmissionError struct {
	Scheduler string // Scheduler name
	JobName   string // Job name
	Output    string // Scheduler output
	Err       error  // Underlying error
}

func (e *SubmissionError) Error() string {
	if e.Output != "" {
		return fmt.Sprintf("%s submission failed for job %s: %v\nOutput: %s",
			e.Scheduler, e.JobName, e.Err, e.Output)
	}
	return fmt.Sprintf("%s submission failed for job %s: %v",
		e.Scheduler, e.JobName, e.Err)
}

func (e *SubmissionError) Unwrap() error {
	return e.Err
}

// ClusterError represents an error querying cluster information
type ClusterError struct {
	Scheduler string // Scheduler name
	Operation string // Operation that failed (e.g., "query GPUs", "get limits")
	Err       error  // Underlying error
}

func (e *ClusterError) Error() string {
	return fmt.Sprintf("%s cluster error during %s: %v",
		e.Scheduler, e.Operation, e.Err)
}

func (e *ClusterError) Unwrap() error {
	return e.Err
}

// DependencyError represents an error with job dependencies
type DependencyError struct {
	JobName       string   // Job name
	MissingDepIDs []string // Missing dependency job IDs
}

func (e *DependencyError) Error() string {
	return fmt.Sprintf("job %s has missing dependencies: %v",
		e.JobName, e.MissingDepIDs)
}

// ScriptCreationError represents an error creating a batch script
type ScriptCreationError struct {
	JobName string // Job name
	Path    string // Script path
	Err     error  // Underlying error
}

func (e *ScriptCreationError) Error() string {
	return fmt.Sprintf("failed to create script for job %s at %s: %v",
		e.JobName, e.Path, e.Err)
}

func (e *ScriptCreationError) Unwrap() error {
	return e.Err
}

// Helper functions for creating errors

// NewParseError creates a new ParseError
func NewParseError(scheduler string, line int, content string, reason string) *ParseError {
	return &ParseError{
		Scheduler: scheduler,
		Line:      line,
		Content:   content,
		Reason:    reason,
	}
}

// NewSubmissionError creates a new SubmissionError
func NewSubmissionError(scheduler string, jobName string, output string, err error) *SubmissionError {
	return &SubmissionError{
		Scheduler: scheduler,
		JobName:   jobName,
		Output:    output,
		Err:       err,
	}
}

// NewClusterError creates a new ClusterError
func NewClusterError(scheduler string, operation string, err error) *ClusterError {
	return &ClusterError{
		Scheduler: scheduler,
		Operation: operation,
		Err:       err,
	}
}

// NewDependencyError creates a new DependencyError
func NewDependencyError(jobName string, missingDepIDs []string) *DependencyError {
	return &DependencyError{
		JobName:       jobName,
		MissingDepIDs: missingDepIDs,
	}
}

// NewScriptCreationError creates a new ScriptCreationError
func NewScriptCreationError(jobName string, path string, err error) *ScriptCreationError {
	return &ScriptCreationError{
		JobName: jobName,
		Path:    path,
		Err:     err,
	}
}

// IsValidationError checks if an error is a ValidationError
func IsValidationError(err error) bool {
	var ve *ValidationError
	return errors.As(err, &ve)
}

// IsParseError checks if an error is a ParseError
func IsParseError(err error) bool {
	var pe *ParseError
	return errors.As(err, &pe)
}

// IsSubmissionError checks if an error is a SubmissionError
func IsSubmissionError(err error) bool {
	var se *SubmissionError
	return errors.As(err, &se)
}

// IsClusterError checks if an error is a ClusterError
func IsClusterError(err error) bool {
	var ce *ClusterError
	return errors.As(err, &ce)
}

// GpuValidationError represents an error when requested GPU is not available
type GpuValidationError struct {
	RequestedType  string    // Requested GPU type
	RequestedCount int       // Requested GPU count
	Suggestions    []string  // Suggested alternative GPUs
	AvailableGpus  []GpuInfo // All available GPUs on cluster
}

func (e *GpuValidationError) Error() string {
	var msg strings.Builder
	msg.WriteString(fmt.Sprintf("requested GPU %s:%d is not available",
		e.RequestedType, e.RequestedCount))

	if len(e.Suggestions) > 0 {
		msg.WriteString("\n\nSuggested alternatives:")
		for i, suggestion := range e.Suggestions {
			msg.WriteString(fmt.Sprintf("\n  %d. %s", i+1, suggestion))
		}
	} else if len(e.AvailableGpus) > 0 {
		msg.WriteString("\n\nAvailable GPUs:")
		for _, gpu := range e.AvailableGpus {
			msg.WriteString(fmt.Sprintf("\n  - %s: %d/%d available (partition: %s)",
				gpu.Type, gpu.Available, gpu.Total, gpu.Partition))
		}
	} else {
		msg.WriteString("\n\nNo GPUs are currently available on the cluster")
	}

	return msg.String()
}

// NewGpuValidationError creates a new GpuValidationError
func NewGpuValidationError(requestedType string, requestedCount int, suggestions []string, availableGpus []GpuInfo) *GpuValidationError {
	return &GpuValidationError{
		RequestedType:  requestedType,
		RequestedCount: requestedCount,
		Suggestions:    suggestions,
		AvailableGpus:  availableGpus,
	}
}

// IsGpuValidationError checks if an error is a GpuValidationError
func IsGpuValidationError(err error) bool {
	var gve *GpuValidationError
	return errors.As(err, &gve)
}
