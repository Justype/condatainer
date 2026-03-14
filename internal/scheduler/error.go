package scheduler

import (
	"errors"
	"fmt"
	"time"
)

// Common errors
var (
	// ErrSchedulerNotFound indicates the scheduler binary was not found
	ErrSchedulerNotFound = errors.New("scheduler binary not found in PATH")

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

// NewScriptCreationError creates a new ScriptCreationError
func NewScriptCreationError(jobName string, path string, err error) *ScriptCreationError {
	return &ScriptCreationError{
		JobName: jobName,
		Path:    path,
		Err:     err,
	}
}

// TimeoutError represents a scheduler command that exceeded the configured timeout
type TimeoutError struct {
	Scheduler string
	Operation string
	Timeout   time.Duration
}

func (e *TimeoutError) Error() string {
	return fmt.Sprintf("%s: %q timed out after %s", e.Scheduler, e.Operation, e.Timeout)
}

