//! Build system for creating overlays

pub mod conda;
pub mod def;
pub mod graph;
pub mod object;
pub mod shell;

pub use conda::CondaBuildObject;
pub use def::DefBuildObject;
pub use graph::BuildGraph;
pub use object::{BuildObject, BuildType};
pub use shell::ShellBuildObject;

use anyhow::Result;

/// Error types for build operations
#[derive(Debug, thiserror::Error)]
pub enum BuildError {
    #[error("Temporary overlay already exists: {0}")]
    TmpOverlayExists(String),
    
    #[error("Build cancelled by user")]
    BuildCancelled,
    
    #[error("Circular dependency detected involving '{0}'")]
    CircularDependency(String),
    
    #[error("Failed to create build object for '{0}': {1}")]
    BuildObjectCreation(String, String),
}

/// Get default number of CPUs for builds
pub fn default_build_ncpus() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
}
