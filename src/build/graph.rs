//! Build graph for dependency resolution and build ordering

use anyhow::Result;
use std::collections::{HashMap, HashSet};
use std::path::Path;

use super::object::BuildObject;
use super::BuildError;
use crate::utils;

/// BuildGraph manages dependency resolution and build ordering
/// 
/// It expands BuildObject items to include all transitive dependencies,
/// detects cycles, produces a topologically-sorted order (dependencies before dependents),
/// and separates items into two ordered lists:
/// - local_builds: BuildObjects to build locally
/// - scheduler_builds: BuildObjects to submit via scheduler
pub struct BuildGraph {
    graph: HashMap<String, Box<dyn BuildObject>>,
    local_builds: Vec<String>,
    scheduler_builds: Vec<String>,
    job_ids: HashMap<String, String>,
    
    images_dir: String,
    tmp_dir: String,
    submit_jobs: bool,
}

impl BuildGraph {
    /// Create a new BuildGraph from a list of BuildObjects
    pub fn new(
        build_objects: Vec<Box<dyn BuildObject>>,
        images_dir: impl AsRef<Path>,
        tmp_dir: impl AsRef<Path>,
        submit_jobs: bool,
    ) -> Result<Self> {
        let mut graph = HashMap::new();
        
        // Seed graph with provided build objects
        for obj in build_objects {
            let name = obj.name_version().to_string();
            if obj.is_installed() {
                utils::print_message(&format!(
                    "Overlay {} is already installed. Skipping.",
                    utils::style_name(&name)
                ));
            }
            graph.insert(name, obj);
        }
        
        let mut bg = Self {
            graph,
            local_builds: Vec::new(),
            scheduler_builds: Vec::new(),
            job_ids: HashMap::new(),
            images_dir: images_dir.as_ref().to_string_lossy().to_string(),
            tmp_dir: tmp_dir.as_ref().to_string_lossy().to_string(),
            submit_jobs,
        };
        
        // Topologically sort the graph
        bg.topological_sort()?;
        
        Ok(bg)
    }
    
    /// Perform topological sort with cycle detection
    /// Separates builds into local and scheduler lists
    fn topological_sort(&mut self) -> Result<()> {
        let mut visiting = HashSet::new();
        let mut visited = HashSet::new();
        let mut order = Vec::new();
        
        // Get all nodes to visit
        let nodes: Vec<String> = self.graph.keys().cloned().collect();
        
        for node in nodes {
            self.visit(&node, &mut visiting, &mut visited, &mut order)?;
        }
        
        // Separate into local and scheduler builds
        for name in order {
            if let Some(obj) = self.graph.get(&name) {
                if obj.requires_scheduler() && self.submit_jobs {
                    self.scheduler_builds.push(name);
                } else {
                    self.local_builds.push(name);
                }
            }
        }
        
        Ok(())
    }
    
    /// Visit a node in the dependency graph (recursive DFS)
    fn visit(
        &mut self,
        node: &str,
        visiting: &mut HashSet<String>,
        visited: &mut HashSet<String>,
        order: &mut Vec<String>,
    ) -> Result<()> {
        if visited.contains(node) {
            return Ok(());
        }
        
        if visiting.contains(node) {
            return Err(BuildError::CircularDependency(node.to_string()).into());
        }
        
        visiting.insert(node.to_string());
        
        // Get or create the build object
        let deps = if let Some(obj) = self.graph.get(node) {
            obj.dependencies().to_vec()
        } else {
            // Try to create build object for dependency
            // TODO: Implement NewBuildObject factory function
            utils::print_warning(&format!("Dependency '{}' not found in graph", node));
            Vec::new()
        };
        
        // Visit dependencies first
        for dep in deps {
            self.visit(&dep, visiting, visited, order)?;
        }
        
        visiting.remove(node);
        visited.insert(node.to_string());
        order.push(node.to_string());
        
        Ok(())
    }
    
    /// Run the build graph (build all objects in order)
    pub fn run(&mut self) -> Result<()> {
        utils::print_message("Starting build process...");
        
        // Build local builds first
        for name in &self.local_builds.clone() {
            if let Some(obj) = self.graph.get(name) {
                if !obj.is_installed() {
                    utils::print_message(&format!(
                        "Building: {}",
                        utils::style_name(name)
                    ));
                    obj.build(true)?;
                }
            }
        }
        
        // Then scheduler builds
        if !self.scheduler_builds.is_empty() {
            utils::print_message(&format!(
                "Submitting {} jobs to scheduler",
                self.scheduler_builds.len()
            ));
            
            for name in &self.scheduler_builds.clone() {
                if let Some(obj) = self.graph.get(name) {
                    if !obj.is_installed() {
                        // TODO: Implement scheduler job submission
                        utils::print_warning(&format!(
                            "Scheduler submission not yet implemented for: {}",
                            name
                        ));
                    }
                }
            }
        }
        
        utils::print_message("Build process completed!");
        Ok(())
    }
    
    /// Get the build order (for debugging/testing)
    pub fn build_order(&self) -> Vec<String> {
        let mut order = Vec::new();
        order.extend(self.local_builds.iter().cloned());
        order.extend(self.scheduler_builds.iter().cloned());
        order
    }
    
    /// Get number of local builds
    pub fn local_build_count(&self) -> usize {
        self.local_builds.len()
    }
    
    /// Get number of scheduler builds
    pub fn scheduler_build_count(&self) -> usize {
        self.scheduler_builds.len()
    }
}
