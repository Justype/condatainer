# Scheduler Implementation Comparison

## Resource Querying Functions

### SLURM (Most Complete - Reference Implementation)
| Function | Purpose | Status |
|----------|---------|--------|
| `getGpuInfo()` | Query GPU types per partition with availability | ✅ **Implemented** |
| `getPartitionLimits()` | Query partition config limits | ✅ **Implemented** |
| `getAvailableResourcesByPartition()` | Query actual available CPUs/memory per partition | ✅ **Implemented** |
| `getMaxNodeResources()` | Query max CPUs/memory per node | ✅ **Implemented** |

**Key Features:**
- GPU info includes partition assignment and available/total counts
- Merges configured limits with actual available resources
- Comprehensive per-partition resource tracking

---

### PBS (Experimental)
| Function | Purpose | Status |
|----------|---------|--------|
| `getGpuInfo()` | Query GPU info with optional queue assignment | ✅ **IMPLEMENTED** |
| `getNodeResources()` | Query node CPUs and memory | ✅ Implemented (refactored) |
| `getQueueLimits()` | Query queue resource limits via `qstat -Qf` | ✅ **Implemented** |
| `getAvailableResourcesByQueue()` | Query actual available CPUs/memory per queue | ✅ **IMPLEMENTED** |

**Recent Changes:**
- ✅ **Added `getGpuInfo()`** - Separate function for GPU detection
- ✅ **Refactored `getNodeResources()`** - Now returns only CPUs and memory (signature changed from returning 3 values to 2)
- ✅ **Queue assignment support** - GPUs include queue info if available from pbsnodes
- ✅ **Added `getAvailableResourcesByQueue()`** - Queries pbsnodes and aggregates resources by queue
- ✅ **Updated `getQueueLimits()`** - Now merges actual available resources into queue limits

---

### LSF (Experimental)  
| Function | Purpose | Status |
|----------|---------|--------|
| `getGpuInfo()` | Query GPU info from bhosts | ✅ **IMPLEMENTED** |
| `getHostResources()` | Query host CPUs and memory | ✅ **FIXED** - Now queries memory |
| `getQueueLimits()` | Query queue resource limits via `bqueues -l` | ✅ **Implemented** |
| `getAvailableResourcesByQueue()` | Query actual available CPUs/memory per queue | ✅ **IMPLEMENTED** |

**Recent Changes:**
- ✅ **Added `getGpuInfo()`** - Queries GPU info via `bhosts -gpu` (LSF 10.1+)
- ✅ **Added `getGpuInfoFromHostDetails()`** - Fallback for older LSF versions
- ✅ **Fixed `getHostResources()`** - Now queries memory from `lshosts`
- ✅ **Added `parseLsfMemoryToMB()`** - Parses LSF memory formats (32G, 65536M, etc.)
- ✅ **Added `getAvailableResourcesByQueue()`** - Queries bhosts/lshosts and aggregates by queue
- ✅ **Updated `getQueueLimits()`** - Now merges actual available resources into queue limits

---

### HTCondor (Experimental)
| Function | Purpose | Status |
|----------|---------|--------|
| `getGpuInfo()` | Query GPU info from condor_status | ✅ **Implemented** |
| `getMaxNodeResources()` | Query max CPUs/memory per node | ✅ **Implemented** |
| `getClusterLimits()` | Query cluster-wide limits (no partitions) | ✅ **Implemented** |
| `getConfigValue()` | Query HTCondor configuration | ✅ **Implemented** |
| Equivalent to `getAvailableResourcesByPartition()` | **N/A** - HTCondor doesn't use partitions | N/A |

**Notes:**
- HTCondor uses matchmaking, not partitions/queues
- Returns single "default" limit representing cluster-wide resources
- Architecture is fundamentally different from SLURM/PBS/LSF

---

## Summary of Missing Functions

### All Critical Features Implemented! ✅

As of Feb 16, 2026, all schedulers now have feature parity with SLURM for resource querying:

#### PBS ✅
- ✅ **`getGpuInfo()`** with queue assignment - Implemented
- ✅ **`getAvailableResourcesByQueue()`** - Implemented
- ✅ **`getQueueLimits()`** merges actual available resources - Implemented

#### LSF ✅
- ✅ **`getGpuInfo()`** - Complete GPU support implemented
- ✅ **`getHostResources()`** with memory querying - Fixed
- ✅ **`getAvailableResourcesByQueue()`** - Implemented
- ✅ **`getQueueLimits()`** merges actual available resources - Implemented

#### HTCondor ✅
- ✅ All applicable features implemented (uses cluster-wide model, not partitions)

### Comparison with SLURM Features

| Feature | SLURM | PBS | LSF | HTCondor |
|---------|-------|-----|-----|----------|
| GPU detection | ✅ Per-partition | ✅ With queue info | ✅ From bhosts | ✅ Per-machine |
| GPU availability tracking | ✅ Available/Total | ✅ Available/Total | ✅ Total (available=total) | ✅ Total only |
| Partition/Queue limits | ✅ Full | ✅ Full | ✅ Full | ✅ Cluster-wide |
| Actual available resources | ✅ Per-partition | ✅ Per-queue | ✅ Per-queue | ✅ Cluster-wide |
| Memory querying | ✅ Per-node | ✅ Per-node | ✅ Per-node | ✅ Per-node |

---

## Implementation Priority

### ✅ All Features Completed (Feb 16, 2026)

**Phase 1: Basic GPU and Memory Support**
1. ✅ **LSF: Add GPU support** - `getGpuInfo()` using `bhosts -gpu` with fallback
2. ✅ **LSF: Fix memory querying** - Updated `getHostResources()` to query memory from `lshosts`
3. ✅ **PBS: Add `getGpuInfo()`** - Extracted GPU info with optional queue assignment from pbsnodes
4. ✅ **PBS: Refactored `getNodeResources()`** - Simplified to return only CPUs and memory

**Phase 2: Available Resources Per Queue (Completed Today)**
5. ✅ **PBS: Add `getAvailableResourcesByQueue()`** - Queries actual per-queue resources from pbsnodes
6. ✅ **LSF: Add `getAvailableResourcesByQueue()`** - Queries actual per-queue resources from bhosts/lshosts
7. ✅ **PBS: Updated `getQueueLimits()`** - Merges available resources and GPU info by queue
8. ✅ **LSF: Updated `getQueueLimits()`** - Merges available resources and GPU info by queue

### Future Enhancements (Optional)
- ⚠️ **HTCondor: Accounting group limits** - Add accounting group-specific resource limits (low priority)
- ⚠️ **LSF: Improved queue-host mapping** - More accurate queue-to-host mapping via LSF API (low priority)
- ⚠️ **PBS: Enhanced GPU type detection** - Parse gpu_model field if available (low priority)
