# CondaTainer: single-file environments and reference data for reproducible analysis on HPC

> Draft for *Bioinformatics* **Application Note** (≤4 pages; ~2,000 words plus one figure).
> Items marked `[TODO]` need a measurement or an author detail before submission.
> §6 describes freezing a writable environment as future work. Revisit once `overlay freeze`
> lands (see `plan/freeze-writable-env.md`).

**Authors.** [TODO: author list] <sup>1</sup>

<sup>1</sup> [TODO: affiliations]

<sup>\*</sup>To whom correspondence should be addressed.

---

## Abstract

**Summary:** A reproducible analysis needs three things preserved together: its code, its software
environment, and the reference data it consumed. Version control solves the first, and Conda and
containers largely solve the second, but neither gives any control over the third, so the genome
build, annotation release, and derived indexes a pipeline depends on go unrecorded. CondaTainer treats
software and reference data alike, packing each tool, environment, or dataset into one file that is
mounted read-only at run time and seen as an ordinary directory, and distributing them through
standard OCI container registries. A project records which of them its scripts use in a small file
committed alongside the code, and one command restores them on a machine that has never run
CondaTainer. Downstream analysis has different needs, since code development is open-ended, so
CondaTainer also provides a writable single-file environment that runs code editors on a compute node.

**Availability and implementation:** CondaTainer is free and open-source under the BSD 3-Clause
licence at https://github.com/Justype/condatainer, with documentation at
https://condatainer.readthedocs.io. It is a single statically linked Go binary for Linux x86-64,
installs without administrative privileges, and requires only Apptainer or Singularity on the cluster.
Recipes and helper scripts are at https://github.com/Justype/cnt-scripts.

**Contact:** [TODO: corresponding e-mail]

**Supplementary information:** Supplementary data are available at *Bioinformatics* online.

---

## 1 Introduction

An analysis is reproducible when someone else can obtain the same result from the same inputs. That
requires three things to survive together: the code, the software environment it ran in, and the
reference data it consumed. Version control addresses the first. The other two receive far less
attention, although a change of aligner patch version or of GENCODE release can alter a result as
readily as a change of code, and neither is recorded by a commit.

Conda, with the Bioconda channel (Grüning *et al.*, 2018), is how most bioinformatics software is
installed. An environment unpacks into tens of thousands of small files, which exhausts a user's inode
quota and loads a shared parallel filesystem such as Lustre with metadata operations, so creating,
copying, or removing one is slow. Containers, distributed per tool by BioContainers (da Veiga Leprevost *et al.*, 2017), avoid
that but fix their contents when they are built. Neither manages reference data: the genome build, the
annotation release, and the indexes derived from them belong to no channel, and neither does licensed
software such as 10x Genomics Cell Ranger. All of it is installed by hand at a cluster path and
referred to by that path, which leaves the components most likely to alter a result outside any
manifest.

Here we present CondaTainer, which manages software and reference data in the way code is already
managed: as independent files, each identified by content, recorded per project, and carried in
version control alongside the analysis that used them.

## 2 Architecture

CondaTainer runs containers with Apptainer (Kurtzer *et al.*, 2017), which needs neither privileges
nor a daemon. A *base image* is one file holding a minimal operating system, and it provides the
container root. Everything else is an *overlay*: one SquashFS file holding one tool, one Conda
environment, one operating-system layer, or one reference dataset, mounted read-only onto the base
when a command runs.

An analysis script declares the overlays it needs, alongside the scheduler directives it already
carries:

```bash
#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=42G
#SBATCH --time=2:00:00

#DEP: star/2.7.11b
#DEP: grch38/star/2.7.11b/gencode49-101

STAR --runThreadN $NCPUS --genomeDir $STAR_INDEX_DIR --readFilesIn r1.fq r2.fq
```

```bash
condatainer run align.sh
```

`condatainer run` reads both kinds of declaration, mounts the named overlays on the base, and submits
the job. The first overlay (`star/2.7.11b`) puts `STAR` on `PATH`. The second (`grch38/star/2.7.11b/gencode49-101`) declares `STAR_INDEX_DIR`, so the script names its reference data instead of hard-coding a cluster path, 
and the same script finds that data on any machine where the overlay is installed. The scheduler directives 
are read from the same script, so the job is submitted with the resources it asks for and `$NCPUS` is set 
from them. Array jobs and job chaining work the same way.

Tool, environment, and dataset overlays each occupy their own location inside the container, so any
combination can be mounted together and none has to be rebuilt when another changes. An
operating-system overlay contributes system paths and stacks on the base.

Packing an overlay as SquashFS has two effects. The overlay is a single file, so an environment that
would consume tens of thousands of inodes consumes one, and it is compressed, so it occupies less
space than the directory it replaces. Overlays are also installed into shared directories rather than
per user, so one copy of a GENCODE index serves everyone in a group.

**Figure 1** summarizes the design.

## 3 Building and distribution

An overlay is either built locally or downloaded, and the two routes are equivalent from the point of
view of an analysis that uses it.

Building starts from a *recipe*: a plain bash script whose metadata is written as comments, covering
build dependencies, the variables the overlay will declare, scheduler resources, and licence and
redistribution terms. Recipes are grouped into collections that resolve in order, so a lab collection
extends or shadows a public one, and each recipe is stored inside the overlay it produced.

Recipes also cover software that no channel may distribute. Cell Ranger is packaged by a short script
that fetches the vendor archive and unpacks it. The result is named, versioned, and used like any
other overlay, and its declared terms prevent it from being published to a public registry.

A Conda environment needs no recipe. `condatainer create -f environment.yml` solves the environment on
node-local storage, which keeps the metadata traffic off the shared filesystem, packs the result, and
stores both the requested environment and the exact package set it resolved to inside the overlay.

Every overlay records how it was built, and an identity is derived from that record rather than from
the payload bytes, because compressed output is not identical across tool versions. The identity names
one exact build.

Distribution uses the same registries that carry container images, such as GitHub Container Registry,
Quay, Harbor, or an institutional Artifactory. An overlay name becomes a repository and version, so a
published collection is browsable like any other set of container packages.

```bash
condatainer registry push grch38/star/2.7.11b/gencode49-101
```

The identity travels with the published overlay, so a consumer can confirm it is the intended build
before transferring tens of gigabytes, and a transfer that large is split into parts that resume after
an interruption. Fetching requires nothing beyond the CondaTainer binary, so a compute node needs no
additional container tooling and no administrative privileges. When a requested overlay is not
installed, CondaTainer queries the configured registries before building it locally, so a reference
index that takes hours to construct is built once and downloaded thereafter.

## 4 Recording an analysis

A project is a directory containing a record that is committed alongside the code.

```bash
condatainer project lock
```

This scans the analysis scripts for the dependencies they declare and resolves each to one exact
overlay. The record holds those identities, the recipes that rebuild each overlay, and the registry
locations each can be fetched from. It contains no payload and no machine-specific path, so it remains
small enough to review as a diff while describing the environment completely.

Two commands then use it. `condatainer project validate` recomputes every identity from the record
alone, with no network and no installed overlay, so a fresh clone can be checked before anything is
downloaded. `condatainer project restore` makes the environment available: overlays already installed
are reused, overlays with a recorded location are downloaded, and the rest are rebuilt from the
stored recipes. Downloads are addressed by content rather than by tag, so what arrives is the build
the project recorded.

`condatainer project push` publishes a project's own overlays to a registry that the project names.
For a published analysis, this keeps the record usable after the collection its software came from has
been reorganized or retired.

## 5 Downstream analysis

Pipelines account for only part of an analysis. Exploration, method development, and figure generation
are open-ended: packages are added as the work proceeds, and the environment is not known in advance.
A read-only overlay does not suit this stage, and large single-cell or imaging datasets require the
work to run on cluster hardware rather than on a personal computer.

CondaTainer provides a writable single-file environment for this stage. Packages are installed into it
as usual, and the state of the work is one file rather than an accumulation of installs across a home
directory.

```bash
condatainer overlay create -s 20G env.img
condatainer exec -w -o env.img bash
```

Helper scripts launch RStudio Server, VS Code, or a graphical desktop on a compute node through the
local scheduler and open a tunnel to it, so interactive work runs on the cluster with the same
overlays the pipeline used.

```bash
condatainer helper rstudio
```

## 6 Future directions

The writable environment is the one part of the workflow that is not yet recordable, because it has no
fixed content while it is being used. Its exact package set can be exported and rebuilt as an
immutable overlay at present, which covers environments assembled from Conda packages. Work in
progress packs the environment directly, so that software installed by other means, including system
packages added under `--fakeroot`, is also captured, allowing a settled interactive environment to be
recorded and shared like any other overlay.

Two limitations apply. A project record fetches by digest, so a downloaded overlay is verified to be
the published bytes, but published overlays are not signed and the registry's own access control is
what establishes who published them. SLURM is the tested scheduler, and PBS Pro, LSF, and HTCondor are
supported experimentally.

## Figure

**Figure 1.** *(A)* Software and reference data are packed as separate single-file overlays, mounted
together on a base image at run time, distributed through an OCI registry, and recorded per project in
a file committed with the code. *(B)* File count, creation time, and on-disk size for representative
Bioconda environments as ordinary Conda prefixes and as CondaTainer overlays `[TODO: measure, 3 to 4
environments]`. *(C)* Time to make a locked RNA-seq project runnable on a machine that had never run
CondaTainer, downloading recorded overlays against rebuilding them from recipes `[TODO: measure]`.
Compression settings and their trade-offs are in Supplementary Table [TODO].

## Acknowledgements

This project used computational resources provided by McMaster University and the Digital Research
Alliance of Canada. [TODO: further acknowledgements]

## Funding

[TODO: funding statement]

*Conflict of Interest:* none declared.

## References

- da Veiga Leprevost, F. *et al.* (2017) BioContainers: an open-source and community-driven framework
  for software standardization. *Bioinformatics*, **33**, 2580–2582.
- Grüning, B. *et al.* (2018) Bioconda: sustainable and comprehensive software distribution for the
  life sciences. *Nat. Methods*, **15**, 475–476.
- Kurtzer, G.M. *et al.* (2017) Singularity: scientific containers for mobility of compute.
  *PLoS ONE*, **12**, e0177459.
- Open Container Initiative (2021) OCI Distribution Specification v1.0.
  https://github.com/opencontainers/distribution-spec
- [TODO: Micromamba/mamba citation; optionally a workflow-manager citation (Nextflow/Snakemake) to
  position the project record against pipeline-level reproducibility]
