# CondaTainer

Use [apptainer](https://apptainer.org)/[conda](https://anaconda.org)/[squashFS](https://github.com/plougher/squashfs-tools) to manage tools for HPC users.

On HPC, conda is not a good solution, because HPC file system is designed for large files, but conda creates many small files. This app uses [micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html) to create a conda environment within an overlay which can be mounted to apptainer/singularity container.

## Setup

Before using CondaTainer, make sure you have apptainer/singularity and squashfs-tools installed and working. If not, please contact your HPC admin.

Do not clone this repo. Just download the [bin/condatainer](bin/condatainer) file and put it in your PATH.

```bash
mkdir -p $SCRATCH/condatainer/bin
wget -O $SCRATCH/condatainer/bin/condatainer 'https://raw.githubusercontent.com/Justype/condatainer/refs/heads/main/bin/condatainer'
chmod +x $SCRATCH/condatainer/bin/condatainer
# Append to PATH
if ! grep -q 'condatainer/bin' ~/.bashrc; then
    echo 'export PATH=$SCRATCH/condatainer/bin:$PATH' >> ~/.bashrc
fi
source ~/.bashrc
```

## CondaTainer command

```bash
# Will create trim-galore, cutadapt, fastqc overlays
condatainer create trim-galore=0.6.10 cutadapt=5.2 fastqc=0.12.1
condatainer create grch38=salmon-1.10.3=gencode44 # create genome index overlay


# Create one overlay with all packages
condatainer create -n sci-rna-seq trim-galore=0.6.10 star=2.7.11b ...
condatainer create -n sci-rna-seq -f sci-rna-seq-packages.yml
condatainer create -p ./sci_rna -f sci-rna-seq-packages.yml # p = prefix
# will create sci_rna.sqf under current folder

condatainer avail # check available local build scripts (please use grep to filter)
condatainer avail grch star 47 # search available build scripts by keywords
condatainer list  # list created overlays

# Run command with specified overlays
condatainer exec -o trim-galore=0.6.10 -o cutadapt=5.2 -o fastqc=0.12.1 command

condatainer check script.sh # check if all required overlays are installed
condatainer check script.sh --auto-install # auto install missing overlays

condatainer run script.sh # run a bash script with required overlays
```

For more help, go to [MANUAL.md](assets/MANUAL.md)

### script example

My script will try to use `#DEP:` comments to get required overlays.

You can use it with SLURM job script as well. If is in a container, `IN_CONDATINER` will be set to `1`.

```bash
#!/bin/bash
#SBATCH --time=1-0:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32GB
##SBATCH --gpus=nvidia_h100_80gb_hbm3_1g.10gb:1
#SBATCH --job-name=script_test
#SBATCH --output=script_test_%j.out
##SBATCH --mail-type=END

#DEP:samtools/1.22
#DEP:bcftools/1.22
#DEP:htslib/1.22
#DEP:grcm39/gtf-gencode/M33
#DEP:./custom_image.img
#CNT --writable-img

if [ -z "$IN_CONDATINER" ]; then
    condatainer run "$0" "$@"
    exit $?
fi

# Do whatever you want with the required tools
samtools --version
bcftools --version
bgzip --version
ls $ANNOTATION_GTF
```

- When runing `condatainer check script.sh`, check if all required overlays are installed.
- When runing `condatainer run script.sh`, auto load required overlays.
- Ambiguous versions is not allowed, e.g., `#DEP:trim-galore` will raise error.

## Overlay Query Steps

1. Check `build-scripts/` folder.
2. Check this git repo [metadata/build-scripts.json](./metadata/build-scripts.json).
3. If not found, try to create from `conda-forge` and `bioconda` channel.

## Naming convention for CondaTainer overlays

extensions:

- `.sqf` suffix for squashfs overlay files (read-only)
- `.img` suffix for ext3 image files (read-write, can be created by `apptainer overlay create`)

SQF naming:

- `conda-env.sqf` is designed for storing a group of conda packages within one overlay.
  - overlay path: `/cnt/conda-env/`
- `package-name--version.sqf` is designed for storing a specific package.
  - overlay path: `/cnt/package-name/version/`
- `assembly-name--data-type--version.sqf` is designed for storing genome data, e.g., `grch38--salmon-index--gencode44.sqf`
  - overlay path: `/cnt/assembly-name/data-type/version/`

IMG naming:

- `name.img` (do not use `--` in the name)
  - overlay path: `/ext3/name/`
  - Please do not use `conda` as the name. My `micromamba` will use /ext3/conda/ for its root prefix.

> [!NOTE]
> Do not change the name of the overlay files after creation. Otherwise, the internal paths will not match. e.g. R libs may use symlinks.

## How does it work?

- The `app-version.sqf` overlay file contains a micromamba conda environment with the specified app installed.
- The abs path of the env is `/cnt/app/version/` within the overlay.
- When running the apptainer container,
  - Overlays are mounted.
  - The `PATH` is modified to include overlays' `bin/`.
  - The order of overlays matters. The later overlays have higher priority in PATH.

```bash
# Make sure you put cutadapt later
# trim-galore conda env will also have cutadapt installed
condatainer exec \
  -o trim-galore=0.6.10 \
  -o cutadapt=5.0 \
  bash -c "fastqc --version && cutadapt --version"

# Will equivalently run:
apptainer exec \
  --env PATH="/cnt/fastqc/0.12.1/bin:/cnt/cutadapt/5.0/bin:/sbin:/bin" \
  -o /path/to/trim-galore--0.6.10.sqf:ro \
  -o /path/to/cutadapt--5.0.sqf:ro \
  images/base_image.sif \
  bash -c "fastqc --version && cutadapt --version"
# Output:
# FastQC v0.12.1
# 5.0
```

Because `trim-galore` depends on `cutadapt`, if you do not specify `-o cutadapt=5.0`, the version installed with `trim-galore` will be used.

```bash
condatainer exec \
  -o trim-galore=0.6.10 \
  bash -c "fastqc --version && cutadapt --version"
# Output:
# FastQC v0.12.1
# 5.2
```

## Q&A

- Q: Can I use GPU within the container?
  - A: Yes. If `nvidia-smi` or `rocm-smi` command is available on the host system, the corresponding apptainer flag `--nv` or `--rocm` will be automatically added when running the container. (But you need to change the building source of base image def from `debian` to `nvidia/cuda` or similar)
- Q: Can I use conda channels other than `conda-forge` and `bioconda`?
  - A: You can use `-f file.yml` to specify a conda environment yaml file when creating overlays.
- Q: Can I modify the base image?
  - A: Yes. You need to place your custom def file at `condatainer/build-scripts/base_image.def` and delete the existing base image `condatainer/images/base_image.sif`. The next time you run `condatainer create`, the base image will be rebuilt.
- Q: Can I use mamba instead of micromamba?
  - A: No. Micromamba is a single binary file and easier to use within the overlay. Mamba requires a full conda installation.
- Q: `srun` or `sbatch`?
  - A: It does not matter only if the apps require interactive input (e.g. cellranger install script). If yes, use `srun --pty bash` to get an interactive shell.

## MEMO

### Apptainer Requirements

[Apptainer repo](https://github.com/apptainer/apptainer)

```bash
sudo apt install uidmap squashfuse gocryptfs fuse-overlayfs
```

### Compression method for squashfs overlays

- Because the old apptainer version does not support zstd compression, `lz4` is used by default for better compatibility.
- Use `zstd` compression level 14 for best speed/size tradeoff for squashfs overlays. 
  - (Use `--lz4` for faster reading speed, but much larger file size, about 2x)
  - `--gzip` is also supported, but not recommended.
- zstd has better compression ratio and IO speed than gzip.

I used cellranger 9.0.1 to benchmark various compression methods for squashfs.

- gzip 1-9 and zstd 1-22 were tested.
- folder size (2.4 GB)
- the results can be found in [squash_bench.csv](assets/squash_bench.csv)

![squash_bench](assets/squash_bench.png)

## References

- [NYU HPC Using Containers on HPC](https://services.rt.nyu.edu/docs/hpc/containers/containers/)
- [NYU HPC Singularity with Conda](https://services.rt.nyu.edu/docs/hpc/containers/singularity_with_conda/)
- [NYU HPC SquashFS and Singularity](https://services.rt.nyu.edu/docs/hpc/containers/squash_file_system_and_singularity/)
- [Compression Method Benchmarks](https://github.com/inikep/lzbench)
- [AllianceCAN Multi-Instance GPU](https://docs.alliancecan.ca/wiki/Multi-Instance_GPU)
- [AllianceCAN CPU RAM GPU ratio](https://docs.alliancecan.ca/wiki/Allocations_and_compute_scheduling#Ratios_in_bundles)
