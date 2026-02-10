## MEMO

### Apptainer

[Apptainer repo](https://github.com/apptainer/apptainer)

### Compression method for squashfs overlays

- Because the old apptainer version does not support zstd compression, `lz4` is used by default for better compatibility.
- Use `zstd` compression level 14 for best speed/size tradeoff for squashfs overlays.
  - (Use `--lz4` for faster reading speed, but much larger file size, about 2x)
  - `--gzip` is also supported, but not recommended.
- zstd has better compression ratio and IO speed than gzip.

I used cellranger 9.0.1 to benchmark various compression methods for squashfs.

- gzip 1-9 and zstd 1-22 were tested.
- folder size (2.4 GB)
- the results can be found in [squash_bench.csv](squash_bench.csv)

![squash_bench](squash_bench.png)
