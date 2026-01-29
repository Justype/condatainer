#!/usr/bin/env bash

SRC_DIR="cellranger-9.0.1"
PROCESSORS=4

# Output CSV
OUT="squash_bench.csv"
echo "algo,level,time_sec,size_mb" > "$OUT"

run_test() {
    local algo="$1"
    local level="$2"
    local extra=""

    local out_file="test-${algo}"
    [[ -n "$level" ]] && out_file="${out_file}-lv${level}"
    out_file="${out_file}.sqf"

    # compression options
    if [[ "$algo" == "gzip" ]]; then
        extra="-comp gzip -Xcompression-level $level"
    elif [[ "$algo" == "zstd" ]]; then
        extra="-comp zstd -Xcompression-level $level"
    else
        echo "Unknown algo: $algo"
        exit 1
    fi

    echo ">>> Running $algo level $level"

    local start=$(date +%s)

    # Perform compression, silence output
    mksquashfs "$SRC_DIR" "$out_file" \
        -processors "$PROCESSORS" \
        -keep-as-directory \
        $extra \
        >/dev/null 2>&1

    local end=$(date +%s)
    local elapsed=$((end - start))

    # File size in MB (no bc needed)
    local size_bytes=$(stat -c%s "$out_file")
    local size_mb=$((size_bytes / 1024 / 1024))

    echo "$algo,$level,$elapsed,$size_mb" >> "$OUT"
}

# --------------------------
# Test gzip (levels 1–9)
# --------------------------
for lv in {1..9}; do
    run_test gzip "$lv"
done

# --------------------------
# Test zstd (levels 1-22)
# --------------------------
for lv in {1..22}; do
    run_test zstd "$lv"
done

echo "Done. Results saved in $OUT"
