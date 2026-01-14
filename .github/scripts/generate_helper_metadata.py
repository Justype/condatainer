#!/usr/bin/env python3
"""Generate helper metadata for scripts grouped by folders under `helpers/`.

Scans each top-level subfolder under `helpers/` (for example
`helpers/slurm`, `helpers/headless`, etc.) and writes JSON to
`metadata/helper-scripts.json` with structure:
{
    "slurm": {
        "script_name": { "path": "helpers/slurm/script_name" },
        ...
    },
    "headless": { ... }
}

Script is safe to run locally or in CI.
"""
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "metadata"
OUT_FILE = OUT_DIR / "helper-scripts.json"

HELPERS_ROOT = ROOT / "helpers"

# Dynamically discover all top-level folders under helpers/ and use the folder
# name as the metadata category key. This allows adding new groups without
# changing this script (e.g., 'slurm', 'headless', ...).
FOLDERS = {}
if HELPERS_ROOT.exists():
    for p in sorted(HELPERS_ROOT.iterdir()):
        if p.is_dir():
            FOLDERS[p.name] = p

EXT_WHITELIST = {".sh", ".py", ".bash", ".md"}


def scan_folder(folder_path: Path):
    out = {}
    if not folder_path.exists():
        return out
    for p in sorted(folder_path.rglob("*")):
        if p.is_file():
            if p.suffix and p.suffix.lower() not in EXT_WHITELIST:
                continue
            # Use repo-relative path
            rel = p.relative_to(ROOT).as_posix()
            key = p.stem
            # If duplicate keys (same stem) exist, make a unique key using path
            if key in out:
                key = rel
            out[key] = {"path": rel}
    return out


def main():
    data = {}
    for k, folder in FOLDERS.items():
        data[k] = scan_folder(folder)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    json_text = json.dumps(data, indent=2, ensure_ascii=False)
    with OUT_FILE.open("w", encoding="utf-8") as f:
        f.write(json_text)

    # Also write a gzipped copy to match generate-metadata.yml style
    try:
        import gzip
        with gzip.open(str(OUT_FILE) + ".gz", "wb") as gz:
            gz.write(json_text.encode("utf-8"))
    except Exception as e:
        print(f"Warning: failed to write gzipped metadata: {e}")

    print(f"Wrote helper metadata to {OUT_FILE} and {OUT_FILE}.gz")


if __name__ == "__main__":
    main()
