#!/usr/bin/env python3
# From https://packagemanager.posit.co/client/#/repos/cran/setup
# From https://packagemanager.posit.co/client/#/repos/bioconductor/setup

# Packages to comment out (can include simple glob suffix '*')
FILTERS = [
    "chromium", "dcraw", "gsfonts", "imagej", "nvidia-cuda-dev",
    "ocl-icd-opencl-dev", "pari-gp", "patch", "saga",
    "software-properties-common", "tcl", "tk", "tk-dev", "tk-table",
    "tesseract-ocr-eng", "texlive", "coinor-*",
]

import sys
import re
from pathlib import Path
from collections import defaultdict

if len(sys.argv) != 2:
    sys.exit("Usage: group_apt_installs.py <input.sh>")

input_path = Path(sys.argv[1]).resolve()
if not input_path.exists():
    sys.exit(f"File not found: {input_path}")

input_base_name = str(input_path.stem)
output_path = input_path.with_suffix(".group.sh")

# system package -> set of R packages
pkg_to_rpkgs = defaultdict(set)

current_rpkg = None

header_re = re.compile(r"^#\s*([A-Za-z0-9_.-]+)\s+requirements:")
apt_re = re.compile(r"^apt-get\s+install\s+-y\s+(.+)$")

with input_path.open() as f:
    for line in f:
        line = line.strip()

        if not line:
            continue

        # Section header
        m = header_re.match(line)
        if m:
            current_rpkg = m.group(1)
            continue

        if current_rpkg is None:
            continue

        # Ignore non-system commands
        if line.startswith("R CMD javareconf"):
            continue
        if line.startswith("curl "):
            continue

        # apt-get install parsing
        m = apt_re.match(line)
        if m:
            for pkg in m.group(1).split():
                pkg_to_rpkgs[pkg].add(current_rpkg)

with output_path.open("w") as out:
    out.write(f"    # Generated using {input_base_name} by group_apt_by_rpkg.py\n")
    # Emit non-lib packages first, then lib* packages last
    non_lib = [p for p in pkg_to_rpkgs.keys() if not p.startswith("lib")]
    lib_pkgs = [p for p in pkg_to_rpkgs.keys() if p.startswith("lib")]

    def is_filtered(pkg_name: str) -> bool:
        for pat in FILTERS:
            if pat.endswith("*"):
                if pkg_name.startswith(pat[:-1]):
                    return True
            else:
                if pkg_name == pat:
                    return True
        return False

    def write_pkg_line(pkg):
        rpkgs = " ".join(sorted(pkg_to_rpkgs[pkg]))
        comment_prefix = "# " if is_filtered(pkg) else ""
        out.write(f"    {comment_prefix}apt-get -y install {pkg} # {rpkgs}\n")

    for pkg in sorted(non_lib):
        write_pkg_line(pkg)
    for pkg in sorted(lib_pkgs):
        write_pkg_line(pkg)

print(f"Generated: {output_path}")
