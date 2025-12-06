#!/usr/bin/env python3
"""
Scan `build-scripts/` and generate a JSON metadata file mapping
normalized spec strings to raw GitHub URLs for the script files.

- App packages: build-scripts/<name>/<version> -> key: "name=version" -> raw URL
- Reference packages: build-scripts/<assembly>/<data_type>/<version> -> key: "assembly=data_type=version" -> raw URL

The generated file is written to `metadata/build-scripts.json`.
"""
import os
import json

REPO = os.getcwd()
BUILD_SCRIPTS = os.path.join(REPO, 'build-scripts')

github_repo = os.environ.get('GITHUB_REPOSITORY')
github_sha = os.environ.get('GITHUB_SHA')
OUT_DIR = os.path.join(REPO, 'metadata')
OUT_FILE = os.path.join(OUT_DIR, 'build-scripts.json')

metadata = {}

if not os.path.isdir(BUILD_SCRIPTS):
    print('No build-scripts directory found; writing empty metadata file.')
    with open(OUT_FILE, 'w') as f:
        json.dump({}, f, indent=2, sort_keys=True)
    raise SystemExit(0)

for script_name in sorted(os.listdir(BUILD_SCRIPTS)):
    if script_name == '0-template':
        continue
    script_relative_dir = os.path.join(BUILD_SCRIPTS, script_name)
    if not os.path.isdir(script_relative_dir):
        continue

    # gather version entries (files or directories)
    versions = [v for v in sorted(os.listdir(script_relative_dir)) if not (v.startswith('template') or v.endswith('_data'))]
    if not versions:
        continue

    first = versions[0]
    # If the first entry is a directory, assume structure is script_name/data_type/version
    if os.path.isdir(os.path.join(script_relative_dir, first)):
        for data_type in versions:
            if data_type.startswith("template") or data_type.endswith('.sh') or data_type.endswith('.py'):
                continue
            data_type_path = os.path.join(script_relative_dir, data_type)
            if not os.path.isdir(data_type_path):
                continue
            sub_versions = sorted([sv for sv in os.listdir(data_type_path) if not sv.startswith('template')])
            for version in sub_versions:
                script_path = os.path.join(data_type_path, version)
                if not os.path.isfile(script_path):
                    continue
                key = f"{script_name}={data_type}={version}"
                link = os.path.relpath(script_path, REPO)

                # Extract metadata from the script file (e.g., #WHATIS, #URL)
                whatis = None
                try:
                    with open(script_path, 'r') as sf:
                        for line in sf:
                            if line.startswith('#WHATIS:'):
                                whatis = line[len('#WHATIS:'):].strip()
                            # stop early if we got both
                            if whatis:
                                break
                except Exception:
                    pass

                metadata[key] = whatis if whatis else "Missing description"
    else:
        # structure is script_name/version
        for version in versions:
            if version.endswith('.sh') or version.endswith('.py') or version.startswith("template"):
                continue
            script_path = os.path.join(script_relative_dir, version)
            if not os.path.isfile(script_path):
                continue
            key = f"{script_name}={version}"
            link = os.path.relpath(script_path, REPO)

            whatis = None
            try:
                with open(script_path, 'r') as sf:
                    for line in sf:
                        if line.startswith('#WHATIS:'):
                            whatis = line[len('#WHATIS:'):].strip()
                        if whatis:
                            break
            except Exception:
                pass

            metadata[key] = whatis if whatis else "Missing description"

# write metadata
os.makedirs(os.path.dirname(OUT_FILE), exist_ok=True)
with open(OUT_FILE, 'w') as f:
    json.dump(metadata, f, indent=2, sort_keys=True)

# write metadata to gzipped version
import gzip
with gzip.open(OUT_FILE + '.gz', 'wt') as f:
    json.dump(metadata, f, indent=2, sort_keys=True)

print(f'Wrote metadata for {len(metadata)} entries to {OUT_FILE}')
