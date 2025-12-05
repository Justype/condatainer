#!/usr/bin/env python3
"""
Scan `build-scripts/` and generate a JSON metadata file mapping
normalized spec strings to raw GitHub URLs for the script files.

- App packages: build-scripts/<name>/<version> -> key: "name=version" -> raw URL
- Reference packages: build-scripts/<assembly>/<data_type>/<version> -> key: "assembly=data_type=version" -> raw URL

The generated file is written to `build-scripts/metadata.json`.

This script expects the environment variables `GITHUB_REPOSITORY` and
`GITHUB_SHA` to be set (the workflow passes them). If they are not set,
it will still attempt to construct relative paths but the raw URLs will
be less useful.
"""
import os
import json

REPO = os.getcwd()
BUILD_SCRIPTS = os.path.join(REPO, 'build-scripts')
OUT_FILE = os.path.join(BUILD_SCRIPTS, 'metadata.json')

github_repo = os.environ.get('GITHUB_REPOSITORY')
github_sha = os.environ.get('GITHUB_SHA')

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
            # skip .sh .py scripts
            if data_type.endswith('.sh') or data_type.endswith('.py'):
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
                if github_repo and github_sha:
                    link = f"https://raw.githubusercontent.com/{github_repo}/{github_sha}/build-scripts/{script_name}/{data_type}/{version}"
                else:
                    link = os.path.relpath(script_path, REPO)

                # Extract metadata from the script file (e.g., #WHATIS, #URL)
                whatis = None
                website = None
                try:
                    with open(script_path, 'r') as sf:
                        for line in sf:
                            if line.startswith('#WHATIS:'):
                                whatis = line[len('#WHATIS:'):].strip()
                            elif line.startswith('#URL:'):
                                website = line[len('#URL:'):].strip()
                            # stop early if we got both
                            if whatis and website:
                                break
                except Exception:
                    pass

                metadata[key] = {
                    'name': script_name,
                    'link': link,
                    'whatis': whatis,
                    'website': website,
                }
    else:
        # structure is script_name/version
        for version in versions:
            if version.endswith('.sh') or version.endswith('.py'):
                continue
            script_path = os.path.join(script_relative_dir, version)
            if not os.path.isfile(script_path):
                continue
            key = f"{script_name}={version}"
            if github_repo and github_sha:
                link = f"https://raw.githubusercontent.com/{github_repo}/{github_sha}/build-scripts/{script_name}/{version}"
            else:
                link = os.path.relpath(script_path, REPO)

            whatis = None
            website = None
            try:
                with open(script_path, 'r') as sf:
                    for line in sf:
                        if line.startswith('#WHATIS:'):
                            whatis = line[len('#WHATIS:'):].strip()
                        elif line.startswith('#URL:'):
                            website = line[len('#URL:'):].strip()
                        if whatis and website:
                            break
            except Exception:
                pass

            metadata[key] = {
                'name': script_name,
                'link': link,
                'whatis': whatis,
                'website': website,
            }

# write metadata
os.makedirs(os.path.dirname(OUT_FILE), exist_ok=True)
with open(OUT_FILE, 'w') as f:
    json.dump(metadata, f, indent=2, sort_keys=True)

print(f'Wrote metadata for {len(metadata)} entries to {OUT_FILE}')
