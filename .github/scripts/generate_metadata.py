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
BUILD_SCRIPTS_DIR = os.path.join(REPO, 'build-scripts')

github_repo = os.environ.get('GITHUB_REPOSITORY')
github_sha = os.environ.get('GITHUB_SHA')
OUT_DIR = os.path.join(REPO, 'metadata')
OUT_FILE = os.path.join(OUT_DIR, 'build-scripts.json')

def get_local_build_scripts():
    """
    Get a name-version to local script path mapping dictionary.
    """
    packages = {}
    if not os.path.isdir(BUILD_SCRIPTS_DIR):
        return packages

    # os.walk mimics 'find' by visiting every subdirectory recursively
    for root, _, files in os.walk(BUILD_SCRIPTS_DIR):
        for filename in files:
            # 1. skip non-build-script files
            if filename.endswith(('.py', '.sh')):
                continue

            full_path = os.path.join(root, filename)

            # 2. skip template files
            if "template" in full_path:
                continue

            # 3. generate key (relative path like 'apps/tool/v1')
            relative_key = os.path.relpath(full_path, BUILD_SCRIPTS_DIR)

            if relative_key.endswith('.def'):
                if relative_key.startswith(('base_image', 'base-overlay')):
                    continue  # skip base scripts
                relative_key = relative_key[:-4]  # remove .def suffix
            
            packages[relative_key] = full_path

    return packages



metadata = {}

if not os.path.isdir(BUILD_SCRIPTS_DIR):
    print('No build-scripts directory found; writing empty metadata file.')

for script_name, path in get_local_build_scripts().items():
    script_metadata = {}
    script_metadata['relative_path'] = os.path.relpath(path, BUILD_SCRIPTS_DIR)
    script_metadata['whatis'] = "Missing description"
    script_metadata['deps'] = []
    script_metadata['sbatch'] = False
    try:
        with open(path, 'r') as sf:
            for line in sf:
                if line.startswith('#WHATIS:'):
                    script_metadata['whatis'] = line[len('#WHATIS:'):].strip()
                if line.startswith('#DEP:'):
                    dep = line[len('#DEP:'):].strip()
                    script_metadata['deps'].append(dep)
                if line.startswith('#SBATCH'):
                    script_metadata['sbatch'] = True
    except Exception:
        pass

    metadata[script_name] = script_metadata

# write metadata
os.makedirs(os.path.dirname(OUT_FILE), exist_ok=True)
with open(OUT_FILE, 'w') as f:
    json.dump(metadata, f, indent=2, sort_keys=True)

# write metadata to gzipped version
import gzip
with gzip.open(OUT_FILE + '.gz', 'wt') as f:
    json.dump(metadata, f, indent=2, sort_keys=True)

print(f'Wrote metadata for {len(metadata)} entries to {OUT_FILE}')
