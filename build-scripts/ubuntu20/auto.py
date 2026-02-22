#!/usr/bin/env python3
"""
Query the Docker Hub to generate ubuntu20/*def files

- R related => query posit/r-base tags
    - For we look for tags like 4.4.3-focal-{arch} => 4.4.3-focal in r4.4.3.def
    python3 build-scripts/ubuntu20/auto.py # to generate new files
    rm build-scripts/ubuntu20/r*.def # to remove generated files
"""

import json
import os
import re
import sys
import urllib.request

SCRIPT_PATH = __file__
BASE_DIR = os.path.dirname(SCRIPT_PATH)

CODE_NAME = "focal"

class RVersions:
    """Discover r-base versions and create corresponding definition files."""

    # expose the code name at class level as well for easy access
    CODE_NAME = CODE_NAME

    # template for each generated definition file; CODE_NAME is inserted
    TPL = "#WHATIS: Ubuntu ${CODE_NAME} R ${VERSION}\n" + \
          "Bootstrap: docker\n" + \
          "From: posit/r-base:${VERSION}-${CODE_NAME}\n\n" + \
          "%environment\n" + \
          "    export LC_ALL=C.UTF-8\n" + \
          "    export LANG=C.UTF-8\n" + \
          "    export PS1='CNT \\w> '\n\n" + \
          "%runscript\n" + \
          "  exec /usr/bin/R \"$@\"\n\n"

    def __init__(self):
        self.base_dir = BASE_DIR

    def query_docker_tags(self, url):
        """Return a list of tag names from the Docker Hub API.
        Follow the "next" pagination links until all tags have been retrieved or an error occurs.
        """
        tags = []
        while url:
            try:
                with urllib.request.urlopen(url, timeout=10) as resp:
                    data = json.load(resp)
            except Exception as e:
                print(f"Warning: failed to query {url}: {e}", file=sys.stderr)
                break

            for entry in data.get("results", []):
                name = entry.get("name")
                if name:
                    tags.append(name)

            url = data.get("next")
        return tags

    def filter_versions(self, tag_list, pattern=None):
        """Extract semver versions from relevant tags.

        By default, search for ``major.minor.patch-{CODE_NAME}`` with an 
            optional architecture suffix (e.g. ``-amd64`` or ``-arm64``).
        """

        if pattern is None:
            pattern = rf"^(\d+\.\d+\.\d+)-{self.CODE_NAME}(?:-[^-]+)?$"

        regex = re.compile(pattern)
        vers = set()
        for t in tag_list:
            m = regex.match(t)
            if m:
                vers.add(m.group(1))
        # sort by numeric components
        def keyfn(v):
            return tuple(int(x) for x in v.split("."))

        return sorted(vers, key=keyfn)

    def run(self):
        tags = self.query_docker_tags("https://hub.docker.com/v2/repositories/posit/r-base/tags?page_size=100")
        if not tags:
            print("No tags retrieved; aborting.")
            return

        versions = self.filter_versions(tags)
        if not versions:
            print("No matching r-base versions found.")
            return

        created = []
        for ver in versions:
            filename = f"r{ver}.def"
            path = os.path.join(self.base_dir, filename)
            if os.path.exists(path):
                # already have a script, skip
                print(f"Skipping existing {filename}")
                continue
            content = self.TPL.replace("${VERSION}", ver).replace("${CODE_NAME}", self.CODE_NAME)
            try:
                with open(path, "w") as of:
                    of.write(content)
                os.chmod(path, 0o755)
                created.append(path)
                print(f"Created {filename}")
            except Exception as e:
                print(f"Failed to write {filename}: {e}", file=sys.stderr)

        if created:
            print(f"Created {len(created)} build-scripts")
        else:
            print("No new build-scripts were created.")

if __name__ == "__main__":
    RVersions().run()
