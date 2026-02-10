#!/usr/bin/env python3

import os

script_path = __file__

def get_assembly():
    return os.path.basename(os.path.dirname(os.path.dirname(script_path)))

def get_new_gencode_versions():
    """
    Try to discover available GENCODE releases for this assembly by
    querying the EBI FTP listing.
    Only return versions that have not yet been created
    """
    import urllib.request
    import re

    base_url = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/'

    versions = []
    if base_url:
        try:
            with urllib.request.urlopen(base_url, timeout=10) as resp:
                txt = resp.read().decode('utf-8', errors='ignore')
            # Look for directory names like release_44/ or release_M36/
            found = set(re.findall(r'release_([0-9M]+)[/"\'>]', txt))
            for f in sorted(found):
                # Also check if gencode.v{version}.primary_assembly.annotation.gtf.gz exists
                try:
                    if not f.startswith('M'):
                        continue
                    f_int = int(f[1:])
                    if f_int < 6:
                        continue

                    if os.path.exists(f"{os.path.dirname(script_path)}/{f}"):
                        continue

                    gtf_url = f"{base_url}release_{f}/gencode.v{f}.transcripts.fa.gz"
                    with urllib.request.urlopen(gtf_url, timeout=5) as gtf_resp:
                        if gtf_resp.status == 200:
                            versions.append(f)
                except Exception:
                    pass
        except Exception as e:
            print(f"Warning: failed to query {base_url}: {e}")

    return versions

def main():
    assembly = get_assembly()
    print(f"Assembly: {assembly}")

    # Read template
    base_dir = os.path.dirname(script_path)
    template_path = os.path.join(base_dir, 'template')
    if not os.path.isfile(template_path):
        print(f"Template not found at {template_path}")
        return
    with open(template_path, 'r') as tf:
        tpl = tf.read()

    gencode_versions = get_new_gencode_versions()
    if not gencode_versions:
        print("No new GENCODE versions found. Nothing to do.")
        return

    created = []
    for gv in gencode_versions:
        # gv might be like '44' or 'M36' — embed as-is
        out_name = f"{gv}"
        out_path = os.path.join(base_dir, out_name)
        if os.path.exists(out_path):
            print(f"Skipping existing {out_name}")
            continue

        content = tpl
        content = content.replace('${ASSEMBLY}', assembly)
        content = content.replace('${GENCODE_VERSION}', str(gv))

        try:
            with open(out_path, 'w') as of:
                of.write(content)
            os.chmod(out_path, 0o755)
            created.append(out_path)
            print(f"Created {out_path}")
        except Exception as e:
            print(f"Failed to write {out_path}: {e}")

    if not created:
        print("No new build-scripts were created.")
    else:
        print(f"Created {len(created)} build-scripts")

if __name__ == "__main__":
    main()