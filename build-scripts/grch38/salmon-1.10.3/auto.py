#!/usr/bin/env python3

import os

script_path = __file__

def get_assembly():
    return os.path.basename(os.path.dirname(os.path.dirname(script_path)))

def get_gencode_versions():
    """Return a sorted list of available gencode version numbers found under
    the assembly's `transcript-gencode/` directory (e.g. ['44', '47'])."""
    assembly_dir = os.path.dirname(os.path.dirname(script_path))
    gtf_dir = os.path.join(assembly_dir, 'transcript-gencode')
    versions = []
    if not os.path.isdir(gtf_dir):
        return versions
    for name in sorted(os.listdir(gtf_dir)):
        if name.startswith('template') or name.endswith('.sh') or name.endswith('.py'):
            continue
        versions.append(name)
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

    gencode_versions = get_gencode_versions()
    if not gencode_versions:
        print("No new GENCODE versions found. Nothing to do.")
        return

    created = []
    for gv in gencode_versions:
        # gv might be like '44' or 'M36' — embed as-is
        out_name = f"gencode{gv}"
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