#!/bin/bash
# Run all auto.py scripts under build-scripts/
set -e

find "$(dirname "$0")/../../build-scripts" -type f -name 'auto.py' | while read -r script; do
    echo "Running $script..."
    python3 "$script"
done

# Run it twice (because some scripts depend on others)
# e.g. star index depends on gtf-gencode

find "$(dirname "$0")/../../build-scripts" -type f -name 'auto.py' | while read -r script; do
    echo "Running $script..."
    python3 "$script"
done
