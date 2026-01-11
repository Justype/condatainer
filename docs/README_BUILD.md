This file describes how to build the Sphinx docs locally.

Prerequisites
- Python 3.10 (matches .readthedocs.yml)
- pip

Quick local build (recommended: use a virtualenv)

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r docs/requirements.txt
# Build HTML
sphinx-build -b html docs docs/_build/html
# Open the built site
xdg-open docs/_build/html/index.html || true
```

For live-reloading during editing:

```bash
sphinx-autobuild docs docs/_build/html --open-browser
```

If your repo contains a `Makefile` in `docs/`, you can instead run:

```bash
cd docs
make html
```

Notes
- `.readthedocs.yml` specifies Python 3.10 and installs `docs/requirements.txt` on RTD.
- If RTD fails during build, check the RTD build log for missing system packages (e.g., libxml2-dev) and add them under the `build` section of `.readthedocs.yml`.
