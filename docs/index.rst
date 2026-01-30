CondaTainer
=========================

`CondaTainer <https://github.com/Justype/condatainer>`_  is an HPC-oriented CLI that streamlines environment and data management by utilizing Apptainer, SquashFS, OverlayFS, and Micromamba.

* **Inode Saver:** Packing 30k+ Conda files into a single image to satisfy inode quotas.
* **Web-App Ready:** Out-of-the-box support for running *RStudio Server*, *code-server* and more on HPC.
* **Environment Isolation:** Supports read-only (`.sqf`) for production and writable (`.img`) for development.
* **Workload Manager Integration:** Native compatibility with Slurm scheduler for batch job submission.

.. toctree::
  :caption: User Guide:
  :maxdepth: 1

  Installation <user_guide/installation>
  Concepts: Overlays <user_guide/concepts>
  Module Overlays <user_guide/module_overlays>
  Writable Workspace Overlays <user_guide/workspace_overlays>
  Read-Only Bundle Overlays <user_guide/bundle_overlays>

.. toctree::
  :caption: HPC Tutorials:
  :maxdepth: 1

  RStudio Server on HPC <tutorials/rstudio-server_on_HPC>
  code-server on HPC <tutorials/code-server_on_HPC>
  VS Code Tunnel on HPC <tutorials/vscode-tunnel_on_HPC>

.. toctree::
  :caption: FAQ:
  :maxdepth: 1
  :glob:

  qa/*

.. toctree::
  :caption: Advanced Usage:
  :maxdepth: 1

  Custom OS Overlays <advanced_usage/custom_os>
  Custom Bundle Overlays <advanced_usage/custom_bundle>

.. toctree::
  :caption: Manuals:
  :hidden:
  :glob:

  manuals/*

.. toctree::
  :caption: Headless Server Tutorials:
  :hidden:

  RStudio Server on Headless <tutorials/rstudio-server_on_server>

