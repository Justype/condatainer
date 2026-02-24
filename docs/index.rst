CondaTainer
=========================

`CondaTainer <https://github.com/Justype/condatainer>`_  is a rootless CLI designed to manage tools, data, and project environments, and seamlessly launch interactive apps on HPC clusters.

* **Web-App Ready:** Launch *RStudio*, *VS Code*, *Jupyter* and more with one command.
* **Unified Management:** Easily organize group-level tools/data and isolate project environments.
* **Inode Saver:** Packing 30k+ Conda files into a single portable image to bypass quota limits.
* **Scheduler Native:** Out-of-the-box integration with *Slurm*, *PBS*, *LSF*, and *HTCondor*.

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

  Helper Scripts Overview <tutorials/helpers_on_HPC>
  RStudio Server <tutorials/rstudio-server_on_HPC>
  RStudio Server with Conda <tutorials/rstudio-server-conda_on_HPC>
  VS Code Server <tutorials/vscode-server_on_HPC>
  VS Code Tunnel <tutorials/vscode-tunnel_on_HPC>
  XFCE4 Desktop <tutorials/xfce4_on_HPC>

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
  :maxdepth: 1
  :glob:

  manuals/*

.. toctree::
  :caption: Other Tutorials:
  :hidden:

  code-server on HPC <tutorials/code-server_on_HPC>
  Helper Scripts on Headless <tutorials/helpers_on_server>
  RStudio Server on Headless <tutorials/rstudio-server_on_server>
