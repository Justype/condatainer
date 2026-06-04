CondaTainer
=========================

`CondaTainer <https://github.com/Justype/condatainer>`_ is a rootless CLI for managing *tools* / *data* / *project environments* and launching apps on HPC — designed for individuals and small teams using institutional or regional compute resources.

* **Web-App Ready:** Launch *RStudio*, *VS Code*, *Jupyter* and more with one command.
* **Unified Management:** Easily organize group-level tools/data and isolate project environments.
* **Inode Saver:** Packing 30k+ Conda files into a single portable image to bypass quota limits.
* **LustreFS Friendly:** Stages conda creation on node-local SSD to reduce network filesystem load.
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
  :caption: Helpers:
  :maxdepth: 1

  Helper Scripts <tutorials/helpers>
  RStudio Server <tutorials/rstudio-server>
  VS Code <tutorials/vscode>
  XFCE4 Desktop <tutorials/xfce4>

.. toctree::
  :caption: FAQ:
  :maxdepth: 1
  :glob:

  qa/*

.. toctree::
  :caption: Advanced Usage:
  :maxdepth: 1

  advanced_usage/*

.. toctree::
  :caption: Manuals:
  :maxdepth: 1
  :glob:

  manuals/*
