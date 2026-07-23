CondaTainer
=========================

`CondaTainer <https://github.com/Justype/condatainer>`_ is a rootless CLI for managing *tools* / *data* / *project environments* and launching apps on HPC — designed for individuals and small teams using institutional or regional compute resources.

* **Web-App Ready:** Launch *RStudio*, *VS Code*, *Jupyter* and more with one command.
* **Unified Management:** Easily organize group-level tools/data and isolate project environments.
* **Inode Saver:** Packing 30k+ Conda files into a single portable image to bypass quota limits.
* **LustreFS Friendly:** Stages conda creation on node-local SSD to reduce network filesystem load.
* **Scheduler Native:** Out-of-the-box integration with HPC scheduler, like *Slurm*.

.. toctree::
  :caption: User Guide:
  :maxdepth: 1

  Installation <user_guide/installation>
  Concepts: Overlays <user_guide/concepts>
  Module Overlays <user_guide/module_overlays>
  Writable Environment Overlays <user_guide/environment_overlays>
  Read-Only Bundle Overlays <user_guide/bundle_overlays>

.. toctree::
  :caption: Helpers:
  :maxdepth: 1
  :glob:

  Helper Scripts <helpers/helpers>
  helpers/*

.. toctree::
  :caption: FAQ:
  :maxdepth: 1
  :glob:

  qa/*

.. toctree::
  :caption: Advanced Usage:
  :maxdepth: 1
  :glob:

  advanced_usage/*

.. toctree::
  :caption: Deployment:
  :maxdepth: 1
  :glob:

  deployment/*

.. toctree::
  :caption: Manuals:
  :maxdepth: 1
  :glob:

  manuals/*
