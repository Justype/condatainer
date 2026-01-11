CondaTainer
=========================

**CondaTainer** is a tool to create and manage lightweight, portable container-like environments using Conda and overlay filesystems. It is designed for use on HPC systems where traditional container technologies may not be available.

By encapsulating Conda environments within **Apptainer/Singularity overlay files**, CondaTainer solves the "small file problem" common on HPC parallel file systems (Lustre/GPFS) while preserving the flexibility of a writable environment.

.. note::
   **CondaTainer** is designed to work without root privileges, leveraging ``fakeroot`` and ``overlayfs`` to allow users to build and manage containers entirely in user-space.

The Tool Suite
--------------

The package consists of two distinct, parallel utilities:

- **CondaTainer**: Wraps environments in overlays using `Apptainer <https://apptainer.org/>`_ and `SquashFS <https://github.com/plougher/squashfs-tools>`_.
- **ModGen**: Automates `Lmod <https://lmod.readthedocs.io/en/latest/>`_/`Environment Modules <https://modules.readthedocs.io/en/latest/>`_ modules creation.

Which one to choose?

- If you care about inode usage and project-level isolation, use **CondaTainer**.
- If you want rstudio-server, igv, and other tools on HPC, use **CondaTainer**.
- If you prefer Lmod and using system available modules, use **ModGen**.

.. toctree::
  :caption: Installation:
  :maxdepth: 1

  CondaTainer <installation/condatainer>
  ModGen <installation/modgen>

.. toctree::
  :caption: User Guide:
  :maxdepth: 1
  
  Concepts: Modules <user_guide/concepts>
  CondaTainer: System Overlays <user_guide/condatainer_system_level>
  CondaTainer: Writable Project Env <user_guide/condatainer_project_level_writable>
  CondaTainer: Read-Only Project Env <user_guide/condatainer_project_level_readonly>
  ModGen: System Modules <user_guide/modgen_system_level>

.. toctree::
  :caption: Tutorials:
  :maxdepth: 1

  RStudio Server on HPC <tutorials/rstudio-server_on_HPC>
  RStudio Server on Headless <tutorials/rstudio-server_on_server>
  code-server on HPC <tutorials/code-server_on_HPC>
  VS Code Tunnel on HPC <tutorials/vscode-tunnel_on_HPC>

.. toctree::
  :caption: FAQ:
  :maxdepth: 1
  :glob:

  qa/*

.. toctree::
  :caption: Manuals:
  :hidden:
  :glob:

  manuals/*
