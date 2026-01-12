CondaTainer
=========================

`CondaTainer <https://github.com/Justype/condatainer>`_ is an HPC-oriented wrapper tool that leverages Apptainer, SquashFS, OverlayFS, and Micromamba to enable the easy management of isolated single-environment files, software packages, reference data, and indices.

.. note::
   * **Solve the Small File Problem:** CondaTainer packs a Conda environment (30k+ inodes) into a **single, portable file**, reducing load on HPC parallel filesystems (Lustre/GPFS).
   * **Seamless Interactive Apps:** It simplifies the deployment of web-based applications, allowing you to easily set up *RStudio Server* and *code-server* without root privileges.

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
  :caption: Advanced Usage:
  :maxdepth: 1
  
  Custom Apptainer DEF <advanced_usage/condatainer_custom_def>
  Custom Build Script <advanced_usage/condatainer_custom_script>

.. toctree::
  :caption: Manuals:
  :hidden:
  :glob:

  manuals/*
