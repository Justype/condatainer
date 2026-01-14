CondaTainer
=========================

`CondaTainer <https://github.com/Justype/condatainer>`_ is an HPC-oriented wrapper tool that leverages Apptainer, SquashFS, and Micromamba to enable the easy management of isolated single-environment files, software packages, reference data, and indices.

* **Inode Saver:** Packing 30k+ Conda files into a single image to satisfy inode quotas.
* **Web-App Ready:** Out-of-the-box support for running *RStudio Server* and *code-server*.
* **Slurm Integration:** Native compatibility with Slurm scheduler for batch job submission.

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

  Custom System Overlays <advanced_usage/condatainer_custom_def>
  Custom Build Script <advanced_usage/condatainer_custom_script>

.. toctree::
  :caption: Manuals:
  :hidden:
  :glob:

  manuals/*

.. toctree::
  :caption: Headless Server Tutorials:
  :hidden:

  RStudio Server on Headless <tutorials/rstudio-server_on_server>

