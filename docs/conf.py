# Configuration file for the Sphinx documentation builder.

import os
import sys
import datetime

# -- Project information -----------------------------------------------------
project = 'CondaTainer'
copyright = f"{datetime.date.today().year}, Justype"
author = 'Justype'

# -- General configuration ---------------------------------------------------
extensions = [
    'myst_parser',
    'sphinx_copybutton',
]

templates_path = ['_templates']
# Ignore original manuals (they are included via wrapper pages)
# Exclude README files so Sphinx does not pick up repository README.* files
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**/README*', 'README*']

# Support Markdown via MyST
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

# MyST settings
myst_enable_extensions = [
    'deflist',
    'html_admonition',
    'html_image',
    'colon_fence',
]
myst_heading_anchors = 4

# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme'
html_logo = '_static/logo_cnt.svg'
html_favicon = '_static/favicon.png'
html_static_path = ['_static']

html_theme_options = {
    'logo_only': True,
    'collapse_navigation': False,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'version_selector': False,
}

# Add GitHub repository info so the theme can link to the project
html_context = {
    'display_github': True,  # Integrate 'Edit on GitHub' links
    'github_user': 'Justype',
    'github_repo': 'condatainer',
    'github_version': 'main',
    'conf_py_path': '/docs/',
}

# Custom CSS
html_css_files = [
    'custom.css',
]