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
html_theme = 'sphinx_book_theme'
html_logo = '_static/logo_cnt.svg'
html_favicon = '_static/favicon.png'
html_static_path = ['_static']
html_title = 'CondaTainer'

html_theme_options = {
    # Repository buttons in the article header (replaces the RTD "Edit on
    # GitHub" html_context block and the old _templates/layout.html icon).
    'repository_url': 'https://github.com/Justype/condatainer',
    'repository_branch': 'main',
    'path_to_docs': 'docs',
    'use_repository_button': True,
    'use_edit_page_button': True,
    'use_issues_button': True,
    'home_page_in_toc': False,
    # Sidebar nav: show top-level entries, allow deeper pages to expand.
    'show_navbar_depth': 1,
    'max_navbar_depth': 4,
    'collapse_navbar': False,
    'logo': {
        'alt_text': 'CondaTainer',
    },
    'search_bar_text': 'Search the docs...',
}

# Sidebar in the theme's default order: 'search-button-field.html' under the
# logo opens the search overlay (Ctrl+K), the same as the mamba docs.
html_sidebars = {
    '**': [
        'navbar-logo.html',
        'search-button-field.html',
        'icon-links.html',
        'sbt-sidebar-nav.html',
    ],
}

# Custom CSS
html_css_files = [
    'custom.css',
]
