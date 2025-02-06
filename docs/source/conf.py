'''
Date: 2024-08-24 22:24:47
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-02-06 16:19:36
Description: 
'''
# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'LazyDock'
copyright = '2024, BHM-Bob'
author = 'BHM-Bob G'

release = '0.11.0'
version = '0.11.0'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    "sphinx.ext.recommonmark",
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'