# Configuration file for Sphinx

project = "RBA_cutsets"
author = "Jérémie Prokob"

extensions = [
    "nbsphinx",  # enables notebook support
]

extensions = [
    "sphinx.ext.duration",
    "sphinx.ext.doctest",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "nbsphinx",
    ]

templates_path = ["_templates"]
exclude_patterns = []

html_theme = "alabaster"  # or "sphinx_rtd_theme" if you prefer
 
