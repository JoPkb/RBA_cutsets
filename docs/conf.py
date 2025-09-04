# Configuration file for Sphinx

project = "My Project"
author = "Your Name"

extensions = [
    "myst_nb",  # enables notebook support
]

templates_path = ["_templates"]
exclude_patterns = []

html_theme = "alabaster"  # or "sphinx_rtd_theme" if you prefer
jupyter_execute_notebooks = "auto" 
