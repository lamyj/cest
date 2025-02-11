import os
import sys

project = "CEST Toolbox"
copyright = "2024, ICube, University of Strasbourg-CNRS"
author = "Julien Lamy"

extensions = ["sphinx.ext.autodoc", "sphinx.ext.napoleon"]

napoleon_preprocess_types = True

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

autodoc_default_options = {
    "member-order": "bysource"
}

html_theme = "furo"
html_static_path = ["_static"]
html_css_files = [
    "css/style.css",
]
html_show_sphinx = False

