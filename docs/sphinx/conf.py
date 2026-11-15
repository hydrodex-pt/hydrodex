# Sphinx configuration for the documentation of the HydróDex Python interface.

import importlib.util
import os
import sys

# The docstrings are read from the compiled hydrodex module. CMake (target
# docs_python) passes the module it built via HYDRODEX_MODULE, which is loaded
# explicitly so that it takes precedence over any installed (e.g. editable)
# version.
module_file = os.environ.get("HYDRODEX_MODULE")
if module_file:
    spec = importlib.util.spec_from_file_location("hydrodex", module_file)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load the hydrodex module from {module_file}")
    module = importlib.util.module_from_spec(spec)
    sys.modules["hydrodex"] = module
    spec.loader.exec_module(module)

import hydrodex

project = "HydróDex"
author = "Giulio Barni, Eric Madge"
copyright = "2026, Giulio Barni, Eric Madge"
version = hydrodex.__version__
release = hydrodex.__version__

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
]

# docstrings are written in NumPy style
napoleon_google_docstring = False
napoleon_numpy_docstring = True

autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "member-order": "bysource",
}

html_theme = "furo"
html_title = f"HydróDex {release} Python API"
