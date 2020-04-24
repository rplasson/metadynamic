# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import re

sys.path.insert(0, os.path.abspath("."))


# for docstring conversion from epytext
# thanks https://stackoverflow.com/questions/10012741/automated-way-to-switch-from-epydocs-docstring-formatting-to-sphinx-docstring-f

re_field = re.compile("@(param|type|rtype|return|raise)")


def fix_docstring(app, what, name, obj, options, lines):
    for i in range(len(lines)):
        lines[i] = re_field.sub(r":\1", lines[i])


def setup(app):
    app.connect("autodoc-process-docstring", fix_docstring)


# -- Project information -----------------------------------------------------

project = "metadynamic"
copyright = "2020, Raphaël Plasson"
author = "Raphaël Plasson"

# The full version, including alpha/beta/rc tags
release = "v1.0.1"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ["sphinx.ext.autodoc"]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "alabaster"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
