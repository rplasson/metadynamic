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
extensions = [
    "sphinx.ext.autodoc",
    "sphinx_autodoc_typehints",
    "sphinxcontrib.autoprogram",
    "cloud_sptheme.ext.autodoc_sections"
]

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
html_theme = "greencloud"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# -- Options for LaTeX output -------------------------------------------------

latex_elements = {
    "papersize": "a4paper",
    "fncychap": "\\usepackage{fncychap}",
    "fontpkg": "\\usepackage{amsmath,amsfonts,amssymb,amsthm}",
    "pointsize": "10pt",
    "preamble": r"""
        \usepackage{amsmath,amsfonts,amssymb,amsthm}
        \usepackage{graphicx}
        \usepackage{ebgaramond}
        \usepackage{fourier-orns}

        \usepackage{pgfornament}

         \def\cinzelfamily{Cinzel-LF}
         \newcommand*\cinzel{\fontfamily{\cinzelfamily}\def\itshape{\fontfamily{CinzelDecorative-LF}\fontshape{n}\selectfont}\selectfont}
         \newcommand*\cinzelblack{\fontfamily{\cinzelfamily}\fontseries{k}\def\itshape{\fontfamily{CinzelDecorative-LF}\fontshape{n}\selectfont}\selectfont}

        \usepackage{color}
        \usepackage{transparent}
        \usepackage{eso-pic}
        \usepackage{lipsum}
        \usepackage[hyphenation, parindent, lastparline, rivers]{impnattypo}
        \usepackage[all]{nowidow}
        \usepackage[activate={true, nocompatibility},babel=true, tracking=true,final]{microtype}

         \definecolor{bl}{rgb}{0.10,0.2,0.6}
         \definecolor{dbl}{rgb}{0.05,0.1,0.3}
         \definecolor{rd}{rgb}{0.6,0.1,0.3}
         \definecolor{drd}{rgb}{0.3,0.05,0.1}
         \definecolor{gn}{rgb}{0.2,0.6,0.1}
         \definecolor{dgn}{rgb}{0.1,0.3,0.05}
    """,
    "maketitle": r"""
    \begin{titlepage}
      \begin{center}
        \pgfornament[color=dgn,width=\textwidth]{88}\\[0.3cm]
        \color{dgn} \Huge \cinzelblack \itshape{} Metadynamic
        \normalfont{}\\
        \decofourright{}\\
        \Large \textsw{\oldstylenums{Université d'Avignon, UMR408}}
        \pgfornament[color=dgn,width=\textwidth]{88} \vfill{}
        \color{dbl} \large \textsc{Documentation and API}\\
        \large by\\
        \huge  \textsw{Raphaël Plasson}\\[0.5cm]
        \large \textit{\oldstylenums{\today}} \vfill{}
        \pgfornament[color=dgn,width=0.75\textwidth]{89}
      \end{center}
    \end{titlepage}
    """,
    "sphinxsetup": """
    hmargin={1.75cm,2.475cm},
    vmargin={2.475cm,3.5cm},
    verbatimwithframe=true,
    TitleColor={rgb}{0.05,0.1,0.3},
    InnerLinkColor={rgb}{0.1,0.3,0.05},
    OuterLinkColor={rgb}{0.3,0.05,0.1}
    """,
}
