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
release = "v1.0.2"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    #    "sphinx_autodoc_typehints",
    "sphinx.ext.inheritance_diagram",
    "sphinxcontrib.autoprogram",
    "cloud_sptheme.ext.autodoc_sections",
    "nbsphinx",
    # "sphinx.ext.imgconverter",
    "sphinxcontrib.rsvgconverter",
]

rsvg_converter_args = ["--format=pdf"]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

inheritance_graph_attrs = dict(fontsize=32, size='"16.0, 20.0"')

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
        \usepackage[T1]{fontenc}
        \usepackage{anyfontsize}
        %\usepackage{amsmath,amsfonts,amssymb,amsthm}
        \usepackage{graphicx}
        \usepackage[cmintegrals,cmbraces]{newtxmath}
        \usepackage{ebgaramond-maths}
        %\usepackage{ascii}
        %\renewcommand*\ttdefault{lcmtt}
        \usepackage[osf,scale=0.85]{sourcecodepro}
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

    % Pygments definitions
    
     \makeatletter
     \def\PY@reset{\let\PY@it=\relax \let\PY@bf=\relax%
         \let\PY@ul=\relax \let\PY@tc=\relax%
         \let\PY@bc=\relax \let\PY@ff=\relax}
     \def\PY@tok#1{\csname PY@tok@#1\endcsname}
     \def\PY@toks#1+{\ifx\relax#1\empty\else%
         \PY@tok{#1}\expandafter\PY@toks\fi}
     \def\PY@do#1{\PY@bc{\PY@tc{\PY@ul{%
         \PY@it{\PY@bf{\PY@ff{#1}}}}}}}
     \def\PY#1#2{\PY@reset\PY@toks#1+\relax+\PY@do{#2}}
     
     \expandafter\def\csname PY@tok@w\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.73,0.73,0.73}{##1}}}
     \expandafter\def\csname PY@tok@c\endcsname{\let\PY@it=\textit\def\PY@tc##1{\textcolor[rgb]{0.25,0.50,0.50}{##1}}}
     \expandafter\def\csname PY@tok@cp\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.74,0.48,0.00}{##1}}}
     \expandafter\def\csname PY@tok@k\endcsname{\let\PY@bf=\textbf\def\PY@tc##1{\textcolor[rgb]{0.00,0.50,0.00}{##1}}}
     \expandafter\def\csname PY@tok@kp\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.00,0.50,0.00}{##1}}}
     \expandafter\def\csname PY@tok@kt\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.69,0.00,0.25}{##1}}}
     \expandafter\def\csname PY@tok@o\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.40,0.40,0.40}{##1}}}
     \expandafter\def\csname PY@tok@ow\endcsname{\let\PY@bf=\textbf\def\PY@tc##1{\textcolor[rgb]{0.67,0.13,1.00}{##1}}}
     \expandafter\def\csname PY@tok@nb\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.00,0.50,0.00}{##1}}}
     \expandafter\def\csname PY@tok@nf\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.00,0.00,1.00}{##1}}}
     \expandafter\def\csname PY@tok@nc\endcsname{\let\PY@bf=\textbf\def\PY@tc##1{\textcolor[rgb]{0.00,0.00,1.00}{##1}}}
     \expandafter\def\csname PY@tok@nn\endcsname{\let\PY@bf=\textbf\def\PY@tc##1{\textcolor[rgb]{0.00,0.00,1.00}{##1}}}
     \expandafter\def\csname PY@tok@ne\endcsname{\let\PY@bf=\textbf\def\PY@tc##1{\textcolor[rgb]{0.82,0.25,0.23}{##1}}}
     \expandafter\def\csname PY@tok@nv\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.10,0.09,0.49}{##1}}}
     \expandafter\def\csname PY@tok@no\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.53,0.00,0.00}{##1}}}
     \expandafter\def\csname PY@tok@nl\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.63,0.63,0.00}{##1}}}
     \expandafter\def\csname PY@tok@ni\endcsname{\let\PY@bf=\textbf\def\PY@tc##1{\textcolor[rgb]{0.60,0.60,0.60}{##1}}}
     \expandafter\def\csname PY@tok@na\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.49,0.56,0.16}{##1}}}
     \expandafter\def\csname PY@tok@nt\endcsname{\let\PY@bf=\textbf\def\PY@tc##1{\textcolor[rgb]{0.00,0.50,0.00}{##1}}}
     \expandafter\def\csname PY@tok@nd\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.67,0.13,1.00}{##1}}}
     \expandafter\def\csname PY@tok@s\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.73,0.13,0.13}{##1}}}
     \expandafter\def\csname PY@tok@sd\endcsname{\let\PY@it=\textit\def\PY@tc##1{\textcolor[rgb]{0.73,0.13,0.13}{##1}}}
     \expandafter\def\csname PY@tok@si\endcsname{\let\PY@bf=\textbf\def\PY@tc##1{\textcolor[rgb]{0.73,0.40,0.53}{##1}}}
     \expandafter\def\csname PY@tok@se\endcsname{\let\PY@bf=\textbf\def\PY@tc##1{\textcolor[rgb]{0.73,0.40,0.13}{##1}}}
     \expandafter\def\csname PY@tok@sr\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.73,0.40,0.53}{##1}}}
     \expandafter\def\csname PY@tok@ss\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.10,0.09,0.49}{##1}}}
     \expandafter\def\csname PY@tok@sx\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.00,0.50,0.00}{##1}}}
     \expandafter\def\csname PY@tok@m\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.40,0.40,0.40}{##1}}}
     \expandafter\def\csname PY@tok@gh\endcsname{\let\PY@bf=\textbf\def\PY@tc##1{\textcolor[rgb]{0.00,0.00,0.50}{##1}}}
     \expandafter\def\csname PY@tok@gu\endcsname{\let\PY@bf=\textbf\def\PY@tc##1{\textcolor[rgb]{0.50,0.00,0.50}{##1}}}
     \expandafter\def\csname PY@tok@gd\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.63,0.00,0.00}{##1}}}
     \expandafter\def\csname PY@tok@gi\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.00,0.63,0.00}{##1}}}
     \expandafter\def\csname PY@tok@gr\endcsname{\def\PY@tc##1{\textcolor[rgb]{1.00,0.00,0.00}{##1}}}
     \expandafter\def\csname PY@tok@ge\endcsname{\let\PY@it=\textit}
     \expandafter\def\csname PY@tok@gs\endcsname{\let\PY@bf=\textbf}
     \expandafter\def\csname PY@tok@gp\endcsname{\let\PY@bf=\textbf\def\PY@tc##1{\textcolor[rgb]{0.00,0.00,0.50}{##1}}}
     \expandafter\def\csname PY@tok@go\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.53,0.53,0.53}{##1}}}
     \expandafter\def\csname PY@tok@gt\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.00,0.27,0.87}{##1}}}
     \expandafter\def\csname PY@tok@err\endcsname{\def\PY@bc##1{\setlength{\fboxsep}{0pt}\fcolorbox[rgb]{1.00,0.00,0.00}{1,1,1}{\strut ##1}}}
     \expandafter\def\csname PY@tok@kc\endcsname{\let\PY@bf=\textbf\def\PY@tc##1{\textcolor[rgb]{0.00,0.50,0.00}{##1}}}
     \expandafter\def\csname PY@tok@kd\endcsname{\let\PY@bf=\textbf\def\PY@tc##1{\textcolor[rgb]{0.00,0.50,0.00}{##1}}}
     \expandafter\def\csname PY@tok@kn\endcsname{\let\PY@bf=\textbf\def\PY@tc##1{\textcolor[rgb]{0.00,0.50,0.00}{##1}}}
     \expandafter\def\csname PY@tok@kr\endcsname{\let\PY@bf=\textbf\def\PY@tc##1{\textcolor[rgb]{0.00,0.50,0.00}{##1}}}
     \expandafter\def\csname PY@tok@bp\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.00,0.50,0.00}{##1}}}
     \expandafter\def\csname PY@tok@fm\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.00,0.00,1.00}{##1}}}
     \expandafter\def\csname PY@tok@vc\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.10,0.09,0.49}{##1}}}
     \expandafter\def\csname PY@tok@vg\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.10,0.09,0.49}{##1}}}
     \expandafter\def\csname PY@tok@vi\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.10,0.09,0.49}{##1}}}
     \expandafter\def\csname PY@tok@vm\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.10,0.09,0.49}{##1}}}
     \expandafter\def\csname PY@tok@sa\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.73,0.13,0.13}{##1}}}
     \expandafter\def\csname PY@tok@sb\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.73,0.13,0.13}{##1}}}
     \expandafter\def\csname PY@tok@sc\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.73,0.13,0.13}{##1}}}
     \expandafter\def\csname PY@tok@dl\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.73,0.13,0.13}{##1}}}
     \expandafter\def\csname PY@tok@s2\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.73,0.13,0.13}{##1}}}
     \expandafter\def\csname PY@tok@sh\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.73,0.13,0.13}{##1}}}
     \expandafter\def\csname PY@tok@s1\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.73,0.13,0.13}{##1}}}
     \expandafter\def\csname PY@tok@mb\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.40,0.40,0.40}{##1}}}
     \expandafter\def\csname PY@tok@mf\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.40,0.40,0.40}{##1}}}
     \expandafter\def\csname PY@tok@mh\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.40,0.40,0.40}{##1}}}
     \expandafter\def\csname PY@tok@mi\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.40,0.40,0.40}{##1}}}
     \expandafter\def\csname PY@tok@il\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.40,0.40,0.40}{##1}}}
     \expandafter\def\csname PY@tok@mo\endcsname{\def\PY@tc##1{\textcolor[rgb]{0.40,0.40,0.40}{##1}}}
     \expandafter\def\csname PY@tok@ch\endcsname{\let\PY@it=\textit\def\PY@tc##1{\textcolor[rgb]{0.25,0.50,0.50}{##1}}}
     \expandafter\def\csname PY@tok@cm\endcsname{\let\PY@it=\textit\def\PY@tc##1{\textcolor[rgb]{0.25,0.50,0.50}{##1}}}
     \expandafter\def\csname PY@tok@cpf\endcsname{\let\PY@it=\textit\def\PY@tc##1{\textcolor[rgb]{0.25,0.50,0.50}{##1}}}
     \expandafter\def\csname PY@tok@c1\endcsname{\let\PY@it=\textit\def\PY@tc##1{\textcolor[rgb]{0.25,0.50,0.50}{##1}}}
     \expandafter\def\csname PY@tok@cs\endcsname{\let\PY@it=\textit\def\PY@tc##1{\textcolor[rgb]{0.25,0.50,0.50}{##1}}}
     
     \def\PYZbs{\char`\\}
     \def\PYZus{\char`\_}
     \def\PYZob{\char`\{}
     \def\PYZcb{\char`\}}
     \def\PYZca{\char`\^}
     \def\PYZam{\char`\&}
     \def\PYZlt{\char`\<}
     \def\PYZgt{\char`\>}
     \def\PYZsh{\char`\#}
     \def\PYZpc{\char`\%}
     \def\PYZdl{\char`\$}
     \def\PYZhy{\char`\-}
     \def\PYZsq{\char`\'}
     \def\PYZdq{\char`\"}
     \def\PYZti{\char`\~}
     % for compatibility with earlier versions
     \def\PYZat{@}
     \def\PYZlb{[}
     \def\PYZrb{]}
     \makeatother
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
