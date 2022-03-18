# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# https://github.com/mcmtroffaes/sphinxcontrib-bibtex/blob/develop/test/roots/test-bibliography_style_label_1/conf.py
from pybtex.style.formatting.unsrt import Style as UnsrtStyle
from pybtex.style.labels import BaseLabelStyle
from pybtex.plugin import register_plugin

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import re
import sys

# sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath(os.path.join('..', '..')))

# -- Project information -----------------------------------------------------

project = "Faceted-Hyper-SARA"
copyright = "2021, P.-A. Thouvenin, A. Dabbech, M. Jiang, A. Jackson, Y. Wiaux, A. Abdulaziz"
author = "P.-A. Thouvenin, A. Dabbech, M. Jiang, A. Jackson, Y. Wiaux, A. Abdulaziz"

# The full version, including alpha/beta/rc tags
release = re.sub('^v', '', os.popen('git describe').read().strip())
# The short X.Y version.
version = release
# release = '0'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinxcontrib.matlab",  # support for Matlab
    "sphinx.ext.coverage",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx_rtd_theme",
    "sphinxcontrib.bibtex",
    "sphinx.ext.mathjax",  # LaTeX support
    "sphinx.ext.autosectionlabel",  # references in the same document 
]

autosummary_generate = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = ".rst"

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "lib/faceted-wavelet-transform",
    "lib/measurement-operator",
]
matlab_src_dir = os.path.abspath(os.path.join("..", ".."))
# primary_domain = "mat"
print("matlab_src_dir: {}".format(matlab_src_dir))

autodoc_mock_imports = ["im", "imager", "lsm", "mqt", "ms", "numpy", "pyfits", "pyrap", "Pyxis", "scipy", "Tigger", "astropy", "casacore", "time"]

# Bibliography
# .. https://sphinxcontrib-bibtex.readthedocs.io/en/latest/usage.html
bibtex_bibfiles = ["strings_all_ref.bib", "biblio.bib"]
bibtex_encoding = "utf-8-sig"
bibtex_default_style = 'mystyle'
# bibtex_reference_style =  # 'author_year'

# a simple label style which uses the bibtex keys for labels
class MyLabelStyle(BaseLabelStyle):

    def format_labels(self, sorted_entries):
        for entry in sorted_entries:
            yield entry.key


class MyStyle(UnsrtStyle):

    default_label_style = MyLabelStyle


register_plugin('pybtex.style.formatting', 'mystyle', MyStyle)


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

html_theme_options = {
    "logo_only": False,
    "display_version": True,
    "prev_next_buttons_location": "bottom",
    "style_external_links": False,
    # "vcs_pageview_mode": "",
    # 'style_nav_header_background': 'white',
    # Toc options
    "collapse_navigation": True,
    "sticky_navigation": True,
    "navigation_depth": 2,
    "includehidden": True,
    "titles_only": False,
}

# Napoleon settings
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = False
napoleon_use_rtype = True
napoleon_preprocess_types = False
napoleon_type_aliases = None
napoleon_attr_annotations = True
