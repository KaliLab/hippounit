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
sys.path.insert(0, os.path.abspath('../..'))

# The master toctree document.
master_doc = 'index'
source_suffix = '.rst'

# -- Project information -----------------------------------------------------

project = 'HippoUnit'
copyright = u'2018, Sara Saray, Shailesh Appukuttan, Szabolcs Káli, Andrew Davison, Christian Rossert'
author = u'Sara Saray, Shailesh Appukuttan, Szabolcs Káli, Andrew Davison, Christian Rossert'

# The full version, including alpha/beta/rc tags
release = '1.3.3'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
import sphinx_automodapi

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx_automodapi.automodapi'
]

autodoc_default_options = {
    'imported-members':True
}
add_module_names = False
# html_sidebars = { '**': ['globaltoc.html', 'relations.html', 'sourcelink.html', 'searchbox.html'] }
html_sidebars = { '**': ['globaltoc.html', 'sourcelink.html', 'searchbox.html'] } # remove next/previous links on sidebar

autosummary_generate = True
numpydoc_show_class_members = False

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'nature'
pygments_style = 'sphinx'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# html_theme_options = {
#     'page_width': 'auto',
#     'body_max_width': 'auto'
# }

# These paths are either relative to html_static_path
# or fully qualified paths (eg. https://...)
html_css_files = [
    'css/custom.css',
]

# ----------------------------------------------------------------------------

def rstjinja(app, docname, source):
    """
    Render our pages as a jinja template for fancy templating goodness.
    """
    # Make sure we're outputting HTML
    if app.builder.format != 'html':
        return
    src = source[0]
    rendered = app.builder.templates.render_string(
        src, app.config.html_context
    )
    source[0] = rendered

def setup(app):
    app.connect("source-read", rstjinja)

import inspect, importlib, sciunit
package_import_name = "hippounit"

tests_module = "{}.tests".format(package_import_name)
module = importlib.import_module(tests_module)
test_classes = [x[0] for x in inspect.getmembers(module,
                    lambda member: inspect.isclass(member)
                                    and not(issubclass(member, sciunit.Test)
                                            and (tests_module in member.__module__)))]

capabilities_module = "{}.capabilities".format(package_import_name)
module = importlib.import_module(capabilities_module)
capability_classes = [x[0] for x in inspect.getmembers(module,
                    lambda member: inspect.isclass(member)
                                    and not(issubclass(member, sciunit.Capability)
                                            and (capabilities_module in member.__module__)))]

scores_module = "{}.scores".format(package_import_name)
module = importlib.import_module(scores_module)
score_classes = [x[0] for x in inspect.getmembers(module,
                    lambda member: inspect.isclass(member)
                                    and not(issubclass(member, sciunit.Score)
                                            and (scores_module in member.__module__)))]


import json
test_json = {}
print(os.getcwd())
if os.path.isfile("VF_test_info.json"):
    with open("VF_test_info.json", "r") as f:
        test_json = json.load(f)

html_context = {
    'test_classes'       : test_classes,
    'capability_classes' : capability_classes,
    'score_classes'      : score_classes,
    'test_json'          : test_json
}

# ----------------------------------------------------------------------------
