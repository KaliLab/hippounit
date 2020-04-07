.. HippoUnit documentation master file, created by
   sphinx-quickstart on Thu Dec 19 13:56:00 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

========================
HippoUnit: Documentation
========================

This package contains validation tests for models of hippocampus,
based on the SciUnit framework and the NeuronUnit package.

Current tests cover somatic behavior and signal propagation and integration in apical dendrites of
hippocampal CA1 pyramidal cell models.

Contents
========
.. toctree::
    :maxdepth: 1
    :numbered:

    page_introduction
    {% if test_json|length != 0 %}
    page_testLibrary
    {% endif %}
    page_implementation
    page_tests
    page_capabilities
    page_modelloader
    page_scores
    page_testing
    page_acknowledgement


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
