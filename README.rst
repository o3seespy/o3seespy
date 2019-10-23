.. image:: https://travis-ci.org/eng-tools/o3seespy.svg?branch=master
   :target: https://travis-ci.org/eng-tools/o3seespy
   :alt: Testing Status

.. image:: https://img.shields.io/pypi/v/o3seespy.svg
   :target: https://pypi.python.org/pypi/o3seespy
   :alt: PyPi version

.. image:: https://coveralls.io/repos/github/eng-tools/o3seespy/badge.svg
   :target: https://coveralls.io/github/eng-tools/o3seespy

.. image:: https://img.shields.io/badge/license-MIT-blue.svg
    :target: https://github.com/eng-tools/o3seespy/blob/master/LICENSE
    :alt: License

.. image:: https://eng-tools.github.io/static/img/ecp-badge.svg
    :target: https://eng-tools.github.io
    :alt: ECP project

.. image:: https://zenodo.org/badge/125842866.svg
   :target: https://zenodo.org/badge/latestdoi/125842866
   :alt: DOI

********
o3seespy
********

Object-orientated native python version of Opensees

Features
========

This package provides:

 * A 'pythonic' version of opensees - in that the input parameters to create opensees objects are all lowercase key value arguments
 * Object-orientated implementation of Opensees. Python is object-orientated and the underlying opensees C++ source code is object orientated.
 * Type checking of before calling C++ opensees code, so that python debugging and errors can be viewed
 * Additional features for using opensees in python:
    - saving and loading data directly from opensees into numpy arrays
    - saving and loading data directly from opensees into json files
    - Save an entire model as a json file
 * All object numbering handled by objects - no need for number tags!



How to Use
==========

`pip install o3seespy`

Examples
--------


Useful material
===============

*

Contributing
============

How do I get set up?
--------------------

1. Run ``pip install -r requirements.txt``


Package conventions
-------------------

* All names should be the same as the opensees tcl version, except:
    - The name should be converted to snake_case for a parameter or function
    - The name should be converted to CamelCase for an Object
    - The name should be converted to ALL_CAPS for static variables
    - If the name matches a python special name (e.g. lambda, in) then it should be adjusted according to the dictionary
    - Objects should be namespaced based on the object type (e.g. element, material)
    - For parameter that are used across many objects (e.g. atmospheric pressure) a standard name should be used



Testing
-------

Tests are run with pytest

* Locally run: ``pytest`` on the command line.

* Tests are run on every push using travis, see the ``.travis.yml`` file


Deployment
----------

To deploy the package to pypi.com you need to:

1. Push to the *pypi* branch. This executes the tests on circleci.com

2. Create a git tag and push to github, run: ``trigger_deploy.py`` or manually:

.. code:: bash

    git tag 0.5.2 -m "version 0.5.2"
    git push --tags origin pypi


Documentation
-------------

Built via Sphinx following: https://codeandchaos.wordpress.com/2012/07/30/sphinx-autodoc-tutorial-for-dummies/

For development mode

 1. cd to docs
 2. Run ``make html``

Docstrings follow numpy convention (in progress): https://numpydoc.readthedocs.io/en/latest/format.html

To fix long_description in setup.py: ``pip install collective.checkdocs``, ``python setup.py checkdocs``
