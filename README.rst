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

Example: Inelastic SDOF
-----------------------


.. code:: python

    import numpy as np
    import matplotlib.pyplot as plt

    import openseespy.opensees as opy
    import o3seespy as opw

    from tests.conftest import TEST_DATA_DIR

    # Load a ground motion
    dt = 0.01
    rec = np.loadtxt(TEST_DATA_DIR + 'test_motion_dt0p01.txt')

    # Define inelastic SDOF
    period = 1.0
    xi = 0.05
    mass = 1.0
    f_yield = 1.5  # Reduce this to make it nonlinear
    r_post = 0.0

    # Initialise OpenSees instance
    osi = opw.OpenseesInstance(dimensions=2, state=0)

    # Establish nodes
    bot_node = opw.node.Node(osi, 0, 0)
    top_node = opw.node.Node(osi, 0, 0)

    # Fix bottom node
    opw.Fix(osi, top_node, opw.static.FREE, opw.static.FIXED, opw.static.FIXED)
    opw.Fix(osi, bot_node, opw.static.FIXED, opw.static.FIXED, opw.static.FIXED)
    # Set out-of-plane DOFs to be slaved
    opw.EqualDOF(osi, top_node, bot_node, [opw.static.Y, opw.static.ROTZ])

    # nodal mass (weight / g):
    opw.Mass(osi, top_node, mass, 0., 0.)

    # Define material
    k_spring = 4 * np.pi ** 2 * mass / period ** 2
    bilinear_mat = opw.uniaxial_material.Steel01(osi, fy=f_yield, e0=k_spring, b=r_post)

    # Assign zero length element, # Note: pass actual node and material objects into element
    opw.element.ZeroLength(osi, bot_node, top_node, mat_x=bilinear_mat, r_flag=1)

    # Define the dynamic analysis
    load_tag_dynamic = 1
    pattern_tag_dynamic = 1

    values = list(-1 * rec)  # should be negative
    opy.timeSeries('Path', load_tag_dynamic, '-dt', dt, '-values', *values)
    opy.pattern('UniformExcitation', pattern_tag_dynamic, opw.static.X, '-accel', load_tag_dynamic)

    # set damping based on first eigen mode
    angular_freq = opy.eigen('-fullGenLapack', 1) ** 0.5
    beta_k = 2 * xi / angular_freq
    opw.rayleigh.Rayleigh(osi, alpha_m=0.0, beta_k=beta_k, beta_k_init=0.0, beta_k_comm=0.0)

    # Run the dynamic analysis
    opw.wipe_analysis(osi)

    # Run the dynamic analysis
    opw.algorithm.Newton(osi)
    opw.system.SparseGeneral(osi)
    opw.numberer.RCM(osi)
    opw.constraint.Transformation(osi)
    opw.integrator.Newmark(osi, gamma=0.5, beta=0.25)
    opw.analysis.Transient(osi)

    opw.test_check.EnergyIncr(osi, tol=1.0e-10, max_iter=10)
    analysis_time = (len(values) - 1) * dt
    analysis_dt = 0.001
    outputs = {
        "time": [],
        "rel_disp": [],
        "rel_accel": [],
        "rel_vel": [],
        "force": []
    }

    # access underlying openseespy commands to control analysis
    while opy.getTime() < analysis_time:

        opy.analyze(1, analysis_dt)
        curr_time = opy.getTime()
        outputs["time"].append(curr_time)
        outputs["rel_disp"].append(opy.nodeDisp(top_node.tag, opw.static.X))
        outputs["rel_vel"].append(opy.nodeVel(top_node.tag, opw.static.X))
        outputs["rel_accel"].append(opy.nodeAccel(top_node.tag, opw.static.X))
        opy.reactions()
        outputs["force"].append(-opy.nodeReaction(bot_node.tag, opw.static.X))  # Negative since diff node
    opy.wipe()
    for item in outputs:
        outputs[item] = np.array(outputs[item])


    plt.plot(outputs['time'], outputs['rel_disp'], label='o3seespy')
    periods = np.array([period])

    # Compare closed form elastic solution
    from eqsig import sdof
    resp_u, resp_v, resp_a = sdof.response_series(motion=rec, dt=dt, periods=periods, xi=xi)
    plt.plot(np.arange(len(rec)) * dt, resp_u[0], ls='--', label='Elastic')
    plt.legend()
    plt.show()

.. image:: ./examples/readme_example.png
  :width: 400
  :alt: Output from example

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

* How should youngs modulus be named?



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
