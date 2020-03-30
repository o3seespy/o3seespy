.. o3seespy documentation master file, created by
   sphinx-quickstart on Sat Mar 28 15:34:41 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

o3seespy Documentation
======================

This is the official documentation of the o3seespy package.

o3seespy - Object-oriented OpenSees in Python

The documentation is built off the
OpenSees tcl and OpenSeesPy documentation and has been amended and added to to cover the additional
functionality of the o3seespy package.

The main difference between o3seespy  and OpenSeesPy is that o3seespy provides a pure python layer
to evaluate OpenSees inputs. This improves debugging and allows parameters to be easily accessed within Python.

Contents
========


* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. toctree::

   o3seespy.command
   o3seespy.tools


Features
========

This package provides:

 1. A 'pythonic' version of OpenSees
    - all input parameters to create OpenSees objects are all lowercase key-value arguments
    - all class objects are CamelCase
    - static string variables defined in ALL_CAPS
    - Where possible the exact name used in the original TCL version has been kept
 2. Fully namespaced package allowing full auto-complete e.g. 'o3.uniaxial_material.Steel01(...)'
 3. Replication of underlying object-oriented C++ source code using Python objects.
 4. Type checking of inputs before calling C++ OpenSees code, so that python debugging and errors can be viewed
 5. In code documentation using python docstrings - can view the documentation within your IDE
 6. Additional features for using OpenSees in python:
    - saving and loading data directly from OpenSees into numpy arrays
    - saving and loading data directly from OpenSees into json files
    - Save an entire model as a json file - allows efficient passing of models between servers
 7. All object numbering handled by objects - no need for number tags!
 8. Additional logic checking of optional inputs


Related packages
================

The eco-system of Python packages that relate to this package are outlined in the figure below.

.. image:: https://eng-tools.github.io/static/img/package-space.svg
    :width: 80%
    :align: center
    :alt: Geotechnical Python packages


Example: Dynamic inelastic SDOF analysis
----------------------------------------


.. code:: python

    import numpy as np
    import matplotlib.pyplot as plt

    import o3seespy as o3

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
    osi = o3.OpenSeesInstance(ndm=2, state=0)

    # Establish nodes
    bot_node = o3.node.Node(osi, 0, 0)
    top_node = o3.node.Node(osi, 0, 0)

    # Fix bottom node
    o3.Fix3DOF(osi, top_node, o3.cc.FREE, o3.cc.FIXED, o3.cc.FIXED)
    o3.Fix3DOF(osi, bot_node, o3.cc.FIXED, o3.cc.FIXED, o3.cc.FIXED)
    # Set out-of-plane DOFs to be slaved
    o3.EqualDOF(osi, top_node, bot_node, [o3.cc.Y, o3.cc.ROTZ])

    # nodal mass (weight / g):
    o3.Mass(osi, top_node, mass, 0., 0.)

    # Define material
    k_spring = 4 * np.pi ** 2 * mass / period ** 2
    bilinear_mat = o3.uniaxial_material.Steel01(osi, fy=f_yield, e0=k_spring, b=r_post)

    # Assign zero length element, # Note: pass actual node and material objects into element
    o3.element.ZeroLength(osi, [bot_node, top_node], mats=[bilinear_mat], dirs=[o3.cc.DOF2D_X], r_flag=1)

    # Define the dynamic analysis
    load_tag_dynamic = 1
    pattern_tag_dynamic = 1

    # Define the dynamic analysis
    acc_series = o3.time_series.Path(osi, dt=dt, values=-1 * rec)  # should be negative
    o3.pattern.UniformExcitation(osi, dir=o3.cc.X, accel_series=acc_series)

    # set damping based on first eigen mode
    angular_freq = o3.get_eigen(osi, solver='fullGenLapack', n=1) ** 0.5
    beta_k = 2 * xi / angular_freq
    o3.rayleigh.Rayleigh(osi, alpha_m=0.0, beta_k=beta_k, beta_k_init=0.0, beta_k_comm=0.0)

    # Run the dynamic analysis
    o3.wipe_analysis(osi)

    # Run the dynamic analysis
    o3.algorithm.Newton(osi)
    o3.system.SparseGeneral(osi)
    o3.numberer.RCM(osi)
    o3.constraints.Transformation(osi)
    o3.integrator.Newmark(osi, gamma=0.5, beta=0.25)
    o3.analysis.Transient(osi)

    o3.test_check.EnergyIncr(osi, tol=1.0e-10, max_iter=10)
    analysis_time = (len(rec) - 1) * dt
    analysis_dt = 0.001
    outputs = {
        "time": [],
        "rel_disp": [],
        "rel_accel": [],
        "rel_vel": [],
        "force": []
    }

    while o3.get_time(osi) < analysis_time:
        o3.analyze(osi, 1, analysis_dt)
        curr_time = o3.get_time(osi)
        outputs["time"].append(curr_time)
        outputs["rel_disp"].append(o3.get_node_disp(osi, top_node, o3.cc.X))
        outputs["rel_vel"].append(o3.get_node_vel(osi, top_node, o3.cc.X))
        outputs["rel_accel"].append(o3.get_node_accel(osi, top_node, o3.cc.X))
        o3.gen_reactions(osi)
        outputs["force"].append(-o3.get_node_reaction(osi, bot_node, o3.cc.X))  # Negative since diff node
    o3.wipe(osi)
    for item in outputs:
        outputs[item] = np.array(outputs[item])


    plt.plot(outputs['time'], outputs['rel_disp'], label='o3seespy')
    periods = np.array([period])

    # Compare closed form elastic solution
    from eqsig import sdof
    resp_u, resp_v, resp_a = sdof.response_series(motion=rec, dt=dt, periods=periods, xi=xi)
    plt.plot(np.arange(len(rec)) * dt, resp_u[0], ls='--', label='Elastic')
    plt.legend()
    plt.savefig('readme_example.png')
    plt.show()













