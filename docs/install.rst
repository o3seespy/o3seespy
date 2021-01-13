Installation
============

To use o3seespy you need to run a python compiled version of OpenSees.
The easiest way to do this is to use the official version: `openseespy`, which is automatically installed
when you install o3seespy with pip.

.. code:: bash

    pip install o3seespy

Running o3seespy on MacOS
-------------------------

Some versions of openseespy are not supported on MacOS. To use o3seespy you have three options:

 1. Use openseespy==3.1.5.12, `pip install openseespy==3.1.5.12`
 2. Compile OpenSees to Python yourself using the `compiling script
<https://github.com/eng-tools/OpenSees/blob/master/MAKES/Makefile.def.MacOS10.15-python>`_
 3. Contact the developers of o3seespy for a compiled version.

If you are using options 2 or 3, you need to follow the steps below for using a self-compiled version of openseespy.


Using a self-compiled version of openseespy with o3seespy
---------------------------------------------------------

o3seespy has a provision to run a module called `custom_openseespy` instead of `openseespy` if `custom_openseespy` is installed.
To do this you need to:

1. compile openseespy to create the Python shared library ('opensees.so' on Mac and Linux, 'opensees.pyd' for Windows),
2. you need to place the shared library file into a folder called 'custom_openseespy'
3. `copy the setup.py file from here<https://github.com/o3seespy/o3seespy/blob/master/docs/setup_py_file_for_custom_openseespy.txt>`_
 make sure you name the file 'setup.py' and put it in the 'custom_openseespy' folder.
4. In your terminal (start your virtual environment if you are using one), and change the directory to the folder
that contains the 'custom_openseespy' folder,
5. Run the following command `pip install -e custom_openseepy`.
6. Now you should be able to run o3seespy and it will use the `custom_openseespy` version.
