=======
History
=======

3.1.0.XX (2020-XX-XX)
--------------------
* Added commands for applying fixities to list of nodes (e.g. `o3.Fix2DOFMulti`), and for equal DOF command
* Added function for generating a grid of nodes `build_regular_node_mesh`
* Added option for compressing the output of a opy file by applying for loops
* Added `add_fixity_to_dof` to try to apply fixity but not fail if fixity already existing
* Added `friction_models` containing all the friction model objects.
* Added truss element objects
* Fixed issue with BeamOnNonlinearWinklerFoundation where fixities were not applied to base node.
* Fixed issue with getting element node tags when there is only one element
* Fixed issues with Contact elements
* Added more solver algorithms



