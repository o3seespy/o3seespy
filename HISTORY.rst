=======
History
=======

Pre-release
-----------

3.4.0.3 (2024-03-28)
--------------------
* Updates for PDMY02 and CYCLIQCPSP
* Added `den` to objects in consistent manner referring to mass density

3.2.0.7 (2022-12-13)
--------------------
* Allow saving results with a prefix for file names
* Improved support for variable transient analysis
* Improved mesh partitioning support for tied boundaries
* Fixed issue with modelling imposed ground motions and multi supports

3.2.0.6 (2022-03-03)
--------------------
* Removed double underscore in parameter names
* Added support for element Rayleigh damping

3.2.0.5 (2021-11-29)
--------------------
* Can now define algorithms and test checks without applying them, and can apply/reapply them to the OpenSeesInstance using the method `.reapply()`
* Adding support for database commands
* Minor changes to way inputs are handled for Hysteretic material
* More support for Explicit solvers
* Fixed issue with send command

3.2.0.3 (2021-07-15)
--------------------
* Added option to generate individual fibres
* Added algorithms `NewtonLineSearch` and `Broyden`
* Can now call `get_node_coords` using the node tag if kwarg `node_as_tag=True`
* Fixed bug where parameters could not be passed to `bcast`
* `results` now has support for parallel processors
* Added missing parameters to `algorithm.NewtonLineSearch` object
* Exposed openseespy as `o3.ops` to be consistent with the `ops` convention. `o3.opy` still available for backwards compatibility
* Added `ModalDamping`
* Added `integrator.GimmeMCK` to allow export for different matrices

3.2.0.2 (2021-05-24)
--------------------
* Fixed issue where the `custom_openseespy` package would not be used if installed. Now used as the default if available.
* Can record node and element output to XML file using `NodeToXML`, `NodesToXML`, `ElementToXML` and `ElementsToXML`
* Added uniaxial materials 'PySimple2', `QzSimple2` and `TzSimple2`
* Added support for Results2D output to handle manual node numbering.
* Added method to results where node tags can be rezeroed if not using incremental node numbers
* `node.build_regular_node_mesh()` can handle manual node numbers for 2D and 3D.
* Added option for handling error when applying fixities to list of nodes (e.g. `o3.Fix2DOFMulti`) where node may be `None`.
* New method on `OpenSeesInstance`, `set_log_file()` records logs to a temporary file.
* Added option for `system.apply_mumps_or(<alternative-solver>)` where tries to apply MUMPS solver and if fails will apply an alternative solver.

3.2.0.1 (2021-01-28)
--------------------
* Updated `gen_shallow_foundation_bnwf` command to include shear elements
* Added `mp.partition` command for automatically partitioning the model

3.2.0.0 (2021-01-10)
---------------------
* Removed reference to `o3seespy` in import statements within the package
* Cleaned deprecated files

3.1.0.35 (2020-12-18)
---------------------
* Added support for 4 node elements
* Fixed issues with using the `ManzariDafalias` material
* Fixed issue with setting node mass when defining node, can now set with the `mass=[<nodal-masses>]` parameter

3.1.0.34 (2020-11-03)
---------------------
* Added support for parallel MP node numbering, can now set node tag on creation by passing `tag=<tag-number>`,
* Added support for all OpenSees tags to support parallel MP numbering by having unique numbers for different processor
  ids, simply pass in `mp=True` into `OpenSeesInstance()`.
* Added `o3.domain_change(osi)` to provide the 'domainChange' command.
* Added `integrator.NewmarkExplicit`

3.1.0.33 (2020-10-08)
---------------------
* Added `senscmds` (beta), `parallelcmds`(beta), `fsicmds`(beta)


3.1.0.27 (2020-09-21)
---------------------
* Added `algorithm.ModifiedNewton`
* Added `Parameter`
* Changed install requirements to support latest version of openseespy on Windows and Linux


3.1.0.26 (2020-07-20)
--------------------
* Added commands for applying fixities to list of nodes (e.g. `o3.Fix2DOFMulti`), and for equal DOF command
* Added function for generating a grid of nodes `build_regular_node_mesh`
* Added option for compressing the output of an opy file by applying for loops for repetitive commands
* Added `add_fixity_to_dof` to try to apply fixity but not fail if fixity already existing
* Added `friction_models` containing all the friction model objects.
* Added truss element objects
* Fixed issue with BeamOnNonlinearWinklerFoundation (alpha status) where fixities were not applied to base node.
* Fixed issue with `get_all_ele_node_tags_as_dict` function when there is only one element
* Fixed issues with Contact elements
* Added more solver algorithms



