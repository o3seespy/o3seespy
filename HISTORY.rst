=======
History
=======

3.2.0.1 (2021-01-28)
---------------------
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



