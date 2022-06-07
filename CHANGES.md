List of all changes in user interface and major internal changes.

***********************************************
# Flow123d version 3.9.0 
(2022-06-06)
Alpha version of 4.x major release, with mostly backward compatible input.

## User interface
* FieldFormula use BParser with some conversion rules to keep FParser syntax, but 
  minor incompatibilities may happen.
* VTK output does not prevent numbering of elements and nodes of the input mesh
* Repporting of HM solver non-convergence

## New features
* Implementation of Native VTK output
* Implementation of GMSH and VTK readers, compatible writers and readers allowing passing outputs as initiali conditions.
* Support for "native" output of presure edge DOFs
* Poroelasticity model with conductivity dependent on the stress tensor.
* Allowed output of a field in several "interpolations", e.g. CellData together with NodeData and NativeData
* FieldFormula can depend on other field within the equation FieldSet
* user fields (test 02/16) 
* Contact non-penetration condition for fractures in the mechanics module.

## Bug fixes
* Poroelasticity bug fixes, extended test suite.
* Year time unit changed to 365.2425 days to be closer to astronomical year.
* velocity reconstruction in zero time step
* partialy fixed sequantial coupling that must work for both steady and unsteady Darcy flow, 
  while the time interpolation of the resulting velocity should be different
* fixed check of the mash compatibility for the boundary FieldFE

## Internals
* BParser used in FieldFormula. Requires SSE2 instructions.
* New assembly algorithm used consistently through the code (without performance regrassion)
* FieldModel introduced.
* FieldModel used in DarcyFlow (pressure head to piezometric head conversion), DG transport Concentration and Heat parameter models, darcy flow to velocity.
* Removed explicit face permutations, replaced by element-nodes premutations with guaranteed face matching
* Optimization of the elements and nodes order for the memory locality.
* removed xprintf
* Boundary Mesh - allow setup of a limited DOF Handler on the boundary
* old asserts completely replaced
* new scheme for building docker images, own build of MPICH with Intel Omnipath support
* Intel and GNU based images, consistent image naming scheme
* VecMPI improved and used in FV transport
* Fixed memory deallocation error in the profiler.
* Use version 4.0.0 of the development images (PETSC 3.17, PERMON 3.17, Armadillo 10.5.2, BDDCML 2.6)
* Fixed default IDs for MSH output format starting now from 1.





***********************************************
# Flow123d version 3.1.0
(2020-12-17)

### New features
* poroelasticity model

### Bug fixes
* ignore unused nodes if these are present in mesh

### Internals
* remove old mixed-hybrid dofhandler
* apply field_fe in equations, mainly passing darcian velocity to transport model, field_fe in reaction term etc.
* improvements in FE classes (FEValues)
* some optimizations, value cache

<!--
Probably even more changes - mirrors current master branch.
-->

***********************************************

# Flow123d version 3.0.4
(2020-06-16)


### Bug fixes
* incorrect values when exporting only pressure_p1 in Darcy Flow model
* check and do not allow duplicit input record keys
* no output of observe data when no equation output fields selected
* when creating region from element ids, check the same dimension of elements
* bug in output when applying DummyDataCache (ended with bad_cast)
* give only warning when mesh contains unused nodes (not belonging to any elements)


# Flow123d version 3.0.3
(2019-09-02)

### New features
* observe data output when running parallel
* removed maximal number of regions limit

### Bug fixes
* fixed several memory leaks (mesh and la)


# Flow123d version 3.0.2
(2019-09-02)

### Bug fixes
* fixed boundary edges indexing in balance to fix fluxes in balance

<!--
# Flow123d version 3.0.1
TODO: do not know what was actually released...

### User interface
* Rename placeholder '${INPUT}' to '$INPUT_DIR$'.
* Improved YAML converter
* FieldElementwise replaced by FieldFE

### Install
* `[windows]` switch from `Docker Toolbox` (using VirtualBox) to `Docker for Windows` (native virtualization)
* `[windows]` update installer, now using native windows installer via NSIS
* `[ linux ]` update installer, can now be installed from console via curl
* update install manual (remove old docs)
* update runtest, better logging support, (verbosity level support)
* docker image hosting is now preferable way to deliver Flow123d


# Flow123d version 3.1.0
* Rename placeholder '${INPUT}' to '$INPUT_DIR$'.
* Improved YAML converter
* FieldElementwise replaced by FieldFE

# Flow123d version 3.0.9
(2019-04-02)

* `[windows]` switch from `Docker Toolbox` (using VirtualBox) to `Docker for Windows` (native virtualization)
* `[windows]` update installer, now using native windows installer via NSIS
* `[ linux ]` update installer, can now be installed from console via curl
* update install manual (remove old docs)
* update runtest, better logging support, (verbosity level support)
* docker image hosting is now preferable way to deliver Flow123d
-->

***********************************************

# Flow123d version 3.0.0
(2017-12-30)

* YAML converter

***********************************************

# Flow123d version 2.2.0
(2017-11-17)

## User interface
* Input YAML file supports including of other YAML files and
  including of tables in CSV files.
* Support for binary and compressed VTK output.
* FieldElementwise support for reading VTK files.
* Summary table for uninitialized fields instead of plenty separate warnings.

### New features
* Solute_AdvectionDiffusion_DG model supports tensor of diffusion.
* Anisotropic automatic choice of DG penalty. Helps for pure diffusion into matrix around
  fractures with high advection.

## Bug fixes
* Fix wrong communication between dimensions when porosities differ
  in heat and solute transport.
* Fix in batch wrappers for Windows. Return to caller after simulation is done.
* Fix of minor error in balance file for DG transport model with sorption.
* Improved options for the linear solver for DG heat and transport.

### Internals:
* Simplification of output classes.
* Simplification in Balance


***********************************************

# Flow123d version 2.1.2
(2017-02-21)

## User interface
* Improved installer script and extended install documentation.

### New features
* Add support for compressed VTK output.
* Better test for initialization of fields.

### Bug fixes
* Fix that runtest hangs when ndiff produced a binary stdout.
* Fix noncompatible P0 algorithm.


### Internals:
* Use python 3 only.
* Own build of mpich preventing error on valgrind use.
* Cleanup of input headers.


***********************************************

# Flow123d version 2.1.0
(2016-12-16)

## User interface

* introduction of tutorials (tests/05_tutorials)
* Docker is used for running Flow123d in platform independent environment.
* Change meaning (and names) of parameters of Sorptions. Closer to common usage.
* Updated documentation.
* Improved format of generated input reference fro both Latex and HTML.



### New features
* FieldTimeFunction - field constant in space and interpolated in time from a table of values.
* Unit conversion for fields.
* Support for binary VTK output.
* Checks for field limits.

### Bug fixes
* Fix unsteady water flow without prescribed timestep.
* Fix configuration of the linear solver as part of the nonlinear solver.
* Fix problem with YAML tag in autoconversion of Abstract.
* Fix search of observe points.
* Fix curious bug in selection of output times specific to individual equations.
* Fix unsteady MH solver.
* Fix observe output for no observation fields.
* Fix some destructors.
* Fix solute balance in zero time.
* Fix init_piezo_head

### Internals:

* reorganization of the integration tests
* use Docker to build and test Flow123d
* replace boost::shared_ptr by std::shared_ptr
* Simplified finishing of input types.

***********************************************

# Flow123d version 2.0.0
(2016-08-31)

### User interface
* Adopt YAML (http://www.yaml.org/start.html) as the format for the principal input file.
  The CON format (derivate of JSON) is still accepted, but obsolete. Only pure JSON format will be supported in the future.

* A tool is provided for conversion from the CON with the structure of the 1.8.2 version to the YAML
  format with the new structure. The tool can also convert YAML files of the version 2.0.0_rc.

  Usage:    ```bin/input_convert.sh path/to/old_file.con```
            ```bin/input_convert.sh path/to/old_file.yaml```  

* Vector valued fields are replaced by "multifields". This allows independent setting for individual components.
  E.g. ``` init_conc = [ 0, 1, {TYPE=FieldFormula, value="x*y"} ]```

* Automatic conversion from Record to Array is implemented as a transpose, e.g.
  ```{ a: [1,2], b: [3,4] }```

  may be converted to

  ```[ {a:1, b:3}, {a:2, b:4} ]. ```

  Useful for multifields.

* Field descriptors in the 'input_fields' list need not to form increasing time series, this is required only
  for the sequence of field descriptors of the single field. E.g. this is valid:
    ```
    input_fields : [
       { time:0.0, region:"ALL", conductivity:1 },
       { time:1.0, region:"ALL", conductivity:2 },
       { time:0.0, region:"ALL", init_pressure: 0.5}        
    ]
    ```
* Removed support for old BCD files.

* Introduction of observation points.

* Independent output times for the balance output and the field output.
  Support for more complex scheme of output times.

* Changes in the structure of the input file:

    * Add obligatory key: `/flow123d_version` with a string marking the version of the file format

    * All boundary fluxes on input are treated as positive if they are increasing the mass balance.
      This is contradictory to the convention used in PDE, but is consistent with treatment of the volume sources.
      Same convention is used in balance output files. The conversion script takes care of
      sign change for constant and formula input fields, but produce an invalid input
      in the case of `FieldElementwise` or `FieldPython` in order to mark these fields that must be
      fixed manually by the user.

    * Combined Neumann + Robin boundary condition. The separate BC types 'neumann' and 'robin' for
      the Darcy flow and the DG transport are replaced by the single BC type 'total_flux'.
      Moreover, DG transport supports BC type 'diffusive_flux' which applies only to
      the diffusion-dispersion boundary flux.


    * Rename couplings to: 'Coupling_OperatorSplitting', 'Coupling_Sequential'
    
    * 'Coupling_Sequantial' have keys 'flow_equation', 'solute_equation', 'heat_equation' instead of
    'primary_equation' and 'secondary_equation'.
    
    * Rename equations to: 'Flow_Darcy_MH', 'Flow_Richards_LMH', Solute_AdvectionDiffusion_DG', 'Solute_Advection_FV', 'Heat_AdvectionDiffusion_DG',


     * Move setting of the Solute_Advection_FV and Solute_AdvectionDiffusion_DG under the
       key 'transport' of the Coupling_OperatorSplitting. This is necessary to allow coupling
       of reactions with the diffusive transport.


     * Modification of the mesh setting. Regions and region sets are now treated equally, the term "region"
       now means the set of elementary regions, i.e. what was the region set up to now. Regions are imported
       from the mesh (named physical group in GMSH), further regions may be created by operations in the
       list: /problem/mesh/regions; which now support operations from both previous lists:
           /problem/mesh/regions and /problem/mesh/sets     
    
     * Remove the key 'r_set' from the field descriptors. Use the key 'region' instead. Moreover,
       list of region labels may be used.
    
     * Flow_Darcy_MH dynamicaly switch between steady and unsteady behavior according to the value of the
     'storativity' field. Zero value (default) force the steady model. Similarly for the Flow_Richards_LMH
     however both fields 'storativity' and 'water_content_saturated' must be zero (default).
    
     * Rename the key 'solver' of the flow equations to 'linear_solver'  and move it under the
      new key 'nonlinear_solver'.
    
     * Separation of output_stream (only for main equations) with specification of resulting format and setting common to
       equation outputs usign the stream. The equation outputs specify "what to output", namely 'output_fields' and 'observe_fields'.
    
     * Rename 'BOUNDARY' and 'IMPLICIT BOUNDARY' region sets to '.BOUNDARY' and '.IMPLICIT_BOUNDARY' respectively.
    
     * Add HTML format of the generated documentation to the structure of the input file.




### New features
* Experimental support of Richards equation (model 'Flow_Richards_LMH'). Can not be used with transport processes yet.
* DG - Reactions operator splitting.
* 'total_flux' combined Neumann and Robin BC for the flow equations.
* 'seepage_face' and 'river' boundary conditions for the flow equations
* 'total_flux' combined Neumann and Robin BC and 'diffusive_flux' the Solute_AdvectionDiffusion_DG equation.
* Coupling of the Solute_AdvectionDiffusion_DG model with the general reaction term via operator splitting.
* Using stable Pade aproximant for Decays and FirstOrderReaction.
* Make whole reaction term models unconditionaly stable (no additional time step constraints).
* Use conservative formulation in Solute_Advection_FV, allow time dependent 'porosity' and 'cross_section' fields.
* Unified messaging and logging.
* Use PETSc 3.6



### Bug fixes
* Fix bug in treatment of the paths given by -i and -o command line arguments.
* Fix calculation of the dispersivity for small velocities in Solute_AdvectionDiffusion_DG.
* Allow arbitrary sources in Solute_Advection_FV, set CFL condition appropriately.
* Fix bug in release 2.0.0_rc in Richards_LMH when using 3D elements.

### Internals:
* Introduce generic input types, allow simplification of the documentation.
* Introduce attributes for input types and record keys, allow passing additional information to exteranal
  applications, e.g. GeoMop.
* Introduce input type Tuple.
* YAML input reader.
* Treat Input::Abstract as an interface. Allow Record derived from more abstracts (interfaces).
  Input::Abstract have no keys and is not derived from the Record anymore.
* Adopt the linear algebra classes from the deal.ii project.
* Memory profiler.
* Consistent usage of the ReferenceElement class, to get numbering and orientation of various geometry entities.
* Display code coverage of unit tests on Jenkins build server interface.
* Input::Factory mechanism is used for constructing objects according to the input. Reduce dependencies of sources.
* Executable wrappers dealing with compatibility of the shared libraries.
* General default values (JSON format).
* Stream based messaging and logging.
* ASSERTS with simple output of involved variables.



************************************************

# Flow123d version 1.8.2
(2015-03-15)

### User interface

* Key DarcyMHOutput/balance_output replaced by DarcyFlowMH/balance.
* Key TransportOperatorSplitting/mass_balance replaced by TransportOperatorSplitting/balance.
* Key SoluteTransport_DG/mass_balance replaced by SoluteTransport_DG/blance.
* Key HeatTransfer_DG/mass_balance replaced by HeatTransfer_DG/balance.
* Unified record Balance for configuration of balance output.
* New key DarcyFlowMH/gravity to set gravity vector.
* All output selections support fields 'region_id' and 'subdomain'.
* TransportOperatorSplitting/substances and TransportDG/substances.
  are arrays of Substance records. Substance record contains substance name and its molar mass.
  Molar mass is then used in the ReactionTerm.
* Key 'molar_mass' removed from Sorption model records.
* Reaction term model LinearReactions is replaced by specialized models FirstOrderReaction and RadioactiveDacay
  Both new models have records with structure similar to the LinearReactions by specialized to
  the particular physical context. In particular RadioactiveDecayProduct allows specification of the decay energy
  with in future will be used as a source in the heat equation.
* Reaction term model PadeAproximant is now variant of 'ode_solver' used by some reaction models.
* Alternative 'ode_solver' is LinearODEAnalytic that was used as default by previous versions.
* TimeGovernor is now autoconstructable from the key 'max_dt'.
* Meaning of TimeGovernor/init_dt has changed. Newly it is suggestion of the very first time step.

### New features
* United implementation of balance calculations, fast even for explicit solvers.
* Balance output readable by table processors.
* SoluteTransport_DG: changes in interpretation of boundary conditions (see manual for current meaning).
* Faster BIHTree for search of non-comaptible neighborings.
* Continuous integration of all install packages.
* Windows installer.
* Automatic creation of the output directory.
* Faster dual porosity solver.


### Bug fixes
* Memory leak in output classes.
* All models now adjust its time steps to match changes of input fields and output times.
  Several bugs fixed in all models. Added dedicated tests.
* Wrong profiler times for the linear solver.
* Infinite loop during processing of the --JSON_machine command line argument.
* Wrong treatment of possible failure of nonlinear solver in sorption. Missing check of input values.
* Binary for 32-bit Windows did not run. (CYGWIN bug)
* Resolved instability of inflow/outflow detection in DG solver.
* Wrong line numbers in some error messages produced by JSON parser.
* Missing error messages on Windows system.
* Bugs in non-compatible 2D-1D assembly in DarcyFlowMH.
* Radioactive decay should depend on the molar weight (conserve amount of nuclei not mass).
* PadeAproximant did not work in parallel.
* Wrong mass balance with non-zero reaction term.
* Memory leak in LinSys_PETSC.
* Fails to read mesh with region ID repeating for regions of different dimension.







***********************************************

# Flow123d version 1.8.1
(2014-05-14)

### Bug fixes

* fix missing error messages under Windows
* fix a memory leak in LinSys_PETSC::gather_solution
* fix weights in BDDC for Darcy flow
* fix wrong behavior when wrong region name is specified



### Other changes
* minor improvements of manual
* add Python support under Windows
* add BDDC support under Windows
* use BDDCML 2.4
* fix and unit tests for BIHTree class
* improved test 23 - dual porosity
* updated petsc_options_help


************************************************

# Flow123d version 1.8.0
(2014-03-31)

### User interface:
* every equation can output any input, computed, or post-processed field,
  new keys: 'output_stream', 'output_fields'
* input fields given in one list under key 'input_fields' instead of separated 'bulk_data' and 'bc_data' lists   
* improved mass balance files (yet not correct with reaction term)
* output of discontinuous linear finite elements
* enhancements in flow123d.sh script; changed syntax, support for MetaCentrum clusters
* support for BDDCML parallel solver in Darcy flow model
* relative references in the input file produce warning reporting used absolute path
* input interface warnings for unknown keys in the input file
* warning for usage of field default values
* Flow123d tries to keep IDs of nodes and elements in output into GMSH format same as in the input mesh;
  but it is not guaranteed to future and/or for parallel case

### Physics:
* transition coefficients in coupling of different dimensions are computed from other parameters,
  'sigma' parameter is additional scaling coefficient with default value 1.0
* concentration sources in Transport DG
* general reaction term for operator splitting; can combine dual porosity, adsorption and reactions
* limited solubility model in adsorptions
* TransportDG supports second and third order finite elements
* TransportDG supports Dirichlet, Neumman, and Robin boundary conditions
* Darcy flow supports non-compatible 2d-1d meshes
* support for heat transfer simulations (based on TransportDG); can not be combined with solute transport yet
* ConvectionTransport adapt its time step to match change times of  the input fields


### Internals:
* transition to git, github
* transition to ci2 test server; YAML description of tests
* new output subsystem
* Fields can read from input independently and can hold short time history
* optimization of mesh reading
* optimization of Schur complement computation
* optimization of TransportDG assembly
* fast interpolation method for nonlinear adsorptions
* output of the input type tree in machine readable JSON format
* use C++11; g++ 4.7 and newer
* transition to petsc 3.4.3
* run_test.sh script use flow123d.sh consistently, can be used also on parallel systems
* general classes for interpolation tables (not used yet)
* AdHocAbstractRecord - abstract record with different kind of documentation
* automatic installation of BDDCML library
* reimplemented BIHTree, faster and more robust
* transition to gtest 1.7.0
* automatic night packages for both Linux and Windows
* replace EqDataBase by more general class FieldSet
