List of all changes in user interface and major internal changes.

***********************************************

#Flow123d version 1.8.1
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

#Flow123d version 1.8.0 
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