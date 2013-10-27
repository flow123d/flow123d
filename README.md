Flow123d
========

Flow123d is a simulator of underground water flow and transport in fractured
porous media. Novelty of this software is support of computations on complex
meshes consisting of simplicial elements of different dimensions. Therefore
we can combine continuum models and discrete fracture network models.

Current version includes mixed-hybrid solver for steady and unsteady Darcy's
flow, finite volume model for transport of several substances. Using operator
splitting we support models for various local processes including dual
porosity, sorption, decays and reactions.

All computations can be run in parallel using PETSc and MPI. However, the
mesh data and output operations are still sequential. Program supports output
into VTK format for visualization and postprocessing in Paraview. 

License
-------

The source code of Flow123d is licensed under GPL3 license. For details look
at files LICENSE and GPL3.

Build Flow123
-------------

### Prerequisities

If you are running Windows, you have to install 'cygwin' for emulation of
POSIX unix environment. Then all work has to be done in the directories under
cygwin, e.g "C:\cygwin\home\user\".

For build you will need some developement software packages. On Windows use
cygwin to install them. On Linux use package manager of your distribution
('apt-get' for Ubuntu and Debian, 'yum' for RedHat and Centos)

Requested packages are: 

* gcc, gcc4-g++, gcc4-fortran     C/C++ compiler and Fortran77 compiler for compilation of BLAS library.
* python, perl, cmake (>2.8)      Scripting languages 
* make, cmake (>2.8), patch       Building tools
* libboost                        General purpose C++ library 

> Flow123d can install Boost automatically during configuration, but it may
> take long.
  
Namely you need development version of the sub-libraries "Program Options" and
"Serialize" in KUbuntu these are in separate packages:

* libboost-program-options-dev,
* libboost-serialize-dev

Under Cygwin you also need to install:  

* openssh (for MPICH)
* some editor that can save Unix like text files (notepad++, vim under cygwin) 

optionaly you may be interested in:

* doxygen, graphviz     source generated documentation and its support tool for diagrams 

Other libraries that can be possibly useful for developers can be found on address:

http://dev.nti.tul.cz/trac/flow123d/wiki/Developement
http://dev.nti.tul.cz/trac/flow123d/wiki/Software

### Install PETSc Library

Flow versions 1.7.x depends on the PETSC library 3.2.0-xx.
You can download this version from:

http://www.mcs.anl.gov/petsc/petsc-as/documentation/installation.html

Assume that you unpack the installation tree to the directory:

    $ cd /home/jb/local/petsc

Change to this directory and set this as your PETSC directory:

    $ export PETSC_DIR=`pwd`

For developelemt you will need at least debuging and production build of the
library. First set a name for the debugging configuration:

    $ export PETSC_ARCH=linux-gcc-dbg

And run the configuration script, for example with following options:

    $ ./config/configure.py --with-debugging=1 --CFLAGS-O=-g --FFLAGS-O=-g \
      --download-mpich=yes --download-parmetis=yes --download-f-blas-lapack=1

This also automagically install BLAS, Lapack, MPICH, and ParMetis so it takes
a while, it can be about 15 min. If everything is OK, you obtain table with
used compilers and libraries. Finally compile PETSC with this configuration:

    $ make all

To test the compilation run:

    $ make test

To obtain PETSC configuration for the production version you can use e.g.

    $ export PETSC_ARCH=linux-gcc-dbg
    $./config/configure.py --with-debugging=0 --CFLAGS-O=-O3 --FFLAGS-O=-O3 \
       --download-mpich=yes --download-parmetis=yes --download-f-blas-lapack=1
    $ make all
    $ make test

> ** Important Notes: **
> 
> * For some reasons if you let PETSc to download and install its own MPICH it
>   overrides your optimization flags for compiler. Workaround is to edit
>   file ${PETSC_DIR}/${PETSC_ARCH}/conf/petscvariables and modify variables
>   XXXFLAGS_O back to the values you wish.
> * For production builds the PETSc configuration should use system wide
>   MPI implementation.
> * You have to compile PETSc with Parmetis support.
> * Configurations mentioned above are minimalistic. 
>   Next we describe several additional configure options which can be useful.

#### Alternatives and Troubleshooting

* Windows: If you use a shell script for PETSC configuration under cygwin,
always check if you use UNIX line ends. It can be specified in the notepad
of Windows 7.
  
* For some Windows versions it may be necessary to set in configuration of PETSC
    --with-timer=nt 

* You may get strange errors during configuration of PETSc, like 

    C:\cygwin\bin\python.exe: *** unable to remap C:\cygwin\bin\cygssl.dll to same address as parent(0xDF0000) != 0xE00000

or other errors usually related to DLL conflicts.
(see http://www.tishler.net/jason/software/rebase/rebase-2.4.2.README)
To fix DLL libraries you should perform:

 1. shutdown all Cygwin processes and services
 2. start 'ash' (do not use bash or rxvt)
 3. execute '/usr/bin/rebaseall' (in the ash window)

Possible problem with 'rebase':
    /usr/lib/cygicudata.dll: skipped because nonexist
    .
    .
    .
    FixImage (/usr/x86_64-w64-mingw32/sys-root/mingw/bin/libgcc_s_sjlj-1.dll) failed with last error = 13

   Solution (ATTENTION, depends on Cygwin version): 
    add following to line 110 in the rebase script:
    -e '/\/sys-root\/mingw\/bin/d'


* For some Windows versions it may be necessary to set in configuration of PETSC
  (TODO: symptoms of this issue ?)
   
    --with-timer=nt 

* By default PETSC will create dynamically linked libraries, which can be shared be more applications. But on some systems 
  (in particular we have problems under Windows) this doesn't work, so one is forced to turn off dynamic linking by:

    --with-shared=0

* If you want only serial version of PETSc (and Flow123d)
  add --with-mpi=0 to the configure command line.

* You can have several PETSC configuration using different PETSC_ARCH.
  
* PETSC use BLAS and few LAPACK functions for its local vector and matrix
  operations. The speed of BLAS and LAPACK have dramatic impact on the overall
  performance. There is a sophisticated implementation of BLAS called ATLAS.
  ATLAS performs extensive set of performance tests on your hardware then make
  an optimized implementation of  BLAS code for you. According to our
  measurements the Flow123d is about two times faster with ATLAS compared to
  usual --download-f-blas-lapack (on x86 architecture and usin GCC).
   
  In order to use ATLAS, download it from ... and follow their instructions.
  The key point is that you have to turn off the CPU throttling. To this end
  install 'cpufreq-set' or `cpu-freq-selector` and use it to set your processor
  to maximal performance:
  
    cpufreq-set -c 0 -g performance
    cpufreq-set -c 1 -g performance

   ... this way I have set performance mode for both cores of my Core2Duo.

   Then you need not to specify any special options, just run default configuration and make. 
   
   Unfortunately, there is one experimental preconditioner in PETSC (PCASA) which use a QR decomposition Lapack function, that is not
   part of ATLAS. Although it is possible to combine ATLAS with full LAPACK from Netlib, we rather provide an empty QR decomposition function
   as a part of Flow123d sources.
   See. HAVE_ATTLAS_ONLY_LAPACK in ./makefile.in

* PETSC provides interface to many useful packages. You can install them 
  adding further configure options:

  --download-superlu=yes         # parallel direct solver
  --download-hypre=yes           # Boomer algebraic multigrid preconditioner, many preconditioners
  --download-spools=yes          # parallel direc solver
  --download-blacs=ifneeded      # needed by MUMPS
  --download-scalapack=ifneeded  # needed by MUMPS
  --download-mumps=yes           # parallel direct solver
  --download-umfpack=yes         # MATLAB solver

For further information about use of these packages see:

http://www.mcs.anl.gov/petsc/petsc-2/documentation/linearsolvertable.html

http://www.mcs.anl.gov/petsc/petsc-as/snapshots/petsc-current/docs/manualpages/PC/PCFactorSetMatSolverPackage.html#PCFactorSetMatSolverPackage
http://www.mcs.anl.gov/petsc/petsc-as/snapshots/petsc-current/docs/manualpages/Mat/MAT_SOLVER_SPOOLES.html#MAT_SOLVER_SPOOLES
http://www.mcs.anl.gov/petsc/petsc-as/snapshots/petsc-current/docs/manualpages/Mat/MAT_SOLVER_MUMPS.html#MAT_SOLVER_MUMPS
http://www.mcs.anl.gov/petsc/petsc-as/snapshots/petsc-current/docs/manualpages/Mat/MAT_SOLVER_SUPERLU_DIST.html
http://www.mcs.anl.gov/petsc/petsc-as/snapshots/petsc-current/docs/manualpages/Mat/MAT_SOLVER_UMFPACK.html

http://www.mcs.anl.gov/petsc/petsc-as/snapshots/petsc-dev/docs/manualpages/PC/PCHYPRE.html

### Build Step 2 - Compile Flow123

Copy file  makefile.in.cmake.template to makefile.in.cmake:

    $ cp makefile.in.cmake.template makefile.in.cmake

Edit file makefile.in.cmake, set PETSC_DIR and PETSC_ARCH variables.

You can specify type of build:

    $ set(CMAKE_BUILD_TYPE debug)

or 

    $ set(CMAKE_BUILD_TYPE release)

or you can directly set flags for C and C++ compiler:
      
    $ set(CC_FLAGS "-O3 -DNODEBUG -pg ")

Then run the compilation by:

    $ make all

This should run configuration and then build process. If the configuration
fails and you need to correct your makefile.in.cmake or other system setting
you have to cleanup all files generated by unsuccessful cmake configuration by:

    $ make clean-all

Try this every if your build doesn't work and you don't know why.

> **Petsc Detection Problem:**
>
> CMake try to detect type of your PETSc installation and then test it
> by compiling and running simple program. However, this can fail if the 
> program has to be started under 'mpiexec'. In such a case, please, set:
>
>    set (PETSC_EXECUTABLE_RUNS YES)
>
> in your makefile.in.cmake file, and perform: `make clean-all; make all`

--------------------------------------------------------------------------------

For further information about program usage see documentation in "doc/" in
particular reference manual "doc/flow_doc". 

