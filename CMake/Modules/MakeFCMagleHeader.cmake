# Set desired compilers through standard CMake variables
# CMAKE_Fortran_COMPILER
# CMAKE_C_COMPILER
#  
#
project(FortranMaglingHeader C CXX Fortran)

include(FortranCInterface)
FortranCInterface_HEADER(FC_Magle.h MACRO_NAMESPACE "FC_SYMBOL") 