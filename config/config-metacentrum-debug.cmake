# Configuration for metacentrum tarkil
# debug build

# main config
set(FLOW_BUILD_TYPE debug)
set(CMAKE_VERBOSE_MAKEFILE on)
set(LIB_ROOT             "/auto/praha1/jan-hybs/libs")


# external libraries
set(PETSC_DIR           "${LIB_ROOT}/petsc/")
set(PETSC_ARCH          "linux-Debug")

set(BDDCML_ROOT         "${LIB_ROOT}/bddcml_blopex/bddcml/Debug")
set(YamlCpp_ROOT_HINT   "${LIB_ROOT}/yamlcpp/Debug")
set(Armadillo_ROOT_HINT "${LIB_ROOT}/armadillo/Debug")


# additional info
set(USE_PYTHON          "yes")
set(PLATFORM_NAME       "linux_x86_64")


# need to load following modules in order to build flow123d
# 
# module add cmake-2.8.12
# module add mpich2
# module load boost-1.56-gcc
# module load python-2.7.10-gcc