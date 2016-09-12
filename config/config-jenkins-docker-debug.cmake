# Configuration for CI server
# debug build

# main config
set(FLOW_BUILD_TYPE debug)
set(CMAKE_VERBOSE_MAKEFILE on)


# external libraries
set(PETSC_DIR           "/usr/lib/petsc/")
set(PETSC_ARCH          "linux-Debug")

set(BDDCML_ROOT         "/usr/lib/bddcml_blopex/bddcml/Debug")
set(YamlCpp_ROOT_HINT   "/usr/lib/yamlcpp/Debug")
set(Armadillo_ROOT_HINT "/usr/lib/armadillo/Debug")


# additional info
set(USE_PYTHON          "yes")
set(PLATFORM_NAME       "linux_x86_64")