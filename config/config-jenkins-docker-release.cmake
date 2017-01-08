# Configuration for CI server
# release build

# main config
set(FLOW_BUILD_TYPE release)
set(CMAKE_VERBOSE_MAKEFILE on)


# external libraries
set(PETSC_DIR           "/usr/lib/petsc/")
set(PETSC_ARCH          "linux-Release")

set(BDDCML_ROOT         "/usr/lib/bddcml_blopex/bddcml/Release")
set(YamlCpp_ROOT_HINT   "/usr/lib/yamlcpp/Release")
set(Armadillo_ROOT_HINT "/usr/lib/armadillo/Release")


# additional info
set(USE_PYTHON          "yes")
set(PLATFORM_NAME       "linux_x86_64")
