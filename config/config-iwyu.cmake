###############################################3
# Configuration for running IWYU (include what you use) tool.
# Assumes IWYU installed manually and accessible from docker.

# Clang must be used for IWYU
set( CMAKE_CXX_COMPILER "clang++-3.8")

set(FLOW_BUILD_TYPE debug)
set(CMAKE_VERBOSE_MAKEFILE on)
set(USE_PYTHON "yes")
set(PLATFORM_NAME "linux_x86_64")



# --------------------------------------

set(PETSC_DIR           "/usr/lib/petsc/")
set(PETSC_ARCH          "linux-Debug")
set(BDDCML_ROOT         "/usr/lib/bddcml_blopex/bddcml/Debug")
set(YamlCpp_ROOT_HINT   "/usr/lib/yamlcpp/Debug")

# Path to the installation of IWYU
# see also 'installation steps' script config/iwyu.sh

set(CMAKE_CXX_INCLUDE_WHAT_YOU_USE "/usr/local/bin/include-what-you-use")
