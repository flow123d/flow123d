# Configuration for CI server
# debug build

# main config
set(FLOW_BUILD_TYPE debug)
set(FLOW_CC_FLAGS "-g -O0 -Werror")
set(CMAKE_VERBOSE_MAKEFILE on)


# external libraries
set(PETSC_DIR              /usr/local/petsc-3.8.3/)
set(BDDCML_ROOT            /usr/local/bddcml-2.5.0/bddcml)
set(Armadillo_ROOT_HINT    /usr/local/armadillo-8.3.4)
set(YamlCpp_ROOT_HINT      /usr/local/yamlcpp-0.5.2)


# additional info
set(USE_PYTHON          "yes")
set(PLATFORM_NAME       "linux_x86_64")
