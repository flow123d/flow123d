# Configuration for CI server
# debug build

# main config
set(FLOW_BUILD_TYPE coverage)
set(CMAKE_VERBOSE_MAKEFILE on)


# external libraries
set(PETSC_DIR              /usr/local/petsc_3.8.3/)
set(BDDCML_ROOT            /usr/local/bddcml_2.6)
set(Armadillo_ROOT_HINT    /usr/local/armadillo_8.3.4)
set(YamlCpp_ROOT_HINT      /usr/local/yamlcpp_0.6.3)


# additional info
set(USE_PYTHON          "yes")
set(PLATFORM_NAME       "linux_x86_64")
