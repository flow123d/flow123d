# Configuration for CI server
# debug build

# main config
set(FLOW_BUILD_TYPE debug)
#set(FLOW_CC_FLAGS "-g -O0 -Werror")
set(CMAKE_VERBOSE_MAKEFILE on)


# external libraries
set(USE_PYTHON          "yes")
set(PETSC_DIR              /usr/local/petsc_3.25.0)
set(BDDCML_ROOT            /usr/local/bddcml_2.6)
set(PERMON_ROOT            /usr/local/permon_3.25.0)
set(Armadillo_ROOT_HINT    /usr/local/armadillo_15.2.4)
set(YamlCpp_ROOT_HINT      /usr/local/yamlcpp_0e6e28d)
set(EIGEN_ROOT             /usr/local/eigen-3.4.0)

# additional info
set(USE_PYTHON          "yes")
set(PLATFORM_NAME       "linux_x86_64")
