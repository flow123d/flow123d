# Configuration for CI server
# release build

# main config
set(FLOW_BUILD_TYPE release)
#set(FLOW_CC_FLAGS "-O3 -DNDEBUG -Werror")
set(CMAKE_VERBOSE_MAKEFILE on)


# external libraries
set(PETSC_DIR              /usr/local/petsc_v3.18.6)
set(BDDCML_ROOT            /usr/local/bddcml_2.6)
set(PERMON_ROOT            /usr/local/permon_3.18.0)
set(Armadillo_ROOT_HINT    /usr/local/armadillo_12.2.0)
set(YamlCpp_ROOT_HINT      /usr/local/yamlcpp_0e6e28d)
set(EIGEN_ROOT             /usr/local/eigen-3.4.0)


# additional info
set(USE_PYTHON          "yes")
set(PLATFORM_NAME       "linux_x86_64")
