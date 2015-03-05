#/bin/bash

make cmake > ../1_make_cmake.log 2>&1
make -j 4 all > ../2_make_all.log 2>&1
make -j 4 -C build_tree/unit_tests gtest_mpi_obj > ../3_make_gtest_mpi.log 2>&1
make -C build_tree/unit_tests/fields -k all-test > ../4_test_fields.log 2>&1
make -C tests/09_flow_steady_time_dep test-all > ../5_flow_09.log 2>&1
