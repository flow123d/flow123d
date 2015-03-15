#/bin/bash

# echo "make cmake"
# make cmake > ../1_make_cmake.log 2>&1

echo "make -j 4 all"
make -j 4 all > ../2_make_all.log 2>&1

echo "make -j 4 -C build_tree/unit_tests gtest_mpi_obj"
make -j 4 -C build_tree/unit_tests gtest_mpi_obj > ../3_make_gtest_mpi.log 2>&1

echo "make -C build_tree/unit_tests/fields -k all-test"
make -C build_tree/unit_tests/fields -k all-test > ../4_test_fields.log 2>&1

echo "make -C tests/09_flow_steady_time_dep test-all"
make -C tests/09_flow_steady_time_dep test-all > ../5_flow_09.log 2>&1


echo "find and copy"

branch=$(git branch | sed -n -e 's/^\* \(.*\)/\1/p')
stamp=$(date +"%d_%H-%M-%S")
find . -type f -name *profiler*\*.log -exec mv {} "../logs/${stamp}_${branch}.json" \;
