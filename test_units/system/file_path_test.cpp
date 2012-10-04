/*
 * file_path_test.cpp
 *
 *  Created on: May 24, 2012
 *      Author: jb
 */



#include <gtest_throw_what.hh>
#include "file_path.hh"

using namespace std;

TEST(FilePath, output_relative) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    EXPECT_DEATH( {FilePath("input/${INPUT}/init.in", FilePath::input_file);},
            "Creating FileName object before set_io_dirs is called.");

    FilePath::set_io_dirs("/main_root", "variant_input", "../output_root");

    FilePath fp = FilePath("output.vtk", FilePath::output_file);
    string str_fp = fp;

    // relative output substitution; conversion to string
    EXPECT_EQ("/main_root/../output_root/output.vtk", str_fp);

    // conversion to string
    EXPECT_EQ("/main_root/../output_root/output.vtk", string(fp));

}

TEST(FilePath, output_absolute) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    FilePath::set_io_dirs("/main_root", "variant_input", "/output_root");

    // relative output substitution; conversion to string
    string str = FilePath("output.vtk", FilePath::output_file);
    EXPECT_EQ("/output_root/output.vtk", str);
}

TEST(FilePath, input) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    FilePath::set_io_dirs("/main_root", "variant_input", "");

    // relative output substitution; conversion to string
    string str = FilePath("subdir/${INPUT}/init.in", FilePath::input_file);
    EXPECT_EQ("/main_root/subdir/variant_input/init.in", str);
}
