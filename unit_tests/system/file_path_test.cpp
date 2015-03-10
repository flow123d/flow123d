/*
 * file_path_test.cpp
 *
 *  Created on: May 24, 2012
 *      Author: jb
 */



#include <flow_gtest.hh>
#include "system/file_path.hh"

#include <boost/filesystem.hpp>

using namespace std;

TEST(FilePath, output_relative) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    // Following was converted to Warning, that can not be simply tested.

    //EXPECT_DEATH( {FilePath("input/${INPUT}/init.in", FilePath::input_file);},
    //        "Creating FileName object before set_io_dirs is called.");

    FilePath::set_io_dirs("./work_dir/xx","/main_root", "variant_input", "../output_root");


    FilePath fp = FilePath("output.vtk", FilePath::output_file);
    string str_fp = fp;

    // relative output substitution; conversion to string
    EXPECT_EQ(FilePath::get_absolute_working_dir()+"/work_dir/output_root/output.vtk", str_fp);

    // conversion to string
    EXPECT_EQ(FilePath::get_absolute_working_dir()+"/work_dir/output_root/output.vtk", string(fp));

}

TEST(FilePath, output_absolute) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    string abs_path = FilePath::get_absolute_working_dir();

    FilePath::set_io_dirs("/work_dir/xx","/main_root", "variant_input", abs_path+"/output_root");

    // relative output substitution; conversion to string
    string str = FilePath("output.vtk", FilePath::output_file);
    EXPECT_EQ(abs_path+"/output_root/output.vtk", str);
}

TEST(FilePath, input) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    FilePath::set_io_dirs("./work_dir/xx","/main_root", "variant_input", "../output_root");

    // relative output substitution; conversion to string
    string str = FilePath("subdir/${INPUT}/init.in", FilePath::input_file);
    EXPECT_EQ("/main_root/subdir/variant_input/init.in", str);
}


TEST(FilePath, create_output_dir) {
    //::testing::FLAGS_gtest_death_test_style = "threadsafe";
    boost::filesystem::remove_all("./work_dir");
    FilePath::set_io_dirs("./work_dir","", "my_input", "my_output");
    EXPECT_TRUE(boost::filesystem::is_directory("./work_dir/my_output"));
    FilePath fp("subdir/some_file.xyz", FilePath::output_file);
    fp.create_output_dir();
    EXPECT_TRUE(boost::filesystem::is_directory("./work_dir/my_output/subdir"));
    boost::filesystem::remove_all("./work_dir");
}
