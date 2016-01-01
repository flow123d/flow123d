/*
 * file_path_test.cpp
 *
 *  Created on: May 24, 2012
 *      Author: jb
 */



#include <flow_gtest.hh>
#include "system/file_path.hh"
#include "system/system.hh"

#include <boost/filesystem.hpp>

using namespace std;

/**
 * Unit test of relative output_dir
 */
TEST(FilePath, output_relative) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    // Following was converted to Warning, that can not be simply tested.

    //EXPECT_DEATH( {FilePath("input/${INPUT}/init.in", FilePath::input_file);},
    //        "Creating FileName object before set_io_dirs is called.");

    FilePath::set_io_dirs("./work_dir/xx", "/main_root", "variant_input", "../output_root");

    { // relative output_dir x relative output substitution
		FilePath fp = FilePath("output.vtk", FilePath::output_file);
		string str_fp = fp;

		// relative output substitution; conversion to string
		EXPECT_EQ(FilePath::get_absolute_working_dir()+"work_dir/output_root/output.vtk", str_fp);

		// conversion to string
		EXPECT_EQ(FilePath::get_absolute_working_dir()+"work_dir/output_root/output.vtk", string(fp));
    }

    // relative output_dir x absolute output substitution, this case is not allowed
    EXPECT_THROW_WHAT( { FilePath(FilePath::get_absolute_working_dir()+"work_dir/xx/output.vtk", FilePath::output_file); },
    		FilePath::ExcAbsOutputPath, "Can not set absolute path" );

}


/**
 * Unit test of absolute output_dir
 */
TEST(FilePath, output_absolute) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    string abs_path = FilePath::get_absolute_working_dir();

    FilePath::set_io_dirs("/work_dir/xx", "/main_root", "variant_input", abs_path+"output_root");

    {
		// relative output substitution; conversion to string
		string str = FilePath("output.vtk", FilePath::output_file);
		EXPECT_EQ(abs_path+"output_root/output.vtk", str);
    }

    {
		// absolute output substitution; conversion to string
		string str = FilePath(abs_path+"output_root/output.vtk", FilePath::output_file);
		EXPECT_EQ(abs_path+"output_root/output.vtk", str);
    }
}


/**
 * Unit test of relative input_dir
 */
TEST(FilePath, input_relative) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    FilePath::set_io_dirs("./work_dir/xx","./main_root", "variant_input", "../output_root");

    { // relative input without substitution; conversion to string
		string str = FilePath("subdir/init.in", FilePath::input_file);
		EXPECT_EQ("./main_root/subdir/init.in", str);
    }

    { // relative input with substitution; conversion to string
		string str = FilePath("subdir/${INPUT}/init.in", FilePath::input_file);
		EXPECT_EQ("./main_root/subdir/variant_input/init.in", str);
    }

    { // absolute input without substitution; conversion to string
		string str = FilePath(FilePath::get_absolute_working_dir()+"subdir/init.in", FilePath::input_file);
		EXPECT_EQ(FilePath::get_absolute_working_dir()+"subdir/init.in", str);
    }

    { // absolute input with substitution; conversion to string
		string str = FilePath(FilePath::get_absolute_working_dir()+"subdir/${INPUT}/init.in", FilePath::input_file);
		EXPECT_EQ(FilePath::get_absolute_working_dir()+"subdir/variant_input/init.in", str);
    }
}


/**
 * Unit test of absolute input_dir
 */
TEST(FilePath, input_absolute) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    string abs_path = FilePath::get_absolute_working_dir();
    FilePath::set_io_dirs("/work_dir/xx", "/main_root", "variant_input", abs_path+"output_root");

    {
        // relative input without substitution; conversion to string
        string str = FilePath("./subdir/init.in", FilePath::input_file);
        EXPECT_EQ("/main_root/./subdir/init.in", str);
    }

    {
        // relative input with substitution; conversion to string
        string str = FilePath("./subdir/${INPUT}/init.in", FilePath::input_file);
        EXPECT_EQ("/main_root/./subdir/variant_input/init.in", str);
    }

    {
        // absolute input without substitution; conversion to string
        string str = FilePath(abs_path+"subdir/init.in", FilePath::input_file);
        EXPECT_EQ(abs_path+"subdir/init.in", str);
    }

    {
        // absolute input with substitution; conversion to string
        string str = FilePath(abs_path+"subdir/${INPUT}/init.in", FilePath::input_file);
        EXPECT_EQ(abs_path+"subdir/variant_input/init.in", str);
    }

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
