/*
 * file_path_test.cpp
 *
 *  Created on: May 24, 2012
 *      Author: jb
 */



#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>
#include "system/file_path.hh"
#include "system/system.hh"

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

using namespace std;


/*
 * Conversion function for correct result of tests under Cygwin.
 * Replace double backslash with single slash.
 */
string convert_path(string path) {
#ifdef FLOW123D_HAVE_CYGWIN
    boost::replace_all(path, "\\", "/");
#endif // FLOW123D_HAVE_CYGWIN
	return path;
}


TEST(BoostFileSystem, base_methods) {
	// current_path
	std::cout << "Current path: " << boost::filesystem::current_path() << std::endl;

	// simple path
	{
		boost::filesystem::path path = boost::filesystem::path(".") / "x";
		string path_str = convert_path( path.string() );
		EXPECT_EQ("./x", path_str);
	}

	// path contains actual dir ('.'), format: "/./x/./y"
	{
		boost::filesystem::path path = boost::filesystem::path(".") / "x" / "." / "y";
		string path_str = convert_path( path.string() );
		EXPECT_EQ("./x/./y", path_str);
	}

	// path contains parent dir ('..'), format: "/./x/../y"
	{
		boost::filesystem::path path = boost::filesystem::path(".") / "x" / ".." / "y";
		string path_str = convert_path( path.string() );
		EXPECT_EQ("./x/../y", path_str);
	}

	// relative path starts without '.'
	{
		boost::filesystem::path path = boost::filesystem::path("x") / "y";
		string path_str = convert_path( path.string() );
		EXPECT_EQ("x/y", path_str);
	}
}



/**
 * Unit test of relative output_dir
 */
TEST(FilePath, output_relative) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    // Following was converted to Warning, that can not be simply tested.

    //EXPECT_DEATH( {FilePath("input/${INPUT}/init.in", FilePath::input_file);},
    //        "Creating FileName object before set_io_dirs is called.");

    boost::filesystem::create_directories("./main_root");
    FilePath::set_io_dirs(".", "./main_root", "variant_input", "../output_root");

    { // relative output_dir x relative output substitution
		FilePath fp = FilePath("output.vtk", FilePath::output_file);
		string str_fp = fp;

		// relative output substitution; conversion to string
		EXPECT_TRUE( (boost::filesystem::current_path() / "output_root/output.vtk") == boost::filesystem::path(str_fp) );

		// conversion to string
		EXPECT_TRUE( (boost::filesystem::current_path() / "output_root/output.vtk") == boost::filesystem::path(string(fp)) );
    }

    // relative output_dir x absolute output substitution, this case is not allowed
    EXPECT_THROW_WHAT( { FilePath(FilePath::get_absolute_working_dir()+"main_root/output.vtk", FilePath::output_file); },
    		FilePath::ExcAbsOutputPath, "Can not set absolute path" );

}


/**
 * Unit test of absolute output_dir
 */
TEST(FilePath, output_absolute) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    string abs_path = FilePath::get_absolute_working_dir();

    FilePath::set_io_dirs(".", "/main_root", "variant_input", abs_path+"output_root");

    {
		// relative output substitution; conversion to string
		string str = FilePath("output.vtk", FilePath::output_file);
		EXPECT_TRUE( (boost::filesystem::current_path() / "output_root/output.vtk") == boost::filesystem::path(str) );
    }

    {
		// absolute output substitution; conversion to string
		string str = FilePath(abs_path+"output_root/output.vtk", FilePath::output_file);
		EXPECT_TRUE( (boost::filesystem::current_path() / "output_root/output.vtk") == boost::filesystem::path(str) );
    }
}


/**
 * Unit test of relative input_dir
 */
TEST(FilePath, input_relative) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    FilePath::set_io_dirs(".","./main_root", "variant_input", "../output_root");

    { // relative input without substitution; conversion to string
		string str = FilePath("subdir/init.in", FilePath::input_file);
		EXPECT_TRUE( (boost::filesystem::current_path() / "main_root/subdir/init.in") == boost::filesystem::path(str) );
    }

    { // relative input with substitution; conversion to string
		string str = FilePath("subdir/${INPUT}/init.in", FilePath::input_file);
		EXPECT_TRUE( (boost::filesystem::current_path() / "main_root/subdir/variant_input/init.in") == boost::filesystem::path(str) );
    }

    { // absolute input without substitution; conversion to string
		string str = FilePath(FilePath::get_absolute_working_dir()+"subdir/init.in", FilePath::input_file);
		EXPECT_TRUE( (boost::filesystem::current_path() / "subdir/init.in") == boost::filesystem::path(str) );
    }

    { // absolute input with substitution; conversion to string
		string str = FilePath(FilePath::get_absolute_working_dir()+"subdir/${INPUT}/init.in", FilePath::input_file);
		EXPECT_TRUE( (boost::filesystem::current_path() / "subdir/variant_input/init.in") == boost::filesystem::path(str) );
    }
}


/**
 * Unit test of absolute input_dir
 */
TEST(FilePath, input_absolute) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    string abs_path = FilePath::get_absolute_working_dir();
    FilePath::set_io_dirs(".", "/main_root", "variant_input", abs_path+"output_root");

    {
        // relative input without substitution; conversion to string
        string str = FilePath("./subdir/init.in", FilePath::input_file);
		EXPECT_TRUE( (boost::filesystem::current_path() / "/main_root/./subdir/init.in") == boost::filesystem::path(str) );
    }

    {
        // relative input with substitution; conversion to string
        string str = FilePath("./subdir/${INPUT}/init.in", FilePath::input_file);
		EXPECT_TRUE( (boost::filesystem::current_path() / "/main_root/./subdir/variant_input/init.in") == boost::filesystem::path(str) );
    }

    {
        // absolute input without substitution; conversion to string
        string str = FilePath(abs_path+"subdir/init.in", FilePath::input_file);
		EXPECT_TRUE( (boost::filesystem::current_path() / "subdir/init.in") == boost::filesystem::path(str) );
    }

    {
        // absolute input with substitution; conversion to string
        string str = FilePath(abs_path+"subdir/${INPUT}/init.in", FilePath::input_file);
		EXPECT_TRUE( (boost::filesystem::current_path() / "subdir/variant_input/init.in") == boost::filesystem::path(str) );
    }

}


TEST(FilePath, create_output_dir) {
    //::testing::FLAGS_gtest_death_test_style = "threadsafe";
    boost::filesystem::remove_all("./my_output");
    FilePath::set_io_dirs(".","", "my_input", "my_output");
    EXPECT_TRUE(boost::filesystem::is_directory("./my_output"));
    FilePath fp("subdir/some_file.xyz", FilePath::output_file);
    fp.create_output_dir();
    EXPECT_TRUE(boost::filesystem::is_directory("./my_output/subdir"));
    boost::filesystem::remove_all("./my_output");
}
