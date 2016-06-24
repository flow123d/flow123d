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

using namespace std;



TEST(BoostFileSystem, base_methods) {
	// current_path
	std::cout << "Current path: " << boost::filesystem::current_path() << std::endl;

	// simple absolute path
	{
		stringstream dir;
		dir << boost::filesystem::current_path().string() << DIR_DELIMITER << "x";
		if ( !boost::filesystem::is_directory(dir.str()) ) {
			boost::filesystem::create_directories(dir.str());
			EXPECT_TRUE( boost::filesystem::is_directory(dir.str()) );
			boost::filesystem::remove(dir.str());
			EXPECT_FALSE( boost::filesystem::is_directory(dir.str()) );
		}
	}

	// absolute path contains actual dir ('.'), format is current_path() + "/./x"
	{
		stringstream dir;
		dir << boost::filesystem::current_path().string() << DIR_DELIMITER << "x";
		stringstream dir_with_actual;
		dir_with_actual << boost::filesystem::current_path().string() << DIR_DELIMITER << "." << DIR_DELIMITER << "x";
		if ( !boost::filesystem::is_directory(dir_with_actual.str()) ) {
			boost::filesystem::create_directories(dir_with_actual.str());
			EXPECT_TRUE( boost::filesystem::is_directory(dir.str()) );
			boost::filesystem::remove(dir_with_actual.str());
			EXPECT_FALSE( boost::filesystem::is_directory(dir.str()) );
		}
	}

	// complicated absolute path in format current_path() + "/x/../y"
	{
		stringstream dir_x;
		dir_x << boost::filesystem::current_path().string() << DIR_DELIMITER << "x";
		stringstream dir_y;
		dir_y << boost::filesystem::current_path().string() << DIR_DELIMITER << "y";
		stringstream dir_y_full;
		dir_y_full << boost::filesystem::current_path().string()
				<< DIR_DELIMITER << "x" << DIR_DELIMITER << ".." << DIR_DELIMITER << "y";

		EXPECT_EQ( boost::filesystem::is_directory(dir_y.str()), boost::filesystem::is_directory(dir_y_full.str()));
		if ( !boost::filesystem::is_directory(dir_y_full.str()) ) {
			bool x_dir_exists = boost::filesystem::is_directory(dir_x.str());
			boost::filesystem::create_directories(dir_y_full.str());
			EXPECT_TRUE( boost::filesystem::is_directory(dir_y.str()) );
			EXPECT_TRUE( boost::filesystem::is_directory(dir_y_full.str()) );
			boost::filesystem::remove(dir_y_full.str());
			if (!x_dir_exists) {
				boost::filesystem::remove(dir_x.str());
			}
			EXPECT_FALSE( boost::filesystem::is_directory(dir_y.str()) );
			EXPECT_FALSE( boost::filesystem::is_directory(dir_y_full.str()) );
		}
	}

	// complicated absolute path in format current_path() + "/x/y"
	{
		stringstream dir;
		dir << boost::filesystem::current_path().string() << DIR_DELIMITER << "x" << DIR_DELIMITER << "y";
		stringstream dir_x;
		dir_x << boost::filesystem::current_path().string() << DIR_DELIMITER << "x";

		if ( !boost::filesystem::is_directory(dir.str()) ) {
			bool x_dir_exists = boost::filesystem::is_directory(dir_x.str());
			boost::filesystem::create_directories(dir.str());
			EXPECT_TRUE( boost::filesystem::is_directory(dir.str()) );
			boost::filesystem::remove(dir.str());
			if (!x_dir_exists) {
				boost::filesystem::remove(dir_x.str());
			}
			EXPECT_FALSE( boost::filesystem::is_directory(dir.str()) );
		}
	}

	// relative path starts with '.'
	{
		stringstream dir;
		dir << "." << DIR_DELIMITER << "x";
		stringstream abs_dir;
		abs_dir << boost::filesystem::current_path().string() << DIR_DELIMITER << "x";
		if ( !boost::filesystem::is_directory(dir.str()) ) {
			boost::filesystem::create_directories(dir.str());
			EXPECT_TRUE( boost::filesystem::is_directory(dir.str()) );
			EXPECT_TRUE( boost::filesystem::is_directory(abs_dir.str()) );
			boost::filesystem::remove(dir.str());
			EXPECT_FALSE( boost::filesystem::is_directory(dir.str()) );
			EXPECT_FALSE( boost::filesystem::is_directory(abs_dir.str()) );
		}
	}

	// relative path starts without '.'
	{
		stringstream dir;
		dir << "x";
		stringstream abs_dir;
		abs_dir << boost::filesystem::current_path().string() << DIR_DELIMITER << "x";
		if ( !boost::filesystem::is_directory(dir.str()) ) {
			boost::filesystem::create_directories(dir.str());
			EXPECT_TRUE( boost::filesystem::is_directory(dir.str()) );
			EXPECT_TRUE( boost::filesystem::is_directory(abs_dir.str()) );
			boost::filesystem::remove(dir.str());
			EXPECT_FALSE( boost::filesystem::is_directory(dir.str()) );
			EXPECT_FALSE( boost::filesystem::is_directory(abs_dir.str()) );
		}
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

    FilePath::set_io_dirs("./work_dir/xx", "./main_root", "variant_input", "../output_root");

    { // relative output_dir x relative output substitution
		FilePath fp = FilePath("output.vtk", FilePath::output_file);
		string str_fp = fp;

		// relative output substitution; conversion to string
		EXPECT_TRUE( (boost::filesystem::current_path() / "main_root/work_dir/output_root/output.vtk") == boost::filesystem::path(str_fp) );

		// conversion to string
		EXPECT_TRUE( (boost::filesystem::current_path() / "main_root/work_dir/output_root/output.vtk") == boost::filesystem::path(string(fp)) );
    }

    // relative output_dir x absolute output substitution, this case is not allowed
    EXPECT_THROW_WHAT( { FilePath(FilePath::get_absolute_working_dir()+"main_root/work_dir/xx/output.vtk", FilePath::output_file); },
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

    FilePath::set_io_dirs("./work_dir/xx","./main_root", "variant_input", "../output_root");

    { // relative input without substitution; conversion to string
		string str = FilePath("subdir/init.in", FilePath::input_file);
		EXPECT_TRUE( boost::filesystem::path("./main_root/subdir/init.in") == boost::filesystem::path(str) );
    }

    { // relative input with substitution; conversion to string
		string str = FilePath("subdir/${INPUT}/init.in", FilePath::input_file);
		EXPECT_TRUE( boost::filesystem::path("./main_root/subdir/variant_input/init.in") == boost::filesystem::path(str) );
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
    FilePath::set_io_dirs("/work_dir/xx", "/main_root", "variant_input", abs_path+"output_root");

    {
        // relative input without substitution; conversion to string
        string str = FilePath("./subdir/init.in", FilePath::input_file);
		EXPECT_TRUE( boost::filesystem::path("/main_root/./subdir/init.in") == boost::filesystem::path(str) );
    }

    {
        // relative input with substitution; conversion to string
        string str = FilePath("./subdir/${INPUT}/init.in", FilePath::input_file);
		EXPECT_TRUE( boost::filesystem::path("/main_root/./subdir/variant_input/init.in") == boost::filesystem::path(str) );
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
    boost::filesystem::remove_all("./work_dir");
    FilePath::set_io_dirs("./work_dir","", "my_input", "my_output");
    EXPECT_TRUE(boost::filesystem::is_directory("./work_dir/my_output"));
    FilePath fp("subdir/some_file.xyz", FilePath::output_file);
    fp.create_output_dir();
    EXPECT_TRUE(boost::filesystem::is_directory("./work_dir/my_output/subdir"));
    boost::filesystem::remove_all("./work_dir");
}
