/*
 * field_speed_test.cpp
 *
 *  Created on: Feb 3, 2013
 *      Author: jb
 */

#include <flow_gtest.hh>

#include "system/global_defs.h"


#ifdef FLOW123D_RUN_UNIT_BENCHMARKS

#include <fstream>

#include "system/tokenizer.hh"
#include "system/file_path.hh"
#include "system/sys_profiler.hh"

#include "boost/lexical_cast.hpp"


static const unsigned int loop_call_count =  500000;
static const unsigned int file_line_count = 1000000;
static const unsigned int line_step_count =  345001;
static const std::string file_line_text = "\"some_text_line\"";


TEST(TokenizerPosition, compare_speed) {
	::testing::FLAGS_gtest_death_test_style = "threadsafe";

	FilePath::set_io_dirs(".", UNIT_TESTS_SRC_DIR, "", ".");
	std::vector<Tokenizer::Position> position_data;

	// create file, fill vector of positions
	{

		ofstream fout( FilePath("tokenizer_speed.txt", FilePath::output_file) );
		for (unsigned int i=0; i<file_line_count; i++) {
			position_data.push_back( Tokenizer::Position(fout.tellp(), (i+1), 0) );
			fout << (i+1) << " " << file_line_text << std::endl;
		}
		fout << flush;
		fout.close();
	}
	EXPECT_EQ(position_data.size(), file_line_count);

	Profiler::initialize();
	FilePath in_file("./system/tokenizer_speed.txt", FilePath::input_file);

	// read data by tokenizer
	{
		Tokenizer tok(in_file);
		unsigned int index;
		unsigned int val;

	    START_TIMER("tokenizer");
	    for (unsigned int i=0; i<loop_call_count; i++) {
			index = (i*line_step_count) % file_line_count;
			tok.set_position( position_data[index] );
			tok.next_line(false);
			val = boost::lexical_cast<unsigned int>(*tok); ++tok;
			EXPECT_EQ(position_data[index].line_counter_, val);
	    }
	    END_TIMER("tokenizer");

	    // test of reading after reaching the EOF
	    tok.set_position( position_data[file_line_count-1] );
	    while (tok.next_line(false)) {}
	    tok.set_position( position_data[0] );
	    tok.next_line(false);
	    val = boost::lexical_cast<unsigned int>(*tok); ++tok;
	    EXPECT_EQ(position_data[0].line_counter_, val);
	}

	// read data same as from binary file
	{
		std::ifstream binary_file( string(in_file).c_str() );
		unsigned int index;
		unsigned int val;

		START_TIMER("binary_file");
		for (unsigned int i=0; i<loop_call_count; i++) {
			index = (i*line_step_count) % file_line_count;
			binary_file.seekg( position_data[index].file_position_ );
			binary_file >> val;
			EXPECT_EQ(position_data[index].line_counter_, val);
		}
		END_TIMER("binary_file");

		binary_file.close();
	}

	Profiler::instance()->output(cout);
	Profiler::uninitialize();
}

#endif // FLOW123D_RUN_UNIT_BENCHMARKS
