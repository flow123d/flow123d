/*
 * field_speed_test.cpp
 *
 *  Created on: Feb 3, 2013
 *      Author: jb
 */

#define TEST_USE_MPI

#include <flow_gtest_mpi.hh>

#ifdef RUN_UNIT_BENCHMARKS

#include <fstream>

#include "system/tokenizer.hh"
#include "system/file_path.hh"
#include "system/sys_profiler.hh"


static const unsigned int loop_call_count = 100000;
static const unsigned int file_line_count = 100;
static const std::string file_line_text = "\"some_text_line\"";

struct TestData {
	Tokenizer::Position pos;
	char val;

	TestData(Tokenizer::Position p, char v)
	: pos(p), val(v) {}
};


TEST(TokenizerPosition, compare_speed) {
	::testing::FLAGS_gtest_death_test_style = "threadsafe";

	FilePath::set_io_dirs(".", UNIT_TESTS_SRC_DIR, "", ".");

	/* // create file
	ofstream fout( FilePath("tokenizer_speed.txt", FilePath::output_file) );
	//fout << "1" << " " << "text" << std::endl;
	for (unsigned int i=0; i<file_line_count; i++) {
		fout << "1" << " " << file_line_text << std::endl;
	}
	fout << flush;
	fout.close(); */

	FilePath tok_file("/fields/simplest_cube_data.msh", FilePath::input_file);

	std::vector<TestData> test_data;

	// prepare positions in file and read test data
	{
	    Tokenizer tok(tok_file);

	    while (!tok.eof()) {
	    	Tokenizer::Position pos = tok.get_position();
	    	tok.next_line(false);
	    	test_data.push_back(TestData(pos, (*tok)[0]));
	    }
	    test_data.pop_back();
	}

	Profiler::initialize();

	{
	    Tokenizer tok(tok_file);
	    unsigned int index;

	    START_TIMER("tokenizer");
	    for (unsigned int i=0; i<loop_call_count; i++) {
			index = (i*101) % test_data.size();
			tok.set_position( test_data[index].pos );
			tok.next_line(false);
			EXPECT_EQ(test_data[index].val, (*tok)[0]);
	    }
	    END_TIMER("tokenizer");
	}

	{
		std::ifstream * binary_file = new ifstream;
		binary_file->open( string(tok_file).c_str(), ios::in | ios::binary );
		unsigned int index;
		char * c;

		START_TIMER("binary_file");
		for (unsigned int i=0; i<loop_call_count; i++) {
			index = (i*101) % test_data.size();
			binary_file->seekg( test_data[index].pos.file_position() );
			binary_file->read(c, 1);
			if (c[0] >= '0' && c[0] <= '9') EXPECT_EQ(test_data[index].val, c[0]);
		}
		END_TIMER("binary_file");

		binary_file->close();
	}

	Profiler::instance()->output(MPI_COMM_WORLD, cout);
	Profiler::uninitialize();
}

#endif // RUN_UNIT_BENCHMARKS
