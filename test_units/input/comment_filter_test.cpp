/*
 * comment_filter_test.cpp
 *
 *  Created on: May 7, 2012
 *      Author: jb
 */

#include <gtest/gtest.h>
#include <fstream>
#include <cerrno>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include "input/comment_filter.hh"

using namespace std;


std::string get_file_contents(const string &filename)
{
  std::ifstream in(filename, std::ios::in | std::ios::binary);
  if (in)
  {
    std::string contents;
    in.seekg(0, std::ios::end);
    contents.resize(in.tellg());
    in.seekg(0, std::ios::beg);
    in.read(&contents[0], contents.size());
    in.close();
    return(contents);
  }
  throw(errno);
}


string filter(const string &in) {
    namespace io = boost::iostreams;

    istringstream in_stream(in);
    ostringstream out_stream;
    io::filtering_istream fin;

    fin.push(Input::uncommenting_filter());
    fin.push(in_stream);

    char c;
    while (fin.get(c)) out_stream.put(c);

    return out_stream.str();
}

TEST(Storage, comment_filter) {


    // This seems to produce correct output, but to be sure we have small inputs compared to correct outputs
    cout << filter(get_file_contents(string(UNIT_TESTS_SRC_DIR) + "/input/comment_filter_test.con"));

    // check comments ending with Windows line ends
    EXPECT_EQ("\n\r\n\r", filter("# comment \n\r\n\r"));

}


