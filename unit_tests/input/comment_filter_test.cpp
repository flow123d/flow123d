/*
 * comment_filter_test.cpp
 *
 *  Created on: May 7, 2012
 *      Author: jb
 */

#include <flow_gtest.hh>
#include <fstream>
#include <cerrno>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include "input/comment_filter.hh"

using namespace std;



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

#include "input/json_spirit/json_spirit.h"

TEST(Storage, comment_filter) {

    // test one line comments
    EXPECT_EQ("abc/def",filter("abc/def")); // no comment
    EXPECT_EQ("abc\ndef", filter("abc// xyz\ndef")); // linux
    EXPECT_EQ("abc\n\rdef", filter("abc// xyz\n\rdef")); // windows
    EXPECT_EQ("abc\ndef", filter("abc// \\\ndef")); // escaping is off in comment
    EXPECT_EQ("abc\"//\"\ndef",filter("abc\"//\"\ndef")); // quotation
    EXPECT_EQ("abc\"\\\"//\"\ndef",filter("abc\"\\\"//\"\ndef")); // escaping in quotation
    EXPECT_EQ("abc\\//\ndef",filter("abc\\//\ndef")); // escaping

    // test multiline comments
    EXPECT_EQ("abc\n\ndef", filter("abc/* \n\n*/def")); // preserve line numbers (linux)
    EXPECT_EQ("abc\n\r\n\rdef", filter("abc/* \n\r\n\r*/def")); // preserve line numbers (linux)
    EXPECT_EQ("abcdef", filter("abc/* \\*/def")); // escaping is off in comment
    EXPECT_EQ("abcdef",filter("abc/* \" * / */def")); // garbage in comment

    // check comments ending with Windows line ends
    EXPECT_EQ("\n\r\n\r", filter("// comment \n\r\n\r"));

    // Visual check of filter output
    std::ifstream in(string( string(UNIT_TESTS_SRC_DIR) + "/input/comment_filter_test.con").c_str(), std::ios::in | std::ios::binary);

    std::string contents;
    in.seekg(0, std::ios::end);
    contents.resize(in.tellg());
    in.seekg(0, std::ios::beg);
    in.read(&contents[0], contents.size());
    cout << filter(contents) << endl;

    // Test parsing filter output as JSON file
    in.seekg(0, std::ios::beg);

    boost::iostreams::filtering_istream filter_in;
    filter_in.push(Input::uncommenting_filter());
    filter_in.push(in);

    json_spirit::mValue main_node;
    try {
        json_spirit::read_or_throw( filter_in, main_node);
    } catch (json_spirit::Error_position &e ) {
        cout << "Error. line: " << e.line_ << " col: " << e.column_ << " reason: " << e.reason_ << endl;
        throw e;
    }

}


