/*
 * tokenizer_test.cpp
 *
 *  Created on: Nov 13, 2012
 *      Author: jb
 */




#include <flow_gtest.hh>
#include <sstream>
#include <string>
#include "system/tokenizer.hh"
#include "system/file_path.hh"

using  namespace std;

string input = R"CODE(0 1
"a b c" ,  a b c

$Elements_
 0  1 2 3
$EndElements 
something
)CODE";


void test_tokenizer(Tokenizer &tok) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    EXPECT_TRUE( tok.eol() );
    EXPECT_EQ(0, tok.line_num());

    // test on the first line
    tok.next_line();
    EXPECT_EQ(1, tok.line_num());
    EXPECT_EQ("0", *tok);
    EXPECT_FALSE( tok.eol() );
    ++tok;
    EXPECT_EQ("1", *tok);
    EXPECT_FALSE( tok.eol() );
    ++tok;
    EXPECT_TRUE( tok.eol() );

    EXPECT_DEATH( { *tok; }, "Missing token, .* position: '2'");

    // tokenizer should not lead to seg fault if we iterate over the end of line
    ++tok;
    EXPECT_TRUE( tok.eol() );

    // second line
    tok.next_line();
    EXPECT_EQ(2, tok.line_num());
    EXPECT_EQ("a b c", *tok);
    ++tok;
    EXPECT_EQ(",", *tok);

    // test skip_to
    tok.skip_to("$Elements");
    EXPECT_EQ(4, tok.line_num());

    // should remain on the same line
    tok.skip_to("$Elements");
    EXPECT_EQ(4, tok.line_num());
    EXPECT_EQ("$Elements_", *tok ); // same line readable from beginning
    tok.next_line();
    tok.next_line();

    EXPECT_TRUE( tok.next_line() );
    EXPECT_EQ("something", *tok);


}


TEST(Tokenizer, from_stream) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    std::stringstream ss(input);
    Tokenizer tok(ss);
    test_tokenizer(tok);
    EXPECT_FALSE( tok.eof() );
    EXPECT_FALSE( tok.next_line() ); // no next line
    EXPECT_TRUE( tok.eof() );
}



TEST(Tokenizer, from_file) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    // test Tokenizer constructed from FilePAth
    FilePath tok_file( string(UNIT_TESTS_SRC_DIR) + "/system/tokenizer_test_input", FilePath::input_file);
    Tokenizer tok(tok_file);
    test_tokenizer(tok);
    // in the file we have removed last empty line to test correct behavior in this case
    EXPECT_EQ(7, tok.line_num());
    EXPECT_TRUE( tok.eof() );
    EXPECT_FALSE( tok.next_line() );
    EXPECT_EQ(7, tok.line_num());

}
