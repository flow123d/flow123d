/*
 * tokenizer.cc
 *
 *  Created on: Nov 9, 2012
 *      Author: jb
 */

#include <string>
#include <fstream>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/trim.hpp>

#include "system/system.hh"
#include "system/tokenizer.hh"
#include "system/file_path.hh"

using namespace std;

Tokenizer::Tokenizer( FilePath &fp)
: f_name_(fp),
  own_stream_(NULL),
  line_counter_(0), position_(0),
  line_tokenizer_(line_,  boost::char_separator<char>("\t \n"))
{
    in_ = own_stream_ = new ifstream( string(fp).c_str() ); // open a stream (?? error checking)

}



Tokenizer::Tokenizer( istream &in)
: f_name_("__anonymous_stream__"),
  own_stream_(NULL),
  in_( &in ),
  line_counter_(0), position_(0),
  line_tokenizer_(line_,  boost::char_separator<char>("\t \n"))
{}



bool Tokenizer::skip_to(const std::string& pattern)
{
    ASSERT( in_->good(), "Tokenizer stream (for file: %s) is not ready for i/o operations. Perhaps missing check about correct open.\n", f_name_.c_str());

    while (! eof() &&  line_.find(pattern)==string::npos ) next_line();
    if (! eof()) {
        set_tokenizer();
        return true;
    }
    return false;
}



bool Tokenizer::next_line() {
    // input assert about remaining tokens
    if ( ! eol() ) xprintf(Warn, "Remaining tokens, file '%s', line '%d', after token #%d.\n", f_name_.c_str(), line_num(), position_);

    line_="";
    while ( ! eof() && line_ == "") { std::getline( *in_, line_); boost::trim( line_ ); line_counter_++; }
    if (! eof() ) {
        set_tokenizer();
        return true;
    }
    return false;
}



const std::string & Tokenizer::operator *() const
{
    if ( eol() ) xprintf(UsrErr, "Missing token, file: '%s', line: '%d', position: '%d'.\n", f_name_.c_str(), line_num(), position_);
    return *tok_;
}



void Tokenizer::set_tokenizer()
{
        line_tokenizer_.assign(line_);
        tok_ = line_tokenizer_.begin();
        position_ = 0;   
}



Tokenizer::~Tokenizer() {
    if (own_stream_ != NULL) delete own_stream_;
}

