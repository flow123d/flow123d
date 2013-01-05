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

Tokenizer::Tokenizer(const FilePath &fp)
: f_name_(fp),
  own_stream_(NULL),
  line_counter_(0), position_(0),
  separator_("\\"," \t","\""),
  line_tokenizer_(line_,  separator_)
{
    in_ = own_stream_ = new ifstream;
    own_stream_->open( string(fp).c_str() );
    // check correct openning
    INPUT_CHECK(! own_stream_->fail(), "Can not open input file '%s'.\n", f_name_.c_str() );
    //own_stream_->exceptions ( ifstream::failbit | ifstream::badbit );

}



Tokenizer::Tokenizer( std::istream &in)
: f_name_("__anonymous_stream__"),
  own_stream_(NULL),
  in_( &in ),
  line_counter_(0), position_(0),
  separator_("\\"," \t","\""),
  line_tokenizer_(line_,  separator_)
{}



bool Tokenizer::skip_to(const std::string& pattern, const std::string &end_search_pattern)
{
    ASSERT( in_->good(), "Tokenizer stream (for file: %s) is not ready for i/o operations. Perhaps missing check about correct open.\n", f_name_.c_str());
    bool end_search= (end_search_pattern.size() > 0);

    while (! eof()) {
        if (line_.find(pattern)!=string::npos ) {
            set_tokenizer();
            return true;
        }
        if ( end_search && line_.find(end_search_pattern)!=string::npos ) return false;
        next_line(false);
    }
    return false;
}



bool Tokenizer::next_line(bool assert_for_remaining_tokens) {
    // input assert about remaining tokens
    if (assert_for_remaining_tokens && (! eol() )) {
        //DBGMSG("Line: '%s'\n", line_.c_str());
        xprintf(Warn, "Remaining tokens, file '%s', line '%d', after token #%d.\n", f_name_.c_str(), line_num(), position_);
    }

    if (eof()) return false; // we are sure that at least one getline will occur

    line_="";
    // skip empty lines
    while ( ! eof() && line_ == "") {
        std::getline( *in_, line_);
        // check failure bits
        if (in_->bad()) xprintf(Err, "Can not read from stream, file: '%s', line: '%d'\n", f_name_.c_str(), line_num());
        boost::trim( line_ );
        line_counter_++;
    }
    if (! in_->fail() ) { // allow only eof state after any getline
        set_tokenizer();
        return true;
    } else {
        DBGMSG("Line: '%s'\n", line_.c_str());
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
        // skip leading separators
        while (! eol() && (*tok_).size()==0 ) {position_++; ++tok_;}

}



string Tokenizer::position_msg() const {
    stringstream ss;
    ss << "token: " << pos() << ", line: " << line_num() << ", in file '" << f_name() << "'";
    return ss.str();
}


Tokenizer::~Tokenizer() {
    if (own_stream_ != NULL) delete own_stream_; // this also close the input file
}

