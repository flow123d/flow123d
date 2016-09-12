/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    tokenizer.cc
 * @brief   
 */

#include <string>
#include <fstream>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/trim.hpp>

#include "system/global_defs.h"
#include "system/system.hh"
#include "system/tokenizer.hh"
#include "system/file_path.hh"

using namespace std;

Tokenizer::Tokenizer(const FilePath &fp)
: f_name_(fp),
  own_stream_(nullptr),
  in_(nullptr),
  comment_pattern_(""),
  position_(0, 0, 0),
  separator_("\\"," \t","\""),
  line_tokenizer_(line_,  separator_)
{
    own_stream_ = new ifstream;
    fp.open_stream(*own_stream_);
    in_ = own_stream_;
}



Tokenizer::Tokenizer( std::istream &in)
: f_name_("__anonymous_stream__"),
  own_stream_(nullptr),
  in_( &in ),
  comment_pattern_(""),
  position_(0, 0, 0),
  separator_("\\"," \t","\""),
  line_tokenizer_(line_,  separator_)
{}


void Tokenizer::set_comment_pattern( const std::string &pattern) {
    comment_pattern_=pattern;
}


bool Tokenizer::skip_to(const std::string& pattern, const std::string &end_search_pattern)
{
	//OLD_ASSERT( in_->good(), "Tokenizer stream (for file: %s) is not ready for i/o operations. Perhaps missing check about correct open.\n", f_name_.c_str());
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
    	WarningOut() << "Remaining tokens, file '" << f_name_ << "', line '" << line_num()
    			<< "', after token #" << position_.line_position_ << "." << std::endl;
    }

    if (eof()) return false; // we are sure that at least one getline will occur

    line_="";
    // skip empty lines
    while ( ! eof() && line_ == "") {
        std::getline( *in_, line_);
        position_.line_counter_++;
        // check failure bits
        if (in_->bad()) xprintf(Err, "Can not read from stream, file: '%s', line: '%d'\n", f_name_.c_str(), line_num());
        boost::trim( line_ );
        // if pattern is set and beginning of line match it
        if (comment_pattern_.size() && 0==line_.compare(0, comment_pattern_.size(), comment_pattern_) ) line_="";
    }
    if (! in_->fail() ) { // allow only eof state after any getline
        set_tokenizer();
        return true;
    }

    return false;
}



const std::string & Tokenizer::operator *() const
{
    if ( eol() ) xprintf(UsrErr, "Missing token, file: '%s', line: '%d', position: '%d'.\n", f_name_.c_str(), line_num(), position_.line_position_);
    return *tok_;
}



void Tokenizer::set_tokenizer()
{
        line_tokenizer_.assign(line_);
        tok_ = line_tokenizer_.begin();
        position_.line_position_ = 0;
        // skip leading separators
        while (! eol() && (*tok_).size()==0 ) {position_.line_position_++; ++tok_;}

}



string Tokenizer::position_msg() const {
    stringstream ss;
    ss << "token: " << pos() << ", line: " << line_num() << ", in file '" << f_name() << "'";
    return ss.str();
}


const Tokenizer::Position Tokenizer::get_position()
{
	position_.file_position_ = in_->tellg();
	return position_;
}


void Tokenizer::set_position(const Tokenizer::Position pos)
{
	in_->clear();
	in_->seekg(pos.file_position_);
	position_ = pos;
}


Tokenizer::~Tokenizer() {
    if (own_stream_ != NULL) delete own_stream_; // this also close the input file
}

