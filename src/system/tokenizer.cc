/*
 * tokenizer.cc
 *
 *  Created on: Nov 9, 2012
 *      Author: jb
 */

#include "system/tokenizer.hh"
#include "boost/algorithm/string/trim.hpp"

Tokenizer::Tokenizer( istream &in)
: in_(in), line_counter_(0), line_tokenizer_(line_,  boost::char_separator<char>("\t \n"))
{}



void Tokenizer::next_line() {
    line_="";
    while ( line_ == "") { std::getline( in_, line_); boost::trim( line_ ); line_counter_++; }
    line_tokenizer_.assign(line_);
    tok_ = line_tokenizer_.begin();
    position = 0;
}

