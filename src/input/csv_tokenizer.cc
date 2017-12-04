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
 * @file    csv_tokenizer.cc
 * @brief
 */

#include "input/csv_tokenizer.hh"
#include "input/reader_internal_base.hh"
#include <boost/lexical_cast.hpp>

using namespace boost;
using namespace Input;


CSVTokenizer::CSVTokenizer(const FilePath &fp, std::string field_separator)
: Tokenizer(fp, Tokenizer::Separator("\\", field_separator, "\""))
{}



CSVTokenizer::CSVTokenizer( std::istream &in, std::string field_separator)
: Tokenizer(in, Tokenizer::Separator("\\", field_separator, "\""))
{}


unsigned int CSVTokenizer::get_n_lines()
{
    unsigned int n_lines = 0; // number of lines
    this->set_position( Tokenizer::Position() );
    while ( !this->eof() ) {
    	this->next_line(false);
    	n_lines++;
    }
    if (this->line().size()==0) n_lines--; // removes last line if it is empty
	return n_lines;
}


void CSVTokenizer::skip_header(unsigned int n_head_lines)
{
    this->set_position( Tokenizer::Position() );
    for (unsigned int i=0; i<n_head_lines; i++) {
    	this->next_line(false);
    }
}


int CSVTokenizer::get_int_val()
{
	try {
		return lexical_cast<int>( *(*this) );
	} catch (bad_lexical_cast &) {
		THROW( ReaderInternalBase::ExcWrongCsvFormat() << ReaderInternalBase::EI_TokenizerMsg(this->position_msg()) );
	}
}


double CSVTokenizer::get_double_val()
{
	try {
		return lexical_cast<double>( *(*this) );
	} catch (bad_lexical_cast &) {
		THROW( ReaderInternalBase::ExcWrongCsvFormat() << ReaderInternalBase::EI_Specification("Wrong double value")
				<< ReaderInternalBase::EI_TokenizerMsg(this->position_msg()) );
	}
}


std::string CSVTokenizer::get_string_val()
{
	try {
		return *(*this);
	} catch (bad_lexical_cast &) {
		THROW( ReaderInternalBase::ExcWrongCsvFormat() << ReaderInternalBase::EI_TokenizerMsg(this->position_msg()) );
	}
}
