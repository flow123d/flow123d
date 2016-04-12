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
 * @file    asserts.cc
 * @brief   Definitions of ASSERTS.
 */


#include "system/asserts.hh"


namespace feal {

/*******************************************************************
 * implementation of AssertException
 */

AssertException::AssertException()
: what_type_msg_("Program Error: Violated assert") {}


void AssertException::print_info(std::ostringstream &out) const
{
	out << std::endl << "> In file: " << file_name_ << "(" << line_ << "): Throw in function " << function_ << std::endl;
	out << "> Expression: \'" << expression_ << "\'" << std::endl;
	if (current_val_.size()) {
		out << "> Values:" << std::endl;
		for (auto val : current_val_) {
			out << "  " << val << std::endl;
		}
	}
}


std::string AssertException::what_type_msg() const {
	return what_type_msg_;
}


/*******************************************************************
 * implementation of Assert
 */

Assert::~Assert() {
	if (!thrown_) {
		// We can't throw exception in destructor, we need use this construction
		std::cerr << exception_.what();
		abort();
	}
}


Assert& Assert::set_context(const char* file_name, const char* function, const int line)
{
	exception_.file_name_ = std::string(file_name);
	exception_.function_ = std::string(function);
	exception_.line_ = line;

	return *this;
}


void Assert::error()
{
	thrown_ = true;
	THROW( exception_ );
}


void Assert::warning()
{
	thrown_ = true;
	exception_.what_type_msg_ = "Warning:";
	std::cerr << exception_.what();
}

} // namespace feal
