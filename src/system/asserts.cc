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
#include "system/logger.hh"


namespace feal {

/*******************************************************************
 * implementation of Exc_assert
 */

Exc_assert::Exc_assert()
: line_(0),
  what_type_msg_("Program Error: Violated assert! ") {}


void Exc_assert::print_info(std::ostringstream &out) const
{
	out << std::endl << "> In file: " << file_name_ << "(" << line_ << "): Throw in function " << function_ << std::endl;
	out << "> Expression: \'" << expression_ << "\'" << "\n";
	if (current_val_.size()) {
		out << "> Values:" << std::endl;
		for (auto val : current_val_) {
			out << "  " << val << std::endl;
		}
	}
}


std::string Exc_assert::what_type_msg() const {
	return what_type_msg_;
}


std::ostringstream &Exc_assert::form_message(std::ostringstream &converter) const {

    converter << "--------------------------------------------------------" << std::endl;
    converter << this->what_type_msg();
    print_info(converter);

    print_stacktrace(converter);
    converter << std::endl << "--------------------------------------------------------" << std::endl;

    return converter;
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


void Assert::error(std::string error_msg)
{
	exception_.what_type_msg_ += error_msg;
	thrown_ = true;
	THROW( exception_ );
}


void Assert::warning(std::string warning_msg)
{
	thrown_ = true;
	std::ostringstream info_str, stack_str;
	exception_.print_info(info_str);
	exception_.print_stacktrace(stack_str);
	WarningOut() << warning_msg << info_str.str() << StreamMask::log << stack_str.str();
}

} // namespace feal
