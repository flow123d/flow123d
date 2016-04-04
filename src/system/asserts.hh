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
 * @file    asserts.hh
 * @brief   Definitions of ASSERTS.
 */

#ifndef ASSERTS_HH
#define ASSERTS_HH

#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include "system/exceptions.hh"

namespace feal {

/**
 * @brief Class defining debugging messages.
 *
 * Allows define assert, warning etc.
 */
class Assert : public ExceptionBase {
public:
	/// Constructor.
	Assert(const std::string& expression)
	: _FEAL_ASSERT_A (*this),
	  _FEAL_ASSERT_B (*this),
	  expression_(expression),
	  thrown_(false),
	  what_type_msg_("Program Error: Violated assert") {}

	/// Copy constructor.
	Assert(const Assert& other)
	: _FEAL_ASSERT_A (*this),
	  _FEAL_ASSERT_B (*this),
	  expression_(other.expression_),
	  file_name_(other.file_name_),
	  function_(other.function_),
	  line_(other.line_),
	  current_val_(other.current_val_),
	  thrown_(other.thrown_),
	  what_type_msg_(other.what_type_msg_) {}

	/// Destructor.
	~Assert() {
		if (!thrown_) this->error();
	}

	Assert& _FEAL_ASSERT_A; ///< clever macro A
	Assert& _FEAL_ASSERT_B; ///< clever macro B

	/// Adds name and value of variable
	template <typename T>
	Assert& add_value(T var_current_val, const char* var_name) {
		std::stringstream ss;
		ss << var_name << " : '" << var_current_val << "'";
		current_val_.push_back(ss.str());

		return *this;
	}

	/// Stores values for printing out line number, function, etc
	Assert& set_context(const char* file_name, const char* function, const int line)
	{
		file_name_ = std::string(file_name);
		function_ = std::string(function);
		line_ = line;

		return *this;
	}

	/// Generate error
	void error()
	{
		thrown_ = true;
		THROW( *this );
	}

	/// Generate warning
	void warning()
	{
		thrown_ = true;
		what_type_msg_ = "Warning:";
		std::cerr << this->what();
	}

	/// Print formatted assert message.
	void print_info(std::ostringstream &out) const override
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

protected:
    /// Override @p ExceptionBase::what_type_msg()
    std::string what_type_msg() const override {
    	return what_type_msg_;
    }

    std::string expression_;                  ///< Assertion expression
	std::string file_name_;                   ///< Actual file.
	std::string function_;                    ///< Actual function.
	int line_;                                ///< Actual line.
	std::vector< std::string > current_val_;  ///< Formated strings of names and values of given variables.
	bool thrown_;                             ///< Flag marked if Assert thrown error, warning, ...
	std::string what_type_msg_;               ///< String representation of message type (Program error, Warning, ...)

};

} // namespace feal

// Must define the macros afterwards
/// Internal clever macro A
#define _FEAL_ASSERT_A(x) _FEAL_ASSERT_OP(x, B)
/// Internal clever macro B
#define _FEAL_ASSERT_B(x) _FEAL_ASSERT_OP(x, A)
/// Internal clever macro recursion
#define _FEAL_ASSERT_OP(x, next) \
    _FEAL_ASSERT_A.add_value((x), #x)._FEAL_ASSERT_ ## next


#ifdef FLOW123D_DEBUG_ASSERTS
/// High-level macro
#define FEAL_ASSERT( expr) \
if ( !(expr) ) \
  feal::Assert( #expr).set_context( __FILE__, __func__, __LINE__)._FEAL_ASSERT_A
#else
#define FEAL_ASSERT( expr)
#endif

/**
 * Sources:
 * http://www.drdobbs.com/cpp/enhancing-assertions/184403745
 * https://gist.github.com/hang-qi/5308285
 * https://beliefbox.googlecode.com/svn-history/r825/trunk/src/core/SmartAssert.h
 */

#endif // ASSERTS_HH
