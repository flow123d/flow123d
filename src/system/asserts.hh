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
 * @brief Helper class.
 *
 * Stores data of assert and allow throws exception.
 */
class AssertException : public ExceptionBase {
	friend class Assert;
public:
	// Default constructor
	AssertException()
	: what_type_msg_("Program Error: Violated assert") {}

	/// Destructor.
	~AssertException() {}

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
	std::string what_type_msg_;               ///< String representation of message type (Program error, Warning, ...)
};

/**
 * @brief Class defining debugging messages.
 *
 * Allows define assert, warning etc. either only for debug mode or for release mode also.
 *
 * Definition of asserts is designed using macros FEAL_ASSERT and DEBUG_ASSERT. First macro
 * is used for both modes, second is only for debug. Definition allows to printout given
 * variables too.
 *
 * Examples of usage:
 *
 * 1) We expect empty stings 's1' and 's2', if condition is not satisfied exception will be
 *    thrown. Condition is indicated in first parentheses, variables designed for printout
 *    follow as (s1)(s2). Each variable must be defined in separate parentheses. The last
 *    step is calling of the appropriate assert type, in this case error(). This assert is
 *    performed for debug and release mode.
 @code
    std::string s1, s2;
    ...
    FEAL_ASSERT(s1.empty() && s2.empty())(s1)(s2).error();
 @endcode
 *
 * 2) This example is same as previous, but assert is performed only for debug mode.
 @code
    DEBUG_ASSERT(s1.empty() && s2.empty())(s1)(s2).error();
 @endcode
 *
 * 3) Example is same as case 1). Assert type error is called automatically if any other is
 *    not listed.
 @code
    FEAL_ASSERT(s1.empty() && s2.empty())(s1)(s2);
 @endcode
 *
 * 4) Example with same condition as all previous but with other type - warning. Any exception
 *    is not thrown, only warning is printed.
 @code
    FEAL_ASSERT(s1.empty() && s2.empty())(s1)(s2).warning();
 @endcode
 */
class Assert {
public:
	/// Constructor.
	Assert(const std::string& expression)
	: _FEAL_ASSERT_A (*this),
	  _FEAL_ASSERT_B (*this),
	  thrown_(false)
	{
		exception_.expression_ = expression;
	}

	/// Copy constructor.
	Assert(const Assert& other)
	: _FEAL_ASSERT_A (*this),
	  _FEAL_ASSERT_B (*this),
	  exception_(other.exception_),
	  thrown_(other.thrown_) {}

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
		exception_.current_val_.push_back(ss.str());

		return *this;
	}

	/// Stores values for printing out line number, function, etc
	Assert& set_context(const char* file_name, const char* function, const int line)
	{
		exception_.file_name_ = std::string(file_name);
		exception_.function_ = std::string(function);
		exception_.line_ = line;

		return *this;
	}

	/// Generate error
	void error()
	{
		thrown_ = true;
		THROW( exception_ );
	}

	/// Generate warning
	void warning()
	{
		thrown_ = true;
		exception_.what_type_msg_ = "Warning:";
		std::cerr << exception_.what();
	}

protected:
    AssertException exception_;               ///< Exception object
	bool thrown_;                             ///< Flag marked if Assert thrown error, warning, ...

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


/// Definition of assert for debug and release mode
#define FEAL_ASSERT( expr) \
if ( !(expr) ) \
  feal::Assert( #expr).set_context( __FILE__, __func__, __LINE__)._FEAL_ASSERT_A

/// Definition of assert for debug mode only
#ifdef FLOW123D_DEBUG_ASSERTS
#define DEBUG_ASSERT( expr) \
if ( !(expr) ) \
  feal::Assert( #expr).set_context( __FILE__, __func__, __LINE__)._FEAL_ASSERT_A
#else
#define DEBUG_ASSERT( expr)
#endif

/**
 * Sources:
 * http://www.drdobbs.com/cpp/enhancing-assertions/184403745
 * https://gist.github.com/hang-qi/5308285
 * https://beliefbox.googlecode.com/svn-history/r825/trunk/src/core/SmartAssert.h
 */

#endif // ASSERTS_HH
