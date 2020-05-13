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
#include <vector>                // for vector
#include "system/exceptions.hh"
#include "system/fmt/posix.h"

namespace feal {

/**
 * @brief Helper class.
 *
 * Stores data of assert and allow throws exception.
 */
class Exc_assert : public ExceptionBase {
	friend class Assert;
public:
	// Default constructor
	Exc_assert();

	/// Destructor.
	~Exc_assert() {}

	/// Print formatted assert message.
	void print_info(std::ostringstream &out) const override;

protected:
    /// Override @p ExceptionBase::what_type_msg()
    std::string what_type_msg() const override;

    /// Override @p ExceptionBase::form_message()
    std::ostringstream &form_message(std::ostringstream &) const override;

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
 * Definition of asserts is designed using macros FEAL_ASSERT and FEAL_ASSERT_DBG. First macro
 * is used for both modes, second is only for debug. Definition allows to printout given
 * variables too. For simplifying are designed shorter names of macros ASSERT and ASSERT_DBG,
 * these names can be used if thea aren't in conflicts with external libraries.
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
    FEAL_ASSERT(s1.empty() && s2.empty())(s1)(s2).error("Both strings must be empty!");
 @endcode
 *
 * 2) Parameter of error method is optional.
 @code
    FEAL_ASSERT(s1.empty() && s2.empty())(s1)(s2).error();
 @endcode
 *
 * 3) This example is same as previous, but assert is performed only for debug mode.
 @code
    FEAL_ASSERT_DBG(s1.empty() && s2.empty())(s1)(s2).error("Both strings must be empty!");
 @endcode
 *
 * 4) Example is same as case 1). Assert type error is called automatically if any other is
 *    not listed. This case is not recommended, rather use explicitly calling of error() method.
 @code
    FEAL_ASSERT(s1.empty() && s2.empty())(s1)(s2);
 @endcode
 *
 * 5) Example with same condition as all previous but with other type - warning. Any exception
 *    is not thrown, only warning is printed. Parameter of warning method is optional.
 @code
    FEAL_ASSERT(s1.empty() && s2.empty())(s1)(s2).warning("Both strings should be empty!");
 @endcode
 *
 * For simplifying we have defined several macros for comparsion of values:
 *  - ASSERT_LT(a, b) or ASSERT_LT_DBG(a, b) ... check if (a < b)
 *  - ASSERT_LE(a, b) or ASSERT_LE_DBG(a, b) ... check if (a <= b)
 *  - ASSERT_GT(a, b) or ASSERT_GT_DBG(a, b) ... check if (a > b)
 *  - ASSERT_GE(a, b) or ASSERT_GE_DBG(a, b) ... check if (a >= b)
 *  - ASSERT_EQ(a, b) or ASSERT_EQ_DBG(a, b) ... check if (a == b)
 *  - ASSERT_PTR(obj) or ASSERT_PTR_DBG(obj) ... check if obj is non-null pointer
 *
 * All macros allow easier declaration of assert. Following example shows declarations of same
 * cases with usage of different macros:
 @code
    ASSERT_LT( idx, arr.size() ).error("Index out of array!");
    FEAL_ASSERT( idx < arr.size() )(idx)(arr.size()).error("Index out of array!");
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
		exception_.frames_to_cut_ = { "feal", "Assert"};
	}

	/// Copy constructor.
	Assert(const Assert& other)
	: _FEAL_ASSERT_A (*this),
	  _FEAL_ASSERT_B (*this),
	  exception_(other.exception_),
	  thrown_(other.thrown_) {}

	/// Destructor.
	~Assert();

	Assert& _FEAL_ASSERT_A; ///< clever macro A
	Assert& _FEAL_ASSERT_B; ///< clever macro B

	/// Adds name and value of variable
	template <typename T>
	Assert& add_value(T var_current_val, const char* var_name) {
		std::stringstream ss;
		ss << var_name << " : ";
		ss << "'" << var_current_val << "'"; // Can throw exception if type T hasn't overloading << operator
		exception_.current_val_.push_back(ss.str());

		return *this;
	}

	/// Stores values for printing out line number, function, etc
	Assert& set_context(const char* file_name, const char* function, const int line);

	/// Generate error with given message
	void error(std::string error_msg = "");

	/// Generate warning with given message
	void warning(std::string warning_msg = "");

protected:
    Exc_assert exception_;               ///< Exception object
	bool thrown_;                        ///< Flag marked if Assert thrown error, warning, ...

};



/**
 * Helper class defined empty code.
 *
 * Usage only in FEAL_ASSERT_DBG macro if FLOW123D_DEBUG_ASSERTS is off.
 */
class AssertNull {
public:
	/// Constructor.
	AssertNull()
	: _FEAL_ASSERT_A (*this),
	  _FEAL_ASSERT_B (*this) {}

	/// Copy constructor.
	AssertNull(FMT_UNUSED const Assert& other)
	: _FEAL_ASSERT_A (*this),
	  _FEAL_ASSERT_B (*this) {}

	/// Destructor.
	~AssertNull() {}

	AssertNull& _FEAL_ASSERT_A; ///< clever macro A
	AssertNull& _FEAL_ASSERT_B; ///< clever macro B

	/// Empty method, only guarantees consistent code
	template <typename T>
	inline AssertNull& add_value(FMT_UNUSED T var_current_val, FMT_UNUSED const char* var_name) {
		return *this;
	}

	/// Empty method, only guarantees consistent code
	inline void error(FMT_UNUSED std::string error_msg = "") {}

	/// Empty method, only guarantees consistent code
	inline void warning(FMT_UNUSED std::string warning_msg = "") {}
};

} // namespace feal



/****************************************************************
 * Definitions of macros
 */

/**
 * Definitions of internal macros
 */
/// Internal clever macro A
#define _FEAL_ASSERT_A(x) _FEAL_ASSERT_OP(x, B)
/// Internal clever macro B
#define _FEAL_ASSERT_B(x) _FEAL_ASSERT_OP(x, A)
/// Internal clever macro recursion
#define _FEAL_ASSERT_OP(x, next) \
    _FEAL_ASSERT_A.add_value((x), #x)._FEAL_ASSERT_ ## next


/**
 * Undefining redefinitions of comparsion macros.
 *
 * You can define FEAL_OVERRIDE_ASSERTS (e.g. in unit tests) and undefine macros with same name in external libraries.
 */
#ifdef ASSERT_LT
    #ifdef FEAL_OVERRIDE_ASSERTS
        #undef ASSERT_LT
    #else
        #warning "ASSERT_LT already defined out of FEAL."
    #endif
#endif

#ifdef ASSERT_LE
    #ifdef FEAL_OVERRIDE_ASSERTS
        #undef ASSERT_LE
    #else
        #warning "ASSERT_LE already defined out of FEAL."
    #endif
#endif

#ifdef ASSERT_GT
    #ifdef FEAL_OVERRIDE_ASSERTS
        #undef ASSERT_GT
    #else
        #warning "ASSERT_GT already defined out of FEAL."
    #endif
#endif

#ifdef ASSERT_GE
    #ifdef FEAL_OVERRIDE_ASSERTS
        #undef ASSERT_GE
    #else
        #warning "ASSERT_GE already defined out of FEAL."
    #endif
#endif

#ifdef ASSERT_EQ
    #ifdef FEAL_OVERRIDE_ASSERTS
        #undef ASSERT_EQ
    #else
        #warning "ASSERT_EQ already defined out of FEAL."
    #endif
#endif


/**
 * Definitions of assert macros
 */
/// Definition of assert for debug and release mode
#define FEAL_ASSERT( expr) \
if ( !(expr) ) \
  feal::Assert( #expr).set_context( __FILE__, __func__, __LINE__)._FEAL_ASSERT_A

/// Definition of assert for debug mode only
#ifdef FLOW123D_DEBUG_ASSERTS
    #define FEAL_ASSERT_DBG( expr) \
    if ( !(expr) ) \
      feal::Assert( #expr).set_context( __FILE__, __func__, __LINE__)._FEAL_ASSERT_A
#else
    #define FEAL_ASSERT_DBG( expr) \
    if ( !(expr) ) \
      feal::AssertNull()._FEAL_ASSERT_A
#endif

/// Definition of comparative assert macro (Less Than)
#define ASSERT_LT(a, b) \
	FEAL_ASSERT(a < b)(a)(b)

/// Definition of comparative assert macro (Less Than) only for debug mode
#define ASSERT_LT_DBG(a, b) \
	FEAL_ASSERT_DBG(a < b)(a)(b)

/// Definition of comparative assert macro (Less or Equal)
#define ASSERT_LE(a, b) \
	FEAL_ASSERT(a <= b)(a)(b)

/// Definition of comparative assert macro (Less or Equal) only for debug mode
#define ASSERT_LE_DBG(a, b) \
	FEAL_ASSERT_DBG(a <= b)(a)(b)

/// Definition of comparative assert macro (Greater Than)
#define ASSERT_GT(a, b) \
	FEAL_ASSERT(a > b)(a)(b)

/// Definition of comparative assert macro (Greater Than) only for debug mode
#define ASSERT_GT_DBG(a, b) \
	FEAL_ASSERT_DBG(a > b)(a)(b)

/// Definition of comparative assert macro (Greater or Equal)
#define ASSERT_GE(a, b) \
	FEAL_ASSERT(a >= b)(a)(b)

/// Definition of comparative assert macro (Greater or Equal) only for debug mode
#define ASSERT_GE_DBG(a, b) \
	FEAL_ASSERT_DBG(a >= b)(a)(b)

/// Definition of comparative assert macro (EQual)
#define ASSERT_EQ(a, b) \
	FEAL_ASSERT(a == b)(a)(b)

/// Definition of comparative assert macro (EQual) only for debug mode
#define ASSERT_EQ_DBG(a, b) \
	FEAL_ASSERT_DBG(a == b)(a)(b)

/// Definition of assert macro checking non-null pointer (PTR)
#define ASSERT_PTR( ptr ) \
	FEAL_ASSERT( (ptr) != nullptr )

/// Definition of assert macro checking non-null pointer (PTR) only for debug mode
#define ASSERT_PTR_DBG( ptr ) \
	FEAL_ASSERT_DBG( (ptr) != nullptr )



/// Allow use shorter versions of macro names if these names is not used with external library
#ifndef ASSERT
#define ASSERT( expr) FEAL_ASSERT( expr)
#endif
#ifndef ASSERT_DBG
#define ASSERT_DBG( expr) FEAL_ASSERT_DBG( expr)
#endif




/**
 * Sources:
 * http://www.drdobbs.com/cpp/enhancing-assertions/184403745
 * https://gist.github.com/hang-qi/5308285
 * https://beliefbox.googlecode.com/svn-history/r825/trunk/src/core/SmartAssert.h
 */

#endif // ASSERTS_HH
