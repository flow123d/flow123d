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
 * @file    input_exception.hh
 * @brief   
 */

#ifndef SRC_INPUT_EXCEPTION_HH_
#define SRC_INPUT_EXCEPTION_HH_

#include "system/exceptions.hh"
#include "system/exc_common.hh"

namespace Input {

/**
 * @brief Base of exceptions due to user input.
 *
 * Base class for "input exceptions" that are exceptions caused by incorrect input from the user
 * not by an internal error.
 *
 * @ingroup exceptions
 */
class Exception : public virtual ExceptionBase
{
public:
    const char * what () const throw ();
    virtual ~Exception() throw () {};
};


/**
 *  Declaration of error info class for passing Input::Address through exceptions.
 *  Is returned by input accessors : Input::Record, Input::Array, etc.
 *
 *  Use case example:
 *  Input::Record input = ...;
 *  string name=input.val("name");
 *  if (name.size() > STR_LIMIT) THROW(ExcToLongStr() << EI_Address( input.address_string() ));
 *
 *  TODO: if Address class is persistent (every copy is self contented, we can use Address instead of std::string.
 *  see also ei_address methods.
 */
TYPEDEF_ERR_INFO( EI_Address, const std::string);


/**
 * @brief Macro for simple definition of input exceptions.
 *
 * Works in the same way as @p DECLARE_EXCEPTION, just define class derived from
 * @p InputException. Meant to be used for exceptions due to wrong input from user.
 *
 * Reports input address provided through EI_Address object, see above.
 *
 * @ingroup exceptions
 */
#define DECLARE_INPUT_EXCEPTION( ExcName, Format)                             \
struct ExcName : public virtual ::Input::Exception {                          \
     virtual void print_info(std::ostringstream &out) const {                 \
         using namespace internal;                                            \
         ::internal::ExcStream estream(out, *this);                           \
         estream Format                                                       \
                  << "\nAt input address: "                                 \
                  << ::Input::EI_Address::val;                              \
         out << std::endl;                                                  \
     }                                                                      \
     virtual ~ExcName() throw () {}                                         \
}

/**
 * Simple input exception that accepts just string message.
 */
DECLARE_INPUT_EXCEPTION(ExcInputMessage, << EI_Message::val );


} // namespace Input

#endif /* SRC_INPUT_EXCEPTION_HH_ */

