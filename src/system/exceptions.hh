/*
 * exceptions.hh
 *
 *  Created on: Apr 10, 2012
 *      Author: jb
 */

/**
 * @file Basic exceptions used in Flow123d.
 *
 * We are using boost::exceptions .
 */

#include <boost/exception/all.hpp>

/**
 * Basic exception for Flow123d's exceptions. When deriving particular exceptions always use virtual inheritance:
 * @code
 *      struct my_exception : virtual flow_excepiton {};
 * @endcode
 */
struct FlowException : virtual std::exception, virtual boost::exception { };

/**
 * @brief Macro to simplify declaration of error_info types.
 *
 * These are  used to pass data through boost exceptions from the throw point to the catch point, or possibly collect
 * various data along stack rewinding when an exception is thrown. Example of usage:
 *
 * @code
 * TYPEDEF_ERR_INFO( ErrorCode, int) // declares type ErrorCode_EI
 *
 * ...
 *
 * throw SomeException() << ErrorCode_EI(10); // here you pass 'int'
 *
 * ...
 *
 * catch (SomeException & exception) {
 *      int const * = boost::get_error_info<ErrorCode_EI>(); // here you get pointer to const 'int'
 * }
 * @endcode
 */
#define TYPEDEF_ERR_INFO(tag, type)    typedef boost::error_info< struct tag, type > tag##_EI


#ifndef EXCEPTIONS_HH_
#define EXCEPTIONS_HH_




#endif /* EXCEPTIONS_HH_ */
