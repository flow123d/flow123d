/*
 * exc_common.hh
 *
 *  Created on: Jan 17, 2014
 *      Author: jb
 */

#ifndef EXC_COMMON_HH_
#define EXC_COMMON_HH_

/**
 * @file
 *
 * Common exceptions.
 */

#include <system/exceptions.hh>

/**
 * Assert exception with an string message.
 */
TYPEDEF_ERR_INFO( EI_Message, std::string);
DECLARE_EXCEPTION( ExcAssertMsg, << "Violated Assert! " << EI_Message::val);

/**
 * General exception with message.
 */
DECLARE_EXCEPTION( ExcMessage, << EI_Message::val);







#endif /* EXC_COMMON_HH_ */
