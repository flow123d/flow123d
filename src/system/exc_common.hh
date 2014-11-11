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
TYPEDEF_ERR_INFO( EI_MPI_Rank, int);
DECLARE_EXCEPTION( ExcAssertMsg, << "[" << EI_MPI_Rank::val << "] "
		                         << "Violated Assert! " << EI_Message::val);

/**
 * General exception with message.
 */
DECLARE_EXCEPTION( ExcMessage, << EI_Message::val);



/**
 * Test of ierr return codes for MPI and PETSc
 */
TYPEDEF_ERR_INFO( EI_ErrCode, int);
DECLARE_EXCEPTION( ExcChkErr, << "[" << EI_ErrCode::val << "] ");
DECLARE_EXCEPTION( ExcChkErrAssert, << "[" << EI_ErrCode::val << "] ");






#endif /* EXC_COMMON_HH_ */
