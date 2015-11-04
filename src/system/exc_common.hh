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
 * @file    exc_common.hh
 * @brief   
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
 * Usage:
 * THROW( ExcMessage() << EI_Message("Some message.") )
 */
DECLARE_EXCEPTION( ExcMessage, << EI_Message::val);



/**
 * Test of ierr return codes for MPI and PETSc
 */
TYPEDEF_ERR_INFO( EI_ErrCode, int);
DECLARE_EXCEPTION( ExcChkErr, << "[" << EI_ErrCode::val << "] ");
DECLARE_EXCEPTION( ExcChkErrAssert, << "[" << EI_ErrCode::val << "] ");






#endif /* EXC_COMMON_HH_ */
