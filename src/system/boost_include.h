/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id: global_defs.h 1888 2012-10-04 19:29:53Z jan.brezina $
 * $Revision: 1888 $
 * $LastChangedBy: jan.brezina $
 * $LastChangedDate: 2012-10-04 21:29:53 +0200 (Čt, 04 říj 2012) $
 *
 * @file   global_defs.h
 * @brief  Global macros to enhance readability and debugging, general constants.
 *
 */

#ifndef BOOST_INCLUDE_H
#define BOOST_INCLUDE_H

/// @}

/**
 *  @brief counting reference pointer - shared_ptr
 *
 *  Allow global usage of shared_ptr. This is replacement of raw pointer to some type(e.g. Distribution *).
 *  It uses counter for number of pointers pointing to the allocated chunk and call destructor and deallocation
 *  when counter is zero. It introduce small memory overhead and time overhead on copy of pointers. Should be used
 *  for objects with moderate number of instances especially if several objects has to share pointer to the same object
 *  (Mesh, Distribution, LocalToGlobalMap) but there is no object naturally responsible for its deallocation or
 *  such an object has shorter life then the object it has created.
 *
 *  Usage is as follows:
 *
 *  @code
 *  boost::shared_ptr<Distribution> my_distr;
 *  my_distr = boost::make_shared<Distribution>(.. parameters of constructor ..)
 *
 *  // use it as any other pointer
 *  my_distr->lsize()
 *
 *  @endcode
 *
 *  Notes:
 *  - Do not forget to specify boost::... to avoid possible conflict with std::shared_ptr if C++11 standard is used.
 *  - The template make_shared could have problems to resolve particular constructor if the types of parameters do not match
 *    exactly, use explicit type cast in such a case.
 */
#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/smart_ptr/make_shared.hpp>


//using boost::shared_ptr;
//using boost::make_shared;




#endif // BOOST_INCLUDE_H
