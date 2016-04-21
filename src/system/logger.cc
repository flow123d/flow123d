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
 * @file    logger.cc
 * @brief
 */


#include "system/logger.hh"
#include "config.h"

#include <time.h>
#include <mpi.h>


Logger::Logger(MsgType type, bool every_process)
: type_(type), every_process_(every_process)
{
	// set actual time
	time_t     now = time(0);
	char buf[80];
    strftime(buf, sizeof(buf) - 1, "%y.%m.%d_%H-%M-%S", localtime(&now));
    date_time_ = std::string(buf);

    // set MPI rank
#ifdef FLOW123D_HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
#else
    mpi_rank_ = -1;
#endif
}


Logger::~Logger()
{}


Logger& Logger::set_context(const char* file_name, const char* function, const int line)
{
	file_name_ = std::string(file_name);
	function_ = std::string(function);
	line_ = line;

	return *this;
}
