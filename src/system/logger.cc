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


/*******************************************************************
 * implementation of MultiTargetBuf
 */


const std::string MultiTargetBuf::msg_type_string(MsgType msg_type)
{
	switch (msg_type) {
		case _warning: return "Warning";
		case _message: return "Message";
		case _log:     return "Log";
		default:       return "Debug";
	}
}


MultiTargetBuf::MultiTargetBuf(MsgType type, bool every_process)
: std::stringbuf(), type_(type), every_process_(every_process), streams_mask_(0)
{
	// set actual time
	time_t     now = time(0);
	char buf[80];
    strftime(buf, sizeof(buf) - 1, "%b %d %Y %X", localtime(&now));
    date_time_ = std::string(buf);

    // set MPI rank
#ifdef FLOW123D_HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
#else
    mpi_rank_ = -1;
#endif
}


int MultiTargetBuf::sync() {
	print_to_stream(std::cout, MultiTargetBuf::mask_cout);
	print_to_stream(std::cerr, MultiTargetBuf::mask_cerr);

	str("");
	return 0;
}


void MultiTargetBuf::set_context(const char* file_name, const char* function, const int line)
{
	this->set_mask();
	file_name_ = std::string(file_name);
	function_ = std::string(function);
	line_ = line;
}


void MultiTargetBuf::set_mask()
{
	if ( !every_process_ && (mpi_rank_ > 0) ) return;
	if (type_ == _warning) streams_mask_ = MultiTargetBuf::mask_cerr;
	else streams_mask_ = MultiTargetBuf::mask_cout;
}


void MultiTargetBuf::print_to_stream(std::ostream& stream, unsigned int mask)
{
	if (streams_mask_ & mask) {
		stream << str();
		stream << "type : " << MultiTargetBuf::msg_type_string(type_) << "\n";
		stream << "mpi_rank : " << mpi_rank_ << "\n";
		stream << "time : " << date_time_ << "\n";
		stream << "code_point : " << file_name_ << "(" << line_ << ")\n";
		stream << "function : " << function_ << "\n";
		stream << std::flush;
	}
}


/*******************************************************************
 * implementation of Logger
 */


Logger::Logger(MsgType type, bool every_process)
: std::ostream( new MultiTargetBuf(type, every_process) )
{}


Logger::~Logger()
{
	delete rdbuf();
	rdbuf(NULL);
}


Logger& Logger::set_context(const char* file_name, const char* function, const int line)
{
	rdbuf()->pubsync();
	dynamic_cast<MultiTargetBuf *>(rdbuf())->set_context(file_name, function, line);
	return *this;
}

/**
 * Solving of output to multiple stream
 * http://www.cplusplus.com/forum/general/54588/
 */
