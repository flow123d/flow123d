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
#include "system/global_defs.h"
#include "config.h"

#include <time.h>
#include <mpi.h>
#include <iomanip>


/*******************************************************************
 * implementation of LoggerFileStream
 */

LoggerFileStream& LoggerFileStream::get_instance() {
	return *instance_;
}


int LoggerFileStream::get_mpi_rank() {
#ifdef FLOW123D_HAVE_MPI
	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	return mpi_rank;
#else
	return 0;
#endif
}

void LoggerFileStream::init(const std::string &log_file_name, bool no_log) {
	ASSERT(instance_ == nullptr && !no_log).error("Recurrent initialization of logger file stream.");

	if (no_log) {
		LoggerFileStream::no_log_ = true;
	} else {
		std::stringstream file_name;
		file_name << log_file_name << "." << LoggerFileStream::get_mpi_rank() << ".log";
		instance_ = new LoggerFileStream( file_name.str().c_str() );
	}
}


LoggerFileStream* LoggerFileStream::instance_ = nullptr;
bool LoggerFileStream::no_log_ = false;


LoggerFileStream::~LoggerFileStream() {
	(*this) << std::flush;
	this->close();
}


LoggerFileStream::LoggerFileStream()
: std::ofstream() {}


LoggerFileStream::LoggerFileStream(const char* filename)
: std::ofstream(filename) {}


/*******************************************************************
 * implementation of MultiTargetBuf
 */


const std::string MultiTargetBuf::msg_type_string(MsgType msg_type)
{
	switch (msg_type) {
		case MsgType::warning: return "Warning";
		case MsgType::message: return "Message";
		case MsgType::log:     return "Log";
		default:               return "Debug";
	}
}


MultiTargetBuf::MultiTargetBuf(MsgType type)
: std::stringbuf(), type_(type), every_process_(false), streams_mask_(0), printed_header_(false)
{
	// set actual time
	time_t     now = time(0);
	char buf[80];
    strftime(buf, sizeof(buf) - 1, "%b %d %Y %X", localtime(&now));
    date_time_ = std::string(buf);

    // set MPI rank
    mpi_rank_ = LoggerFileStream::get_mpi_rank();
}


int MultiTargetBuf::sync() {
	if (!streams_mask_) return 0;

	std::string segment;
	bool first_segment = true;

	formated_output_.str("");
	formated_output_ << std::setfill(' ');
	if (!printed_header_) {
		formated_output_ << " -  -" << std::setw(8) << "";
		formated_output_ << "{ type : \"" << MultiTargetBuf::msg_type_string(type_) << "\", ";
		formated_output_ << "mpi_rank : \"" << mpi_rank_ << "\", ";
		formated_output_ << "time : \"" << date_time_ << "\",\n";
		formated_output_ << std::setw(15) << "" << "code_point : \"";
		formated_output_ << file_name_ << "(" << line_ << "), " << function_ << "\" }\n";
	}
	std::istringstream istream(str());
	while(std::getline(istream, segment)) {
		if (first_segment) {
			formated_output_ << std::setw(4) << "" << "- ";
			first_segment = false;
		} else {
			formated_output_ << std::setw(6) << "";
		}
		formated_output_ << segment << "\n";
	}

	print_to_stream(std::cout, MultiTargetBuf::mask_cout);
	print_to_stream(std::cerr, MultiTargetBuf::mask_cerr);
	if (LoggerFileStream::get_instance() != nullptr)
		print_to_stream(LoggerFileStream::get_instance(), MultiTargetBuf::mask_file);

	printed_header_ = true;
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


void MultiTargetBuf::every_proc()
{
	every_process_ = true;
	this->set_mask();
}


void MultiTargetBuf::set_mask()
{
	if ( !every_process_ && (mpi_rank_ > 0) ) return;

	switch (type_) {
	case MsgType::warning:
		if (LoggerFileStream::no_log_)
			streams_mask_ = MultiTargetBuf::mask_cerr;
		else
			streams_mask_ = MultiTargetBuf::mask_cerr | MultiTargetBuf::mask_file;
		break;
	case MsgType::message:
		if (LoggerFileStream::no_log_)
			streams_mask_ = MultiTargetBuf::mask_cout;
		else
			streams_mask_ = MultiTargetBuf::mask_cout | MultiTargetBuf::mask_file;
		break;
#ifndef FLOW123D_DEBUG
	case MsgType::debug: // for release build
		streams_mask_ = 0;
		break;
#endif
	default: //MsgType::log + MsgType::debug (only for debug build)
		if (LoggerFileStream::no_log_)
			streams_mask_ = 0;
		else if (LoggerFileStream::is_init())
			streams_mask_ = MultiTargetBuf::mask_file;
		else
			streams_mask_ = MultiTargetBuf::mask_cerr;
		break;
	}

}


void MultiTargetBuf::print_to_stream(std::ostream& stream, unsigned int mask)
{
	if (streams_mask_ & mask) {
		stream << formated_output_.str() << std::flush;
	}
}


/*******************************************************************
 * implementation of Logger
 */


Logger::Logger(MultiTargetBuf::MsgType type)
: std::ostream( new MultiTargetBuf(type) )
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


Logger& Logger::every_proc()
{
	rdbuf()->pubsync();
	dynamic_cast<MultiTargetBuf *>(rdbuf())->every_proc();
	return *this;
}

/**
 * Solving of output to multiple stream
 * http://www.cplusplus.com/forum/general/54588/
 */
