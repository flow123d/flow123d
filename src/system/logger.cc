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
#include <iomanip>


/*******************************************************************
 * implementation of LoggerOptions
 */

LoggerOptions& LoggerOptions::get_instance() {
	return *instance_;
}


LoggerOptions* LoggerOptions::instance_ = new LoggerOptions();


LoggerOptions::LoggerOptions()
: mpi_rank_(-1), no_log_(false), init_(false) {}


LoggerOptions::~LoggerOptions() {
	if (file_stream_.is_open()) {
		file_stream_ << std::flush;
		file_stream_.close();
	}
}


int LoggerOptions::get_mpi_rank() {
	return mpi_rank_;
}


int LoggerOptions::setup_mpi(MPI_Comm comm) {
	ASSERT(!init_).error("Setup MPI must be performed before setting logger file.");

	return MPI_Comm_rank(comm, &mpi_rank_);
}


void LoggerOptions::set_log_file(std::string log_file_base) {
	ASSERT(!init_).error("Recurrent initialization of logger file stream.");

	if (log_file_base.size() == 0) { // empty string > no_log
		no_log_ = true;
	} else {
		int mpi_rank = LoggerOptions::get_mpi_rank();
		if (mpi_rank == -1) { // MPI is not set, random value is used
			std::random_device rd;
			std::mt19937 gen(rd());
			std::uniform_int_distribution<int> dis(0, 999999);
			mpi_rank = dis(gen);
			WarningOut() << "Unset MPI rank, random value '" << mpi_rank << "' of rank will be used." << std::endl;
		}
		std::stringstream file_name;
		file_name << log_file_base << "." << mpi_rank << ".log";
		file_stream_.open( file_name.str().c_str(), std::ofstream::out );
	}
	init_ = true;
}


void LoggerOptions::reset() {
	mpi_rank_ = -1;
	no_log_ = false;
	init_ = false;
	if (file_stream_.is_open()) {
		file_stream_ << std::flush;
		file_stream_.close();
	}
}


/*******************************************************************
 * implementation of MultiTargetBuf
 */


const std::string MultiTargetBuf::msg_type_string(MsgType msg_type)
{
	switch (msg_type) {
		case MsgType::warning: return "WARNING.";
		case MsgType::message: return "MESSAGE.";
		case MsgType::log:     return "LOG.";
		default:               return "DEBUG.";
	}
}


TimePoint MultiTargetBuf::start_time = TimePoint();


MultiTargetBuf::MultiTargetBuf(MsgType type)
: std::stringbuf(), type_(type), every_process_(false), streams_mask_(0), printed_header_(false)
{
	// set actual time
	TimePoint t = TimePoint();
	date_time_ = TimePoint::format_hh_mm_ss(t-MultiTargetBuf::start_time);

    // set MPI rank
    mpi_rank_ = LoggerOptions::get_instance().get_mpi_rank();
}


int MultiTargetBuf::sync() {
	if (!this->in_avail()) return 0;
	ASSERT_DBG(this->streams_mask_).error("Mask of logger is not set.");

	std::string segment;
	bool first_segment = true;
	std::istringstream istream(str());
	while(std::getline(istream, segment)) {
		segments_.push_back(segment);
	}

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
	for (auto seg : segments_) {
		if (first_segment) {
			formated_output_ << std::setw(4) << "" << "- ";
			first_segment = false;
		} else {
			formated_output_ << std::setw(6) << "";
		}
		formated_output_ << seg << "\n";
	}

	print_to_screen(std::cout, MultiTargetBuf::mask_cout);
	print_to_screen(std::cerr, MultiTargetBuf::mask_cerr);
	if (LoggerOptions::get_instance().is_init())
		print_to_file(LoggerOptions::get_instance().file_stream_, MultiTargetBuf::mask_file);

	// Marks printed header, clean class members
	printed_header_ = true;
	str("");
	segments_.clear();

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
		if (LoggerOptions::get_instance().no_log_)
			streams_mask_ = MultiTargetBuf::mask_cerr;
		else
			streams_mask_ = MultiTargetBuf::mask_cerr | MultiTargetBuf::mask_file;
		break;
	case MsgType::message:
		if (LoggerOptions::get_instance().no_log_)
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
		if (LoggerOptions::get_instance().no_log_)
			streams_mask_ = 0;
		else if (LoggerOptions::get_instance().is_init())
			streams_mask_ = MultiTargetBuf::mask_file;
		else
			streams_mask_ = MultiTargetBuf::mask_cerr;
		break;
	}

}


void MultiTargetBuf::print_to_screen(std::ostream& stream, unsigned int mask)
{
	if (streams_mask_ & mask) {
		bool header_line = false;
		stream << setfill(' ');

		// print header (once time)
		if (!printed_header_) {
			stream << date_time_ << " ";
			if (every_process_) { // rank
				stringstream rank;
				rank << "[" << mpi_rank_ << "]";
				stream << setiosflags(ios::left) << std::setw(5) << rank.str();
			} else {
				stream << std::setw(5) << "";
			}

			if (type_ != MsgType::message) { // type of message (besides message)
				stream << MultiTargetBuf::msg_type_string(type_) << "\n";
			} else {
				header_line = true;
			}
		}

		// print message
		for (auto segment : segments_) {
			if (header_line) {
				header_line = false;
			} else {
				stream << std::setw(18) << "";
			}
			stream << segment << "\n";
		}

		stream << std::flush;
	}
}


void MultiTargetBuf::print_to_file(std::ofstream& stream, unsigned int mask)
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
