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


/// Helper function, use for shorten the code point path
std::string cmn_prefix( std::string a, std::string b ) {
    if( a.size() > b.size() ) std::swap(a,b) ;
    return std::string( a.begin(), std::mismatch( a.begin(), a.end(), b.begin() ).first ) ;
}


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
			WarningOut() << "Unset MPI rank, random value '" << mpi_rank << "' of rank will be used.\n";
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
 * implementation of Logger
 */


Logger::Logger(MsgType type)
: type_(type), every_process_(false), streams_mask_(0)
{
	// set actual time
	TimePoint t = TimePoint();
	date_time_ = TimePoint::format_hh_mm_ss(t-Logger::start_time);

    // set MPI rank
    mpi_rank_ = LoggerOptions::get_instance().get_mpi_rank();
}


Logger::~Logger()
{
	// print output to streams
	print_to_screen(std::cout, cout_stream_, Logger::cout_mask);
	print_to_screen(std::cerr, cerr_stream_, Logger::cerr_mask);
	if (LoggerOptions::get_instance().is_init())
		print_to_file(LoggerOptions::get_instance().file_stream_, Logger::file_mask);
}


const std::string Logger::msg_type_string(MsgType msg_type, bool full_format)
{
	if (full_format) {
		switch (msg_type) {
			case MsgType::warning: return "WARNING.";
			case MsgType::message: return "MESSAGE.";
			case MsgType::log:     return "LOG.";
			default:               return "DEBUG.";
		}
	} else {
		switch (msg_type) {
			case MsgType::warning: return "Wrn";
			case MsgType::message: return "Msg";
			case MsgType::log:     return "Log";
			default:               return "Dbg";
		}
	}
}


TimePoint Logger::start_time = TimePoint();


Logger& Logger::set_context(const char* file_name, const char* function, const int line)
{
	file_name_ = std::string(file_name);
	function_ = std::string(function);
	line_ = line;
	this->set_mask();

	return *this;
}


Logger& Logger::every_proc()
{
	every_process_ = true;
	this->set_mask();

	return *this;
}

void Logger::set_mask()
{
	if ( !every_process_ && (mpi_rank_ > 0) ) return;

	switch (type_) {
	case MsgType::warning:
		if (LoggerOptions::get_instance().no_log_)
			streams_mask_ = Logger::cerr_mask;
		else
			streams_mask_ = Logger::cerr_mask | Logger::file_mask;
		break;
	case MsgType::message:
		if (LoggerOptions::get_instance().no_log_)
			streams_mask_ = Logger::cout_mask;
		else
			streams_mask_ = Logger::cout_mask | Logger::file_mask;
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
			streams_mask_ = Logger::file_mask;
		else
			streams_mask_ = Logger::cerr_mask;
		break;
	}

}


void Logger::print_to_screen(std::ostream& stream, std::stringstream& scr_stream, unsigned int mask)
{
	if (streams_mask_ & mask) {
		bool header_line = false;
		stream << setfill(' ');

		// print header
		stream << date_time_ << " ";
		if (every_process_) { // rank
			stringstream rank;
			rank << "[" << mpi_rank_ << "]";
			stream << setiosflags(ios::left) << std::setw(5) << rank.str();
		} else {
			stream << std::setw(5) << "";
		}

		if (type_ != MsgType::message) { // type of message (besides Message)
			stream << msg_type_string(type_) << "\n";
		} else {
			header_line = true;
		}

		// print message
		std::string segment;
		std::istringstream istream(scr_stream.str());
		while(std::getline(istream, segment)) {
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


void Logger::print_to_file(std::ofstream& stream, unsigned int mask)
{
	if (streams_mask_ & mask) {
		stream << setfill(' ');

		// print header
		stream << "- -" << std::setw(13) << "" << "[ ";
		stream << msg_type_string(type_, false);
		if (every_process_) { // add 'E' (only for every_proc) + print rank
			stream << "E, ";
			if (mpi_rank_ >= 0) stream << setiosflags(ios::right) << std::setw(4) << mpi_rank_;
			else stream << "null";
		} else {
			stream << " , null";
		}
		stream << ", \"" << date_time_ << "\"";
	    // if constant FLOW123D_SOURCE_DIR is defined, we try to erase it from beginning of each CodePoint's filepath
	    #ifdef FLOW123D_SOURCE_DIR
	        string common_path = cmn_prefix(string(FLOW123D_SOURCE_DIR), file_name_);
	        file_name_.erase (0, common_path.size());
	    #endif
		stream << ", \"" << file_name_ << "\", " << line_ << ", \"" << function_ << "\"";
		stream << " ]\n";

		// print message
		std::string segment;
		std::istringstream istream(file_stream_.str());
		stream << "  - |" << "\n";
		while(std::getline(istream, segment)) {
			stream << std::setw(4) << "" << segment << "\n";
		}

		stream << std::flush;
	}
}
