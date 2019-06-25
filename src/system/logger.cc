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
#include "system/logger_options.hh"
#include "system/global_defs.h"
#include "system/file_path.hh"
#include "config.h"

#include <iomanip>


/// Helper function, use for shorten the code point path
std::string cmn_prefix( std::string a, std::string b ) {
    if( a.size() > b.size() ) std::swap(a,b) ;
    return std::string( a.begin(), std::mismatch( a.begin(), a.end(), b.begin() ).first ) ;
}


/*******************************************************************
 * implementation of StreamMask
 */

StreamMask StreamMask::cout = StreamMask(0b00000001);
StreamMask StreamMask::cerr = StreamMask(0b00000010);
StreamMask StreamMask::log  = StreamMask(0b00000100);


StreamMask StreamMask::operator &(const StreamMask &other)
{
   return StreamMask(this->mask_ & other.mask_);
}

StreamMask StreamMask::operator |(const StreamMask &other)
{
   return StreamMask(this->mask_ | other.mask_);
}

int StreamMask::operator()(void)
{
	return mask_;
}


/*******************************************************************
 * implementation of Logger
 */


Logger::Logger(MsgType type)
: type_(type), every_process_(false), line_(0)
{
	// set actual time
	date_time_ = LoggerOptions::format_hh_mm_ss();

    // set MPI rank
    mpi_rank_ = LoggerOptions::get_instance().get_mpi_rank();
}


Logger::~Logger()
{
	// print output to streams
	print_to_screen(std::cout, cout_stream_, StreamMask::cout);
	print_to_screen(std::cerr, cerr_stream_, StreamMask::cerr);
	if (LoggerOptions::get_instance().is_init())
		print_to_file(LoggerOptions::get_instance().file_stream_, this->file_stream_, StreamMask::log);
}


const std::string Logger::msg_type_string(MsgType msg_type, bool full_format)
{
	static std::vector<std::string> type_names = {"WARNING.", "MESSAGE.", "LOG.", "DEBUG.", "ERROR",
			                                      "Wrn", "Msg", "Log", "Dbg", "Err"};
	int type_idx = msg_type;

	if (full_format) {
		return type_names[type_idx];
	} else {
		return type_names[type_idx + 5];
	}
}


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
	case MsgType::error:
		if (LoggerOptions::get_instance().no_log_)
			streams_mask_ = StreamMask::cerr;
		else
			streams_mask_ = StreamMask::cerr | StreamMask::log;
		break;
	case MsgType::message:
		if (LoggerOptions::get_instance().no_log_)
			streams_mask_ = StreamMask::cout;
		else
			streams_mask_ = StreamMask::cout | StreamMask::log;
		break;
	case MsgType::log:
        if (LoggerOptions::get_instance().no_log_)
            streams_mask_ = StreamMask();
        else if (LoggerOptions::get_instance().is_init())
            streams_mask_ = StreamMask::log;
        else
            streams_mask_ = StreamMask::cout;
        break;
#ifdef FLOW123D_DEBUG
    case MsgType::debug: // for debug build
        if (LoggerOptions::get_instance().no_log_)
            streams_mask_ = StreamMask();
        else if (LoggerOptions::get_instance().is_init())
            streams_mask_ = StreamMask::log | StreamMask::cout;
        else
            streams_mask_ = StreamMask::cout;
        break;
#else
	case MsgType::debug: // for release build
		streams_mask_ = StreamMask();
		break;

#endif
	default:
	    ASSERT(false);
	}

	full_streams_mask_ = full_streams_mask_ | streams_mask_;

}



void Logger::print_to_screen(std::ostream& stream, std::stringstream& scr_stream, StreamMask mask)
{
	if ( full_streams_mask_() & mask() ) {
		stream << setfill(' ');

		// print header, if method returns true, message continues on the same line and first line of message
		// doesn't need indentation in following while cycle
		std::stringstream message_stream;
		bool header_line = this->print_screen_header(message_stream);

		// print message
		std::string segment;
		std::istringstream istream(scr_stream.str());
		while(std::getline(istream, segment)) {
			if (header_line) {
				header_line = false;
			} else {
			    message_stream << std::setw(18) << "";
			}
			message_stream << segment << "\n";
		}

		stream << message_stream.str();
		stream << std::flush;
	}
}


void Logger::print_to_file(std::ofstream& stream, std::stringstream& file_stream, StreamMask mask)
{
	if ( full_streams_mask_() & mask() ) {
		stream << setfill(' ');

		// print header
		std::stringstream message_stream;
		this->print_file_header(message_stream);

		// print message
		std::string segment;
		std::vector<std::string> segments;
		std::istringstream istream(file_stream.str());
		while(std::getline(istream, segment)) {
			segments.push_back(segment);
		}
		if (segments.size() > 1) {
		    message_stream << "  - |" << "\n";
			for (auto seg : segments)
			    message_stream << std::setw(4) << "" << seg << "\n";
		} else if (segments.size() == 1) {
		    message_stream << "  - " << segments[0] << "\n";
		}

		stream << message_stream.str();
		stream << std::flush;
	}
}


std::string Logger::compact_file_name(std::string file_name)
{
    // if constant FLOW123D_SOURCE_DIR is defined, we try to erase it from beginning of each CodePoint's filepath
    #ifdef FLOW123D_SOURCE_DIR
        string common_path = cmn_prefix(string(FLOW123D_SOURCE_DIR), file_name);
        file_name.erase (0, common_path.size());
    #endif
    return file_name;
}

bool Logger::print_screen_header(std::stringstream& stream)
{
	stream << date_time_ << " ";
	if (every_process_) { // rank
		stringstream rank;
		rank << "[" << mpi_rank_ << "]";
		stream << setiosflags(ios::left) << std::setw(5) << rank.str();
	} else {
		stream << std::setw(5) << "";
	}

    stream << msg_type_string(type_);
    if (type_ == MsgType::debug) {
        stream << "(" << compact_file_name(file_name_) << ":" << line_ << ")";
    }
	return true;
}


void Logger::print_file_header(std::stringstream& stream)
{
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
	stream << ", \"" << compact_file_name(file_name_) << "\", " << line_ << ", \"" << function_ << "\"";
	stream << " ]\n";
}



/**
 * implementation of operators
 */

Logger &operator<<(Logger & log, StreamMask mask)
{
	// set mask
	log.streams_mask_ = mask;
	log.full_streams_mask_ = log.full_streams_mask_ | log.streams_mask_;

	return log;
}


Logger &operator<<(Logger & log, std::ostream & (*pf) (std::ostream &) )
{
    if ( (log.streams_mask_ & StreamMask::cout)() ) pf(log.cout_stream_);
    if ( (log.streams_mask_ & StreamMask::cerr)() ) pf(log.cerr_stream_);
    if ( (log.streams_mask_ & StreamMask::log)() ) pf(log.file_stream_);
    return log;
}
