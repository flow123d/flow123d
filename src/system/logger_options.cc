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
 * @file    logger_options.cc
 * @brief
 */


#include <mpi.h>

#include "system/logger_options.hh"
#include "system/logger.hh"
#include "system/global_defs.h"
#include "system/time_point.hh"

using namespace std;

/*******************************************************************
 * implementation of LoggerOptions
 */

LoggerOptions& LoggerOptions::get_instance() {
	return *instance_;
}



std::string LoggerOptions::format_hh_mm_ss() {
	TimePoint t = TimePoint();
	double seconds = t-LoggerOptions::start_time;
	ASSERT(seconds > -numeric_limits<double>::epsilon())(seconds).error("Formating of negative time.");

	unsigned int h,m,s,ms;
	unsigned int full_time = (int)(seconds * 1000); // in first step in miliseconds

	ms = full_time % 1000;
	full_time /= 1000;
	s = full_time % 60;
	full_time /= 60;
	m = full_time % 60;
	h = full_time / 60;

	stringstream ss;
	if (h<10) ss << "0";
	ss << h << ":";
	if (m<10) ss << "0";
	ss << m << ":";
	if (s<10) ss << "0";
	ss << s << ".";
	if (ms<100) ss << "0";
	if (ms<10) ss << "0";
	ss << ms;

	return ss.str();
}


TimePoint LoggerOptions::start_time = TimePoint();


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


void LoggerOptions::set_init() {
    ASSERT(!init_).error("Recurrent initialization of logger file stream.");
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
