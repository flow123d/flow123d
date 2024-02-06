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
 * @file    logger_options.hh
 * @brief
 */

#ifndef LOGGER_OPTIONS_HH_
#define LOGGER_OPTIONS_HH_

#include <fstream>
#include <mpi.h>
#include <string>   // for string
class TimePoint;



/**
 * Helper class defined logger output file and flags for setting of logger.
 *
 * Use singleton design pattern.
 *
 * Setting of logger is ensured by two methods: setup_mpi and set_log_file. Both methods are
 * optional.
 *
 *  - setup_mpi sets actual rank of procces according to given MPI communicator. If setting
 *    is not performed, rank has default value -1. For correct functioning this method must
 *    be called before set_log_file.
 *
 *  - set_log_file allows set base name of logger output file (prefix). Name of logger file
 *    is created in format '<log_file_base>.<MPI_rank>.log'. If MPI rank is not set, it's
 *    generated random value (this option is not recommended). Method allows turn off logging
 *    if parameter log_file_base is set to empty string. If set_log_file method is not called,
 *    all logger messages are redirected to screen output.
 *
 * Example of complete initialization of logger:
 *
 @code
   std::string log_file_prefix;
   int mpi_rank;
   // ... set value of log_file_prefix
   MPI_Comm_rank(comm, &mpi_rank);
   LoggerOptions::get_instance().setup_mpi(mpi_rank);
   LoggerOptions::get_instance().set_log_file(log_file_prefix);
 @endcode
 */
class LoggerOptions
{
public:
    /// Initialization flag of Logger.
    enum InitFlag {
        uninitialize,
        no_log,
        initialize
    };

    /// Getter of singleton instance object
    static LoggerOptions& get_instance();

    /**
     * Return actual time from the beginning of application runtime in format HH:MM:SS.SSS
     */
    static std::string format_hh_mm_ss();

    /// Returns number of actual process, if MPI is not supported returns -1.
	inline int get_mpi_rank() const {
		return mpi_rank_;
	}

	/// Set rank of actual process.
	void set_mpi_rank(int mpi_rank);

    /// Reset MPI rank and log file name
	void reset();

    /// Check if singleton instance object is initialize.
	inline LoggerOptions::InitFlag init_flag()
	{ return init_flag_; }

	/// Set \p init_ flag.
	void set_stream(std::string abs_path);

	/// Create unique log file name
	std::string log_file_name(std::string log_file_base);

	/// Destructor
	~LoggerOptions();
private:
	/// Forbidden constructor
	LoggerOptions();

	/// Start time of program, allows you to specify the actual time of program (see \p format_hh_mm_ss method)
	static TimePoint start_time;

	/// Singleton instance
	static LoggerOptions* instance_;

	/// Stream for storing logger messages to file.
	std::ofstream file_stream_;

	/**
	 * @brief Actual process number
	 *
	 * Default value is set to -1 and indicates that MPI is not set
	 */
	int mpi_rank_;

	/// Flag sign if logger is initialized
	InitFlag init_flag_;

	friend class Logger;
};


#endif /* LOGGER_OPTIONS_HH_ */
