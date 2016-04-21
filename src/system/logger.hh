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
 * @file    logger.hh
 * @brief
 */

#ifndef LOGGER_HH_
#define LOGGER_HH_


#include <iostream>



/// Enum of types of Logger class messages.
typedef enum MsgType {
	_warning, _message, _log, _debug
} MsgType;


/**
 * @brief Class for storing log messages.
 *
 */
class Logger : public std::ostream {
public:
	/// Constructor.
	Logger(MsgType type, bool every_process = false);

	/// Stores values for printing out line number, function, etc
	Logger& set_context(const char* file_name, const char* function, const int line);

	/// Destructor.
	~Logger();

private:
	MsgType type_;                 ///< Type of message.
	bool every_process_;           ///< Flag marked if log message is printing for all processes or only for zero process.
	std::string file_name_;        ///< Actual file.
	std::string function_;         ///< Actual function.
	int line_;                     ///< Actual line.
	std::string date_time_;        ///< Actual date and time.
	int mpi_rank_;                 ///< Actual process (if MPI is supported)
};

/// Macro defining one record of log
#define LOG(...) \
	Logger( __VA_ARGS__ ).set_context( __FILE__, __func__, __LINE__)


#endif /* LOGGER_HH_ */
