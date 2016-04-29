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
#include <sstream>



/// Enum of types of Logger class messages.
typedef enum MsgType {
	_warning, _message, _log, _debug
} MsgType;


/**
 * Helper class
 */
class MultiTargetBuf : public std::stringbuf
{
public:
	/// Return string value of given MsgType
	static const std::string msg_type_string(MsgType msg_type);

	/// Constructor
	MultiTargetBuf(MsgType type, bool every_process);

	/// Stores values for printing out line number, function, etc
	void set_context(const char* file_name, const char* function, const int line);
protected:
	virtual int sync();
private:
	static const unsigned int mask_cout = 0b00000001;
	static const unsigned int mask_cerr = 0b00000010;

	/// Set @p streams_mask_ according to the tzpe of message.
	void set_mask();

	/// Print formated message to given stream if mask corresponds with @p streams_mask_.
	void print_to_stream(std::ostream& stream, unsigned int mask);

	MsgType type_;                        ///< Type of message.
	bool every_process_;                  ///< Flag marked if log message is printing for all processes or only for zero process.
	std::string file_name_;               ///< Actual file.
	std::string function_;                ///< Actual function.
	int line_;                            ///< Actual line.
	std::string date_time_;               ///< Actual date and time.
	int mpi_rank_;                        ///< Actual process (if MPI is supported)
	int streams_mask_;                    ///< Mask of logger, specifies streams
	bool printed_header_;                 ///< Flag marked message header was printed (in first call of sync method)
	std::ostringstream formated_output_;  ///< Helper stream, store message during printout to individual output
};


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

};

/// Macro defining one record of log
#define LOG(...) \
	Logger( __VA_ARGS__ ).set_context( __FILE__, __func__, __LINE__)


#endif /* LOGGER_HH_ */
