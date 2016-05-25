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
#include <fstream>
#include <vector>
#include <mpi.h>
#include "system/time_point.hh"



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
   // ... set value of log_file_prefix
   LoggerOptions::get_instance().setup_mpi(MPI_COMM_WORLD);
   LoggerOptions::get_instance().set_log_file(log_file_prefix);
 @endcode
 */
class LoggerOptions
{
public:
    /// Getter of singleton instance object
	static LoggerOptions& get_instance();

    /// Returns number of actual process, if MPI is not supported returns -1.
	int get_mpi_rank();

	/// Set rank of actual process by MPI communicator.
	int setup_mpi(MPI_Comm comm);

    /// Initialize instance object in format 'log_file_base.process.log'.
	void set_log_file(std::string log_file_base);

    /// Reset MPI rank and log file name
	void reset();

    /// Check if singleton instance object is initialize.
	inline bool is_init()
	{ return init_; }

	/// Destructor
	~LoggerOptions();
private:
	/// Forbidden constructor
	LoggerOptions();

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

	/// Turn off logger file output
	bool no_log_;

	/// Flag sign if logger is initialized by set_log_file method
	bool init_;

	friend class Logger;
};


/**
 * @brief Class for storing logger messages.
 *
 * Allow define different levels of log messages and distinguish output streams
 * for individual leves. These output streams are -
 *  - standard console output (std::cout)
 *  - standard error output (std::cerr)
 *  - file output (LoggerOptions class)
 *
 * Logger distinguishes four type of levels -
 *  - warning: printed to standard error and file output
 *  - message: printed to standard console and file output
 *  - log: printed to file output
 *  - debug: printed to file output (level is used only in debug mode)
 *
 * File output is optional. See \p LoggerOptions for setting this output stream.
 *
 * <b>Example of Logger usage:</b>
 *
 * For individual levels are defined macros -
 *  - MessageOut()
 *  - WarningOut()
 *  - LogOut()
 *  - DebugOut()
 * that ensure display of actual code point (source file, line and function).
 *
 * Logger message is created by using an operator << and allow to add each type
 * that has override this operator. Message is terminated with manipulator
 * std::endl. Implicitly logger message is printed only in processor with rank
 * zero. If necessary printed message for all process, it provides a method
 * every_proc().
 *
 * Examples of logger messages formating:
 *
 @code
   MessageOut() << "End of simulation at time: " << secondary_eq->solved_time() << std::endl;
   WarningOut() << "Unprocessed key '" << key_name << "' in Record '" << rec->type_name() << "'." << std::endl;
   LogOut() << "Write output to output stream: " << this->_base_filename << " for time: " << time << std::endl;
   DebugOut() << "Calling 'initialize' of empty equation '" << typeid(*this).name() << "'." << std::endl;
 @endcode
 *
 * Logger message can be created by more than one separate message (std::endl
 * manipulator can be used multiple times):
 *
 @code
   MessageOut() << "Start time: " << this->start_time() << std::endl << "End time: " << this->end_time() << std::endl;
 @endcode
 *
 * In some cases message can be printed for all processes:
 *
 @code
   MessageOut().every_proc() << "Size distributed at process: " << distr->lsize() << std::endl;
 @endcode
 *
 */
class Logger : public std::ostream {
public:
	/// Enum of types of Logger messages.
	enum MsgType {
		warning = 0,
		message = 1,
		log = 2,
		debug = 3
	};

	/// Return string value of given MsgType in full or shorter format (e.g. "WARNING" of "Wrn")
	static const std::string msg_type_string(MsgType msg_type, bool full_format = true);

	/// Constructor.
	Logger(MsgType type);

	/// Stores values for printing out line number, function, etc
	Logger& set_context(const char* file_name, const char* function, const int line);

	/// Set flag every_process_ to true
	Logger& every_proc();

	/// Destructor.
	~Logger();

    // treat manipulators
	Logger & operator<<(Logger & (*pf) (Logger &) )
    {
        return pf(*this);
    }

private:
	static const unsigned int cout_mask = 0b00000001;
	static const unsigned int cerr_mask = 0b00000010;
	static const unsigned int file_mask = 0b00000100;
	static TimePoint start_time;

	/// Set @p streams_mask_ according to the type of message.
	void set_mask();

	/// Print formated message to given screen stream if mask corresponds with @p streams_mask_.
	void print_to_screen(std::ostream& stream, std::stringstream& scr_stream, unsigned int mask);

	/// Print formated message to given file stream if mask corresponds with @p streams_mask_.
	void print_to_file(std::ofstream& stream, unsigned int mask);

	std::stringstream cout_stream_;       ///< Store messages printed to cout output stream
	std::stringstream cerr_stream_;       ///< Store messages printed to cerr output stream
	std::stringstream file_stream_;       ///< Store messages printed to file

	MsgType type_;                        ///< Type of message.
	bool every_process_;                  ///< Flag marked if log message is printing for all processes or only for zero process.
	std::string file_name_;               ///< Actual file.
	std::string function_;                ///< Actual function.
	int line_;                            ///< Actual line.
	std::string date_time_;               ///< Actual date and time.
	int mpi_rank_;                        ///< Actual process (if MPI is supported)
	int streams_mask_;                    ///< Mask of logger, specifies streams

	template <class T>
	friend Logger &operator<<(Logger & log, const T & x);
};


template <class T>
Logger &operator<<(Logger & log, const T & x)
{
	if (log.streams_mask_ & Logger::cout_mask) log.cout_stream_ << x;
	if (log.streams_mask_ & Logger::cerr_mask) log.cerr_stream_ << x;
	if (log.streams_mask_ & Logger::file_mask) log.file_stream_ << x;
    return log;
}



/// Internal macro defining universal record of log
#define _LOG(type) \
	Logger( type ).set_context( __FILE__, __func__, __LINE__)
/// Macro defining 'message' record of log
#define MessageOut() \
	_LOG( Logger::MsgType::message )
/// Macro defining 'warning' record of log
#define WarningOut() \
	_LOG( Logger::MsgType::warning )
/// Macro defining 'log' record of log
#define LogOut() \
	_LOG( Logger::MsgType::log )
/// Macro defining 'debug' record of log
#define DebugOut() \
	_LOG( Logger::MsgType::debug )




#endif /* LOGGER_HH_ */
