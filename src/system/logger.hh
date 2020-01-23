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
#include <algorithm>            // for forward
#include <string>               // for string
#include <sstream>
#include <vector>
#include <string>


#include "system/fmt/format.h"
#include "system/exc_common.hh"



/**
 * @brief Helper class, store mask specifying streams
 *
 * Defines masks of all used streams as static methods and allows combining and comparing masks using
 * the overloaded operators.
 */
class StreamMask {
public:
	/// Empty constructor
	StreamMask()
	: mask_(0) {}

	/// Constructor set \p mask_ value
	StreamMask(int mask)
	: mask_(mask) {}

	/// Predefined mask of std::cout output
	static StreamMask cout;

	/// Predefined mask of std::cerr output
	static StreamMask cerr;

	/// Predefined mask of log file output
	static StreamMask log;

	// Overload & operator
	StreamMask operator &(const StreamMask &other);

	// Overload | operator
	StreamMask operator |(const StreamMask &other);

	// Overload () operator
	int operator()(void);

private:
	int mask_;
};


/**
 * @brief Class for storing logger messages.
 *
 * Allow define different levels of log messages and distinguish output streams
 * for individual leves. These output streams are -
 *  - standard console output (std::cout)
 *  - standard error output (std::cerr)
 *  - file output (see \p LoggerOptions class)
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
   MessageOut() << "End of simulation at time: " << secondary_eq->solved_time() << "\n";
   WarningOut() << "Unprocessed keys '" << keys_vec << "'." << "\n";
   LogOut() << "Write output to output stream: " << this->_base_filename << " for time: " << time << "\n";
   DebugOut() << "Calling 'initialize' of empty equation '" << typeid(*this).name() << "'." << "\n";
 @endcode
 *
 * See that output of vectors of printable objects is supported.
 * Implementation of Logger does not support manipulator std::endl. Please, use "\n" instead.
 * New line character "\n" can be used multiple times.
 *
 * Logger allow using fmtlib functionality for simpler formatting of message:
 *
 @code
   MessageOut() << fmt::format("Start time: {}\nEnd time: {}\n", this->start_time(), this->end_time());
   MessageOut().fmt("Start time: {}\nEnd time: {}\n", this->start_time(), this->end_time());
 @endcode
 *
 * Messages are printed only on the zero MPI process by default.
 * Parallel output on all processes must be required explicitly:
 *
 @code
   MessageOut().every_proc() << "Size distributed at process: " << distr->lsize() << "\n";
 @endcode
 *
 *
 */
class Logger : public std::ostream {
public:
	/// Enum of types of Logger messages.
	enum MsgType {
		warning = 0,
		message = 1,
		log = 2,
		debug = 3,
		error = 4
	};

	/// Return string value of given MsgType in full or shorter format (e.g. "WARNING" of "Wrn")
	static const std::string msg_type_string(MsgType msg_type, bool full_format = true);

	/// Constructor.
	Logger(MsgType type);

	/// Stores values for printing out line number, function, etc
	Logger& set_context(const char* file_name, const char* function, const int line);

	/// Set flag every_process_ to true
	Logger& every_proc();

	/**
	 * @brief Allow use functionality of fmtlib for formating message.
	 *
	 * See examples in description of Logger class.
	 */
	template<class... T>
	Logger& fmt(T&&... t)
	{
	    try {
	        return *this << fmt::format(std::forward<T>(t)...);
	    } catch (const fmt::FormatError & e) {
	        THROW(ExcMessage() << EI_Message("FormatError: " + std::string(e.what())));
	    }
	}

	/// Destructor.
	~Logger();

    
private:
	/// Set @p streams_mask_ according to the type of message.
	void set_mask();

	/// Print formated message to given screen stream if mask corresponds with @p streams_mask_.
	void print_to_screen(std::ostream& stream, std::stringstream& scr_stream, StreamMask mask);

	/// Print formated message to given file stream if mask corresponds with @p streams_mask_.
	void print_to_file(std::ofstream& stream, std::stringstream& file_stream, StreamMask mask);

	/// Print header to screen stream, helper method called from \p print_to_screen.
	bool print_screen_header(std::stringstream& scr_stream);

	/// Print header to file stream, helper method called from \p print_to_file.
	void print_file_header(std::stringstream& file_stream);

	/// Return compact (relative) path to the given source file.
	std::string compact_file_name(std::string file_name);

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
	StreamMask streams_mask_;             ///< Mask of logger, specifies streams in actual time into which to be written
	StreamMask full_streams_mask_;        ///< Mask of logger, specifies all streams into which to be written in logger message

	// Generic printing.
	template <class T>
	friend Logger &operator<<(Logger & log, const T & x);
	// Vector printing support.
	template <class T>
	friend Logger &operator<<(Logger & log, const std::vector<T> & vec);
	// Parametric stream modificator (.
	friend Logger &operator<<(Logger & log, std::ostream & (*pf) (std::ostream &) );
	// Stream mask modificator.
	friend Logger &operator<<(Logger & log, StreamMask mask);
};


Logger &operator<<(Logger & log, StreamMask mask);


Logger &operator<<(Logger & log, std::ostream & (*pf) (std::ostream &) );


template <class T>
Logger &operator<<(Logger & log, const std::vector<T> & vec)
{
    for (T const& c : vec)
        log << c << " ";
	return log;
}


template <class T>
Logger &operator<<(Logger & log, const T & x)
{
	if ( (log.streams_mask_ & StreamMask::cout)() ) log.cout_stream_ << x;
	if ( (log.streams_mask_ & StreamMask::cerr)() ) log.cerr_stream_ << x;
	if ( (log.streams_mask_ & StreamMask::log )() ) log.file_stream_ << x;
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

/**
 * Print variable name and value.
 * Usage:
 * DebugOut() << print_var(x) << print_var(y)
 */
#define print_var(var) \
    std::string(#var) << "=" << (var) << ", "




#endif /* LOGGER_HH_ */
