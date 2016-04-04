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
 * @file    stack_trace.hh
 * @brief
 */

#ifndef STACK_TRACE_HH_
#define STACK_TRACE_HH_

#include <iostream>
#include <cstring>
#include <vector>


/**
 * @brief Class representing stacktrace of exceptions.
 *
 * Store array of backtrace frames and size of this array, method @p print allows to print
 * formated stacktrace into given stream.
 *
 * @ingroup exceptions
 */
class StackTrace
{
public:
	/// Default constructor, fill stacktrace
	StackTrace();

	/// Copy constructor
	StackTrace(const StackTrace &other);

	/// Destructor
	~StackTrace();

    /// Prints formated stacktrace into given stream @p out.
    void print(std::ostream &out, std::vector<std::string> frames_to_cut = std::vector<std::string>()) const;

private:

    /// Array of backtrace frames returned by glibc backtrace_symbols.
    char ** frames_;

    /// Size of stacktrace table - number of frames.
    int n_frames_;
};


#endif /* STACK_TRACE_HH_ */

