/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file xio.h
 * @brief  I/O functions with filename storing, able to track current line in opened file. All standard
 *         stdio functions working with files (not stdin, stdout, stderr) should be replaced
 *         by their equivalents from XIO library.
 *
 * @date 6.4.2010
 * @author Jiri Jenicek
 *
 */

/*
 *   
 *   Commented lines contains functions of stdio library that are not used in Flow123d.
 *   Their equivalents are currently not implemented in XIO library.
 *
 */

#ifndef XIO_H_
#define XIO_H_

#include "global_defs.h"
#include "system/system.hh"

#include <map>

// Size of buffer for text input lines.
#define LINE_SIZE 65536

//! @brief XFILE structure holds additional info to generic FILE
/// @{
typedef struct xfile {
    char * filename;  ///< file name in the time of opening
    char * mode;      ///< opening mode
    int    lineno;    ///< last read line (only for text files)
} XFILE;
//! @}


/**
 * Base class of XIO library.
 *
 * The class implements a singleton pattern.
 * Stores data of file mapping and debug output for XIO function.
 */
class Xio {
public:
	/// mapping of ptr to regular file structure to extended structure
	typedef map< FILE *, XFILE * > XFILEMAP;

	/// return instance
	static Xio *get_instance();
	/// initialize XIO library
	static void init();

	/// Enable/Disable XIO debug output for EACH XIO function call
	void     set_verbosity(int verb);
	/// Get current XIO debug verbosity level
	int      get_verbosity();
	/// Get XIO mapping instance
	XFILEMAP &get_xfile_map();

private:
	// Singleton instance
	static Xio *instance;

	/**
	 * Constructor.
	 *
	 * Initialize XIO library.
	 */
	Xio();

	/// mapping instance
	XFILEMAP xfiles_map_;
	/// internal XIO debug: print info at each XIO function
	int verbosity_;

};

//! @brief XIO library extensions
/// @{
char * xio_getfname( FILE * f );     ///< Get file name from file stream
char * xio_getfmode( FILE * f );     ///< Get file mode from file stream
int    xio_getlinesread( FILE * f ); ///< Get number of read lines from stream
char * xio_getfulldescription( FILE * f ); ///< Get pointer to string with full file description
//! @}

//! @brief File access
/// @{
FILE * xfopen( const char * fname, const char * mode );   ///< Open file (function)
FILE * xfopen( const std::string &fname, const char * mode );
int    xfclose( FILE * stream );                          ///< Close file (function)
int    xfflush( FILE * f );                               ///< Flush stream (function)
FILE * xfreopen( const char * filename, const char * mode, FILE * stream );
//! @}

//! @brief Formatted input/output
/// @{
int     xfprintf( FILE * out, const char * fmt, ... );  ///< Write formatted output to stream (function)
int     xfscanf( FILE * in, const char * fmt, ... );    ///< Read formatted data from stream (function)
//! @}

//! @brief Character input/output
/// @{
char *  xfgets( char *s, int n, FILE *in );        ///< Get string from stream (function)
int     xfgetc( FILE * f );                        ///< Get character from stream (function)
int     xgetc( FILE * f );                         ///< Get character from stream (function)
int     xungetc( int c, FILE * f );                ///< Unget character from stream (function)

//! @}

//! @brief Direct input/output
/// @{
size_t xfread( void * ptr, size_t size, size_t count, FILE * stream );        ///< Read block of data from stream (function)
size_t xfwrite( const void * ptr, size_t size, size_t count, FILE * stream ); ///< Write block of data to stream (function)
//! @}

//! @brief File positioning
/// @{
void xrewind( FILE * f );     ///< Set position indicator to the beginning (function)
//! @}

//! @brief Error-handling
/// @{
int xfeof ( FILE * f );       ///< Check End-of-File indicator (function)
//! @}

#endif /* XIO_H_ */
