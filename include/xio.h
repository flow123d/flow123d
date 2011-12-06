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

#include <cstdio>
#include "system/system.hh"

//! @brief XIO library extensions
/// @{
void   xio_init( void );             ///< XIO library initialization
char * xio_getfname( FILE * f );     ///< Get file name from file stream
char * xio_getfmode( FILE * f );     ///< Get file mode from file stream
int    xio_getlinesread( FILE * f ); ///< Get number of read lines from stream
char * xio_getfulldescription( FILE * f ); ///< Get pointer to string with full file description
void   xio_setverbose( int verb );   ///< Enable/Disable XIO debug output for EACH XIO function call
int    xio_getverbose( void );       ///< Get current XIO debug verbosity level
//! @}

//! @brief File access
/// @{
FILE * xfopen( const char * fname, const char * mode );   ///< Open file (function)
FILE * xfopen( const std::string &fname, const char * mode );
int    xfclose( FILE * stream );                          ///< Close file (function)
int    xfflush( FILE * f );                               ///< Flush stream (function)
FILE * xfreopen( const char * filename, const char * mode, FILE * stream );

//#define setbuf xsetbuf     NOT_USED ///< Set stream buffer (function)
//#define setvbuf xsetvbuf   NOT_USED ///< Change stream buffering (function)
//! @}

//! @brief Formatted input/output
/// @{
int     xfprintf( FILE * out, const char * fmt, ... );  ///< Write formatted output to stream (function)
int     xfscanf( FILE * in, const char * fmt, ... );    ///< Read formatted data from stream (function)

//#define vfprintf xvfprintf  NOT_USED ///< Write formatted variable argument list to stream (function)
//! @}

//! @brief Character input/output
/// @{
char *  xfgets( char *s, int n, FILE *in );        ///< Get string from stream (function)
int     xfgetc( FILE * f );                        ///< Get character from stream (function)
int     xgetc( FILE * f );                         ///< Get character from stream (function)
int     xungetc( int c, FILE * f );                ///< Unget character from stream (function)

//#define fputc xfputc            NOT_USED        ///< Write character to stream (function)
//#define fputs xfputs            NOT_USED        ///< Write string to stream (function)
//#define putc xputc              NOT_USED        ///< Write character to stream (function)

//! @}

//! @brief Direct input/output
/// @{
size_t xfread( void * ptr, size_t size, size_t count, FILE * stream );        ///< Read block of data from stream (function)
size_t xfwrite( const void * ptr, size_t size, size_t count, FILE * stream ); ///< Write block of data to stream (function)
//! @}

//! @brief File positioning
/// @{
void xrewind( FILE * f );     ///< Set position indicator to the beginning (function)

//#define fgetpos xfgetpos   NOT_USED ///< Get current position in stream (function)
//#define fseek xfseek       NOT_USED ///< Reposition stream position indicator (function)
//#define fsetpos xfsetpos   NOT_USED ///< Set position indicator of stream (function)
//#define ftell xftell       NOT_USED ///< Get current position in stream (function)
//! @}

//! @brief Error-handling
/// @{
int xfeof ( FILE * f );       ///< Check End-of-File indicator (function)

//#define clearerr xclearerr NOT_USED ///< Clear error indicators (function)
//#define ferror xferror     NOT_USED ///< Check error indicator (function)
//#define perror xperror     NOT_USED ///< Print error message (function)
//! @}

#endif /* XIO_H_ */
