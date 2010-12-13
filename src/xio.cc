/*!
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file xio.cc
 * @brief  I/O functions with filename storing, able to track current line in opened file. All standard
 *         stdio functions working with files (not stdin, stdout, stderr) should be replaced
 *         by their equivalents from XIO library.
 *
 *         Created on: 6.4.2010
 *         Author: Jiri Jenicek
 *
 */

//TODO Better error handling ( perror()? strerror()? )

#include <string.h>
#include <strings.h>
#include <stdarg.h>
#include <errno.h>
#include <limits.h>

#include <iostream>
#include <map>
#include <algorithm>
#include <iterator>

#include "xio.h"

using namespace std;

//! @brief XFILE structure holds additional info to generic FILE
/// @{
typedef struct xfile {
    char * filename;  ///< file name in the time of opening
    char * mode;      ///< opening mode
    int    lineno;    ///< last read line (only for text files)
} XFILE;
//! @}

typedef map< FILE *, XFILE * > XFILEMAP; ///< mapping of ptr to regular file structure to extended structure

static XFILEMAP xfiles_map;    ///< mapping instance
static int xio_verbosity = 0;  ///< internal XIO debug: print info at each XIO function

//! @brief basic definitions
/// @{
static XFILE xstdin  = {strdup("stdin"),strdup("r"),0};
static XFILE xstdout = {strdup("stdout"),strdup("w"),0};
static XFILE xstderr = {strdup("stderr"),strdup("w"),0};
//! @}

static XFILE * xio_getfptr( FILE * f );

#define XIO_WARN(f) xprintf(Warn, "File pointer '%p' not in xfiles_map. Opened with regular fopen() or already closed?\n", (f) )
#define XIO_PRINT_INFO(f) printf( "XIO: In function '%s', %s\n", __func__, xio_getfulldescription( f ) )
#define XIO_DEBUG(f) do { if ( xio_verbosity > 0 ) XIO_PRINT_INFO(f); } while (0)

/*!
 * @brief XIO library initialization
 */
void xio_init( void )
{
    xfiles_map[stdin]  = &xstdin;
    xfiles_map[stdout] = &xstdout;
    xfiles_map[stderr] = &xstderr;
}

/*!
 * @brief Get file name from pointer to FILE structure.
 * @param[in] f pointer to FILE structure
 * @return pointer to file name if OK, NULL if file stream is not known
 */
char * xio_getfname( FILE * f )
{
    XFILE * xf;
    char * rs = NULL;

    xf = xio_getfptr(f);
    if ( xf )
    {
        rs = xf->filename;
    }
    else
    {
        XIO_WARN(f);
    }

    return rs;
}

/*!
 * @brief Get file mode from file stream
 * @param[in] f pointer to FILE structure
 * @return pointer to file opening mode if OK, NULL if file stream is not known
 */
char * xio_getfmode( FILE * f )
{
    XFILE * xf;
    char * rs = NULL;

    xf = xio_getfptr(f);
    if ( xf )
    {
        rs = xf->mode;
    }
    else
    {
        XIO_WARN(f);
    }

    return rs;
}

/*!
 * @brief Get number of lines that were completely read from file since fopen() or rewind()
 * @param[in] f pointer to FILE structure
 * @return number of lines read if OK, -1 if file stream is not known
 */
int xio_getlinesread( FILE * f )
{
    XFILE * xf;
    int lines = -1;

    xf = xio_getfptr(f);
    if ( xf )
    {
        lines = xf->lineno;
    }
    else
    {
        XIO_WARN(f);
    }

    return lines;
}

/*!
 * @brief Get pointer to string with full file description
 * @param[in] f pointer to FILE structure
 * @return pointer to string with description, null terminated, no LF
 */
char * xio_getfulldescription( FILE * f )
{
    const char prn_format_long[] = "FILE ptr %p: x->name '%s', x->mode '%s', x->line '%d'";
    const char prn_format_short[] = "FILE ptr %p: unknown FILE pointer";
    static char * rs = NULL;
    static int maxlen = 0;
    int len;
    XFILE * xf;

    if (!rs)
    {
        maxlen = LINE_SIZE;
        rs = (char *)xmalloc( maxlen );
    }

    xf = xio_getfptr(f);
    if ( xf )
    {
        len = snprintf( rs, maxlen, prn_format_long, f, xf->filename, xf->mode, xf->lineno );
        if ( len >= maxlen ) //string overflow?
        {
           maxlen = len + 1;
           rs = (char *)xrealloc( rs, maxlen );
           snprintf( rs, maxlen, prn_format_long, f, xf->filename, xf->mode, xf->lineno );
        }
    }
    else
    {
        len = snprintf( rs, maxlen, prn_format_short, f );
        if ( len >= maxlen ) //string overflow?
        {
           maxlen = len + 1;
           rs = (char *)xrealloc( rs, maxlen );
           snprintf( rs, maxlen, prn_format_short, f );
        }
    }

    return rs;
}

/*!
 * @brief Internal XIO locator
 * @param[in] f pointer to FILE structure
 * @return pointer to XFILE structure if OK, NULL if file stream is not known
 */
static XFILE * xio_getfptr( FILE * f )
{
    XFILE * xf = NULL;

    if ( xfiles_map.find(f) != xfiles_map.end() )
    {
        xf = xfiles_map[f];
    }

    return xf;
}

/*!
 * @brief Enable/Disable XIO debug output for EACH XIO function call
 * @param[in] verb 0 to disable (default), positive int to enable
 */
void xio_setverbose( int verb )
{
    xio_verbosity = verb;
}

/*!
 * @brief Get current XIO debug verbosity level
 * @return 0 as disabled, positive int as enabled
 */
int xio_getverbose( void )
{
    return xio_verbosity;
}

/*!
  * @brief fopen() with error handling and filename store
  * @param[in] fname filename to open
  * @param[in] mode  opening mode
  * @return same as ISO C fopen()
 */
FILE *xfopen( const char *fname, const char *mode )
{
    XFILE * xf;
    FILE *rc;

    F_ENTRY;

    ASSERT(!( (fname == NULL) || (mode == NULL) ),"NULL pointer as argument of function xfopen()\n");
    xprintf(MsgLog,"Opening file: '%s'\n", fname);
    rc = fopen( fname, mode );
    ASSERT(!( rc == NULL ),"Cannot open file %s with permissions %s\n", fname, mode );

    //store file name and file opening mode
    xf = (XFILE *)xmalloc(sizeof(XFILE));
    xf->filename = (char *)xmalloc(strlen(fname)+1);
    strcpy( xf->filename, fname );
    xf->mode = (char *)xmalloc(strlen(mode)+1);
    strcpy(xf->mode, mode);
    xf->lineno = 0;
    xfiles_map[rc] = xf;

    XIO_DEBUG( rc );

    return(rc);
}

/*!
 * @brief Flush file stream
 * @param[in,out] f pointer to FILE structure
 * @return same as ISO C fflush()
 */
int xfflush( FILE * f )
{
    XFILE * xf;

    ASSERT(!(f == NULL),"NULL as input argument\n");

    XIO_DEBUG( f );

    xf = xio_getfptr(f);
    if ( !xf )
    {
        XIO_WARN(f);
    }

    return fflush ( f );
}

/*!
 * @brief FCLOSE WITH ERROR HANDLING
 * @param[in,out] stream pointer to FILE structure
 * @return same as ISO C fclose()
 */
int xfclose( FILE *stream )
{
    XFILE * xf;
    int rc;

    F_ENTRY;

    ASSERT(!( stream == NULL ),"NULL pointer as argument of function xfclose()\n");

    XIO_DEBUG( stream );

    rc = fclose( stream );

    INPUT_CHECK(!( rc == EOF ),"Cannot close file %s\n", xio_getfname(stream) );

    if ( rc != EOF )
    {
        xf = xio_getfptr(stream);
        if ( xf )
        {
            xfiles_map.erase(stream);
            xfree( xf->filename );
            xfree( xf->mode );
            xfree( xf );
        }
        else
        {
            XIO_WARN(stream);
        }
    }

    return(rc);
}

/*!
 * @brief Reopen stream with different file or mode
 * @param[in]     filename  Name of the file to be opened
 * @param[in]     mode      File access mode
 * @param[in,out] stream    Pointer to a FILE object that identifies the stream to be reopened
 * @return                  Same as ISO C freopen()
 */
FILE * xfreopen( const char * filename, const char * mode, FILE * stream )
{
    XFILE * xf;
    FILE *rc;

    F_ENTRY;
    ASSERT(!( (mode == NULL) || (stream == NULL)),"Wrong arguments\n");

    rc = freopen( filename, mode, stream );

    ASSERT(!( rc == NULL ),"Cannot reopen file %s with permissions %s\n", filename, mode );

    //store file name and file opening mode
    xf = xio_getfptr(rc);
    if (xf)
    {
        //replace filename, if enough space
        if ( strlen(filename) > strlen(xf->filename) )
            xf->filename = (char *)xrealloc(xf->filename, strlen(filename)+1 );
        strcpy(xf->filename, filename);

        //replace filemode, if enough space
        if ( strlen(mode) > strlen(xf->mode) )
            xf->mode = (char *)xrealloc( xf->mode, strlen(mode)+1 );
        strcpy(xf->mode, mode);

        xf->lineno = 0;
    }
    else
    {
        //new file
        xf = (XFILE *)xmalloc(sizeof(XFILE));
        xf->filename = (char *)xmalloc(strlen(filename));
        strcpy(xf->filename, filename);
        xf->mode = (char *)xmalloc(strlen(mode));
        strcpy(xf->mode, mode);
        xf->lineno = 0;
        xfiles_map[rc] = xf;
    }

    XIO_DEBUG( rc );

    return(rc);
}

/*!
 * @brief FPRINTF WITH ERROR HANDLING
 */
int xfprintf( FILE *out, const char *fmt, ... )
{
    va_list argptr;
    int rc;

    F_ENTRY;

    ASSERT(!( (out == NULL) || (fmt == NULL) ),"NULL pointer as argument of function xfprintf()\n");
    va_start( argptr, fmt );
    rc = vfprintf( out, fmt, argptr );
    va_end( argptr );
    INPUT_CHECK(!( rc == EOF ),"Cannot write to file %s\n", xio_getfname(out) );

    return rc;
}

/*!
 * @brief FSCANF WITH ERROR HANDLING
 */
int xfscanf( FILE *in, const char *fmt, ... )
{
    va_list  argptr;
    int rc;

    F_ENTRY;

    ASSERT(!( (in == NULL) || (fmt == NULL) ),"NULL pointer as argument of function xfscanf()\n");
    va_start( argptr , fmt );
    rc = vfscanf( in, fmt, argptr );
    va_end( argptr );
    INPUT_CHECK(!( (rc == EOF) || (rc == 0) ),"Cannot read from file %s\n", xio_getfname(in) );

    return(rc);
}

/*!
 * @brief getc() with error handling and line count
 * @param[in,out] f pointer to FILE structure
 * @return same as ISO C getc()
 */
int xgetc( FILE * f )
{
    int rc;

    F_ENTRY;
    ASSERT(!(f == NULL), "NULL file\n");

    rc = xfgetc( f );
    XIO_DEBUG( f );

    return( rc );
}

/*!
 * @brief fgetc() with error handling and line count
 * @param[in,out] f pointer to FILE structure
 * @return same as ISO C fgetc()
 */
int xfgetc( FILE * f )
{
    int rc;
    XFILE * xf;

    F_ENTRY;
    ASSERT(!(f == NULL), "NULL file\n");

    rc = fgetc( f );

    //update line count
    xf = xio_getfptr( f );
    if ( xf )
    {
        if ( rc == '\n' )
            xf->lineno++;
    }
    else
    {
        XIO_WARN( f );
    }

    XIO_DEBUG( f );

    INPUT_CHECK(!( (rc == EOF) && (!feof(f))),"Cannot read from file '%s'\n", xio_getfname( f ) );

    return(rc);
}

/*!
 * @brief ungetc() with error handling and line count
 * @param[in]     c     character to push back
 * @param[in,out] f pointer to FILE structure with INPUT stream
 * @return same as ISO C ungetc()
 */
int xungetc( int c, FILE * f )
{
    int rc;
    XFILE * xf;

    F_ENTRY;
    ASSERT(!(f == NULL), "NULL file\n");

    rc = ungetc( c, f );

    //update line count
    xf = xio_getfptr( f );
    if ( xf )
    {
        if ( ( rc == '\n' ) && ( xf->lineno > 0 ) )
            xf->lineno--;
    }
    else
    {
        XIO_WARN( f );
    }

    XIO_DEBUG( f );

    INPUT_CHECK(!(rc == EOF), "Unable to push back character '%c' to file '%s'.\n", c, xio_getfname(f) );
    return( rc );
}

/*!
 * @brief Changes the name of the file or directory specified by oldname to newname
 * @param[in] oldname name of the file to be renamed and/or moved
 * @param[in] newname new name
 * @return same as ISO C rename()
 */
int xrename ( const char * oldname, const char * newname )
{
    int rc;

    F_ENTRY;
    ASSERT(!(( oldname == NULL) || (newname == NULL)), "NULL file name\n");

    rc = rename( oldname, newname );

    ASSERT(!(rc == 0),"Failed when renaming '%s' to '%s'\n", oldname, newname );
    return( rc );
}

/*!
 * @brief Read block of data from stream, handle errors
 * @param[out]    ptr    Pointer to a block of memory with a minimum size of (size*count) bytes.
 * @param[in]     size   Size in bytes of each element to be read.
 * @param[in]     count  Number of elements, each one with a size of size bytes.
 * @param[in,out] stream Pointer to a FILE object that specifies an input stream.
 * @return same as ISO C fread()
 */
size_t xfread( void * ptr, size_t size, size_t count, FILE * stream )
{
    size_t rc;

    F_ENTRY;
    ASSERT(!( (ptr == NULL) || ( stream == NULL) ),"Incorrect arguments\n");

    rc = fread( ptr, size, count, stream );

    XIO_DEBUG( stream );

    INPUT_CHECK(!( (rc != count) && (!feof(stream))),"Cannot read from file '%s'\n", xio_getfname(stream) );

    return(rc);
}

/*!
 * @brief Write block of data to stream, handle errors
 * @param[in]     ptr    Pointer to a block of memory with a minimum size of (size*count) bytes.
 * @param[in]     size   Size in bytes of each element to be written.
 * @param[in]     count  Number of elements, each one with a size of size bytes.
 * @param[in,out] stream Pointer to a FILE object that specifies an output stream.
 * @return same as ISO C fwrite()
 */
size_t xfwrite( const void * ptr, size_t size, size_t count, FILE * stream )
{
    size_t rc;

    F_ENTRY;
    ASSERT(!( (ptr == NULL) || (stream == NULL) ),"Incorrect arguments\n");

    rc = fwrite( ptr, size, count, stream );

    XIO_DEBUG( stream );

    INPUT_CHECK(!( (rc != count) ),"Cannot write to file '%s'\n", xio_getfname(stream) );

    return(rc);
}

/*!
 * @brief FGETS WITH ERROR HANDLING and line count
 * @param[out] s     Pointer to an array of chars where the string read is stored
 * @param[in] n      Maximum number of characters to be read (including the final null-character)
 * @param[in,out] in Pointer to a FILE
 * @return           same as ISO C fgets()
 */
char *xfgets( char *s, int n, FILE *in )
{
    XFILE * xf;
    char *rc = NULL;

    F_ENTRY;

    ASSERT(!( (s == NULL) || (in == NULL) ),"Incorrect arguments of function xfgets()\n");
    rc = fgets( s, n, in );

    //update line count
    xf = xio_getfptr(in);
    if ( xf )
    {
        //check if complete line is read
        if ( rc )
        {
            char c;
            c = s[strlen(s)-1];
            if ( c == '\n' )
            {
                xf->lineno++;
            }
        }
    }
    else
    {
        XIO_WARN(in);
    }

    XIO_DEBUG( in );

    INPUT_CHECK(!( (rc == NULL) && (!feof(in))),"Cannot read from file '%s'\n", xio_getfname(in) );

    return(rc);
}

/*!
 * @brief Rewind file, handle line count
 * @param[in,out] f pointer to FILE structure
 */
void xrewind( FILE * f )
{
    XFILE * xf;

    ASSERT(!(f == NULL),"NULL file argument in xrewind()\n");

    xf = xio_getfptr(f);
    if ( xf )
    {
        xf->lineno = 0;
    }
    else
    {
        XIO_WARN(f);
    }

    rewind( f );

    XIO_DEBUG( f );
}

/*!
 * @brief Check END OF FILE
 * @param[in] f pointer to FILE structure
 * @return same as ISO C feof()
 */
int xfeof ( FILE * f )
{
    XFILE * xf;
    int rc;

    xf = xio_getfptr(f);
    if ( !xf )
    {
        XIO_WARN(f);
    }

    rc = feof ( f );

    XIO_DEBUG( f );

    return rc;
}
