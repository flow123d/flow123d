#define MAXBUFF 1024
#define _VERSION_ "0.1.1"

// GCC accept only __func__
#define __FUNC__        __func__

#ifdef HAVE_WINDOWS
  #include <dir.h>
#else
// for a linux system we assume glic lib
// // with support of ISOC99 functions
// // #define _ISOC99_SOURCE
// //     #define _BSD_SOURCE
// //       #include <stdio.h>
#include <unistd.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

//   // MAXPATH depends on filesystem, if following makes problems
//   // "pathconf" approach form glibc should be used
#define MAXPATH 200
//     //             #include <sys/stat.h>
//     //             #include <errno.h>
//     //               #define DLL_EXPORT
//#define strcmpi strcasecmp
//       //                   #define PATH_SEP "/"
//       //
#endif 
