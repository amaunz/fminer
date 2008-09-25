/* src/config.h.  Generated from config.h.in by configure.  */
/* src/config.h.in.  Generated from configure.in by autoheader.  */

/* Where the data files are located */
#define BABEL_DATADIR "/usr/local/share/openbabel"

/* The version of Open Babel */
#define BABEL_VERSION "2.1.1"

/* Used to export symbols for DLL / shared library builds */
#if defined(WIN32)
 #if defined(USING_OBDLL) // e.g. in src/main.cpp
  #define EXTERN __declspec(dllimport) extern
 #else
  #define EXTERN __declspec(dllexport) extern
 #endif
#else //Everything else (behaviour as original)
 #define EXTERN extern 
#endif


/* Define to 1 if the system has the type `clock_t'. */
#define HAVE_CLOCK_T 1

/* Define to 1 if you have the <conio.h> header file. */
/* #undef HAVE_CONIO_H */

/* Define to 1 if you have the <ctype.h> header file. */
#define HAVE_CTYPE_H 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the <fstream> header file. */
#define HAVE_FSTREAM 1

/* Define to 1 if you have the <fstream.h> header file. */
#define HAVE_FSTREAM_H 1

/* Define to 1 if you have the <iconv.h> header file. */
#define HAVE_ICONV_H 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the <iostream> header file. */
#define HAVE_IOSTREAM 1

/* Define to 1 if you have the <iostream.h> header file. */
#define HAVE_IOSTREAM_H 1

/* Define to 1 if you have the `m' library (-lm). */
#define HAVE_LIBM 1

/* Define to 1 if you have the `z' library (-lz). */
#define HAVE_LIBZ 1

/* Define to 1 if you have the <math.h> header file. */
#define HAVE_MATH_H 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the `rint' function. */
#define HAVE_RINT 1

/* Define to 1 if you have the `snprintf' function. */
#define HAVE_SNPRINTF 1

/* Define to 1 if you have the `sranddev' function. */
/* #undef HAVE_SRANDDEV */

/* Define to 1 if you have the <sstream> header file. */
#define HAVE_SSTREAM 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdio.h> header file. */
#define HAVE_STDIO_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the `strncasecmp' function. */
#define HAVE_STRNCASECMP 1

/* Define to 1 if you have the <strstream> header file. */
#define HAVE_STRSTREAM 1

/* Define to 1 if you have the <strstream.h> header file. */
/* #undef HAVE_STRSTREAM_H */

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/time.h> header file. */
#define HAVE_SYS_TIME_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <time.h> header file. */
#define HAVE_TIME_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if you have the <zlib.h> header file. */
#define HAVE_ZLIB_H 1

/* The file extension used for shared modules */
#define MODULE_EXTENSION ".so"

/* Used to export symbols for DLL / shared library builds */
#if defined(WIN32)
 #if defined(USING_OBDLL) // e.g. in src/main.cpp
  #define OBAPI __declspec(dllimport)
 #else
  #define OBAPI __declspec(dllexport)
 #endif
#else //Everything else (behaviour as original)
 #define OBAPI 
#endif


/* Used to export symbols for DLL / shared library builds */
#if defined(WIN32)
 #if defined(USING_OBDLL) // e.g. in src/main.cpp
  #define OBCOMMON __declspec(dllimport)
 #else
  #define OBCOMMON __declspec(dllexport)
 #endif
#else //Everything else (behaviour as original)
 #define OBCOMMON 
#endif


/* Used to export symbols for DLL / shared library builds */
#if defined(WIN32)
 #if defined(USING_OBDLL) // e.g. in src/main.cpp
  #define OBCONV __declspec(dllimport)
 #else
  #define OBCONV __declspec(dllexport)
 #endif
#else //Everything else (behaviour as original)
 #define OBCONV
#endif


/* Used to export symbols for DLL / shared library builds */
#if defined(WIN32)
 #if defined(USING_OBDLL) // e.g. in src/main.cpp
  #define OBERROR __declspec(dllimport)
 #else
  #define OBERROR __declspec(dllexport)
 #endif
#else //Everything else (behaviour as original)
 #define OBERROR 
#endif


/* Used to export symbols for DLL / shared library builds */
#if defined(WIN32)
 #if defined(USING_OBDLL) // e.g. in src/main.cpp
  #define OBFPTR __declspec(dllimport)
 #else
  #define OBFPTR __declspec(dllexport)
 #endif
#else //Everything else (behaviour as original)
 #define OBFPTR 
#endif


/* Define to the address where bug reports for this package should be sent. */

/* Define to the full name of this package. */

/* Define to the full name and version of this package. */

/* Define to the one symbol short name of this package. */

/* Define to the version of this package. */

/* set if scandir needs a const */
#define SCANDIR_CONST const

/* set if scandir needs a const */
#define SCANDIR_T (int (*)(const dirent *))


#if !HAVE_SNPRINTF
extern "C" int snprintf( char *, size_t, const char *, /* args */ ...);
#endif


/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Define to 1 if you can safely include both <sys/time.h> and <time.h>. */
#define TIME_WITH_SYS_TIME 1

/* Define to 1 if your processor stores words with the most significant byte
   first (like Motorola and SPARC, unlike Intel and VAX). */
/* #undef WORDS_BIGENDIAN */

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif
