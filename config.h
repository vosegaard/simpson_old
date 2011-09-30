/* #undef ENABLE_MINUIT */
#define SECURE_SERVER 1
/* #undef f2cFortran */
#define HAVE_GETLINE 1
#define REVERSEBYTES 1
/* #undef PIVOT_PHASE */
#define HAVE_F_EXIT 1
/* #undef ENABLE_SIMMOL */
/* #undef ENABLE_OPENGL */
#define HAVE_POW_DI 1

/* Disable networking */
/* #undef DISABLE_NETWORK */

/* Compile with OpenGL support */
/* #undef ENABLE_OPENGL */

/* Compile with Tk support */
#define ENABLE_TK 0

//#define MPI
#ifdef MPI //include header file for MPI setup
#include <mpi.h>
#endif

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you don't have `vprintf' but do have `_doprnt.' */
/* #undef HAVE_DOPRNT */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the `f2c' library (-lf2c). */
/* #undef HAVE_LIBF2C */

/* Define to 1 if you have the `for' library (-lfor). */
/* #undef HAVE_LIBFOR */

/* Define to 1 if you have the `ftn' library (-lftn). */
/* #undef HAVE_LIBFTN */

/* Define to 1 if you have the `Futil' library (-lFutil). */
/* #undef HAVE_LIBFUTIL */

/* Define to 1 if you have the `g2c' library (-lg2c). */
#define HAVE_LIBG2C 1

/* Define to 1 if you have the `m' library (-lm). */
#define HAVE_LIBM 1

/* Define to 1 if you have the `m_4sqrt' library (-lm_4sqrt). */
/* #undef HAVE_LIBM_4SQRT */

/* Define to 1 if you have the `ots' library (-lots). */
/* #undef HAVE_LIBOTS */

/* Define to 1 if you have the `Ufor' library (-lUfor). */
/* #undef HAVE_LIBUFOR */

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the `strstr' function. */
#define HAVE_STRSTR 1

/* Define to 1 if you have the `strtod' function. */
#define HAVE_STRTOD 1

/* Define to 1 if you have the `strtol' function. */
#define HAVE_STRTOL 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if you have the `vprintf' function. */
#define HAVE_VPRINTF 1

/* Name of package */
#define PACKAGE "simpson"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT ""

/* Define to the full name of this package. */
#define PACKAGE_NAME "configure.in"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "configure.in minuit/minuit/minuit/d506cm.inc"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "configure-in"

/* Define to the version of this package. */
#define PACKAGE_VERSION "minuit/minuit/minuit/d506cm.inc"

/* Make the code multithreaded */
/* #undef PARALLEL */

/* Define as the return type of signal handlers (`int' or `void'). */
#define RETSIGTYPE void

/* Enable serverlogging */
/* #undef SERVER_LOGGING */

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1


/* Compile GD */
#define WITH_GD 1

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

/* Define to `unsigned' if <sys/types.h> does not define. */
/* #undef size_t */
