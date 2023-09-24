#if !defined(__POLYRAY_IO_DEFS)
#define __POLYRAY_IO_DEFS

#include <stdio.h>
#include <stdarg.h>

#define MAX_FILE_DEPTH       (8)

/* Central file open routine - this will automatically search
   the directories in the environment variable for places that the
   file might reside */
extern FILE *PathFileOpen(char *environ, char *str, char *options);

#if defined( VISUALC )
#define RTCDECL __cdecl
#else
#define RTCDECL
#endif

/* Message and Error reporting */
void RTCDECL message(char *fmt, ...);  /* Simply display the message      */
void RTCDECL status(char *fmt, ...);   /* Display message to stderr       */
void RTCDECL warning(char *fmt, ...);  /* Display "WARNING: ..."          */
void RTCDECL error(char *fmt, ...);    /* Display "ERROR: ...", then exit */
void RTCDECL fatal(char *fmt, ...);    /* Display "FATAL: ...", then exit */

extern void SetInputFile(char *);
extern void SetMessageLog(char *);
extern void include_file_action(char *environ, char *name);
extern int yywrap(void);
extern int status_flag;   /* Set to 0 to eliminate status output */
extern int warnings_flag; /* Set to 0 to eliminate warning output */
extern int errors_flag;   /* Set to 0 to eliminate error output */

#endif /* __POLYRAY_IO_DEFS */

