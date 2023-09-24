#include <stdlib.h>
#include <string.h>
#include "io.h"
#include "memory.h"

int status_flag = 1;    /* By default, print all status messages */
int warnings_flag = 1;  /* By default, print all warnings */
int errors_flag = 1;    /* By default, print all errors */
#define EXIT_NOW(val) { SetMessageLog(NULL); free_all_memory(); exit(val); }

static FILE *message_log = stderr; /* File to write all messages */

typedef struct file_table {
   FILE *file;       /* Handle for this file */
   char *name;       /* Name of this file */
   int   line;       /* Current line in this file */
   } file_table;

static char *CurrentFileName = NULL;
static file_table File_List[MAX_FILE_DEPTH];
static int File_Name_Depth = 0;

/* To work with LEX & YACC, we need to be able to manipulate "yyin" and
   "yylineno" */
extern FILE *yyin;
extern int yylineno;

FILE *
PathFileOpen(char *environ, char *str, char *options)
{
   char *ename, *path, tenviron[128], filename[256];
   FILE *file;

   /* If we don't have an environment string, then
      just use the file name as given. */
   if (environ == NULL)
      return fopen(str, options);
   else {
      /* If it's in the local directory, then we use that */
      if ((file = fopen(str, options)) != NULL)
         return file;

      /* Check all places in the environment variable */
      ename = getenv(environ);

      /* Make local to prevent overwriting something important */
      sprintf(tenviron, "%s", ename);

      /* Step through all the directories listed in the environment
         variable. */
      for (path = strtok(ename, " ;");
           path != NULL;
           path = strtok(NULL, " ;")) {
         sprintf(filename, "%s/%s", path, str);
         if ((file = fopen(filename, options)) != NULL)
            return file;
         }

      /* Desperate last try in the current directory */
      return fopen(str, options);
      }
}

void
SetMessageLog(char *str)
{
   if (message_log != NULL &&
       message_log != stderr)
      fclose(message_log);
   if (str == NULL)
      message_log = stderr;
   else {
      if ((message_log = fopen(str, "w")) == NULL)
         error("Cannot open %s\n", str);
      }
}

void
SetInputFile(char *str)
{
   yylineno = 1;
   if (str == NULL) {
      yyin = stdin;
      CurrentFileName = "<stdin>";
      }
   else {
      CurrentFileName = str;
      if ((yyin = fopen(str, "r")) == NULL)
         error("Can't open file '%s'\n", str);
      }
   File_List[0].file = yyin;
   File_List[0].name = str;
   File_List[0].line = 0;
   File_Name_Depth = 1;
}

int
yywrap()
{
   if (File_Name_Depth > 0) {
      File_Name_Depth--;
      if (File_List[File_Name_Depth].file != stdin &&
          fclose(File_List[File_Name_Depth].file) != 0)
         error("Failed to close file '%s'\n",
               File_List[File_Name_Depth].name);
      if (File_Name_Depth == 0)
         return 1;
      free(File_List[File_Name_Depth].name);
      yyin            = File_List[File_Name_Depth-1].file;
      CurrentFileName = File_List[File_Name_Depth-1].name;
      yylineno        = File_List[File_Name_Depth].line;
      return 0;
      }
   else {
      yyin = stdin;
      return 1;
   }
}

void
include_file_action(char *environ, char *name)
{
   int len;
   char *strtmp;

   if (File_Name_Depth >= MAX_FILE_DEPTH)
      warning("Include files nested too deep\n");

   /* Now that we have the name of the include file, we
      stash it away for when the EOF is reached. */
   File_List[File_Name_Depth].file = PathFileOpen(environ, name, "rt");
   if (File_List[File_Name_Depth].file == NULL)
      warning("Failed to open include file: '%s'\n", name);
   else {
      len = strlen(name);
      strtmp = malloc((len + 1) * sizeof(char));
      if (strtmp == NULL)
         error("Failed to allocate space for a string constant\n");
      memcpy(strtmp, name, len);
      strtmp[len] = '\0';
      File_List[File_Name_Depth].line = yylineno;
      File_List[File_Name_Depth].name = strtmp;
      yyin = File_List[File_Name_Depth].file;
      CurrentFileName = strtmp;
      yylineno = 1;
      File_Name_Depth++;
      }
}

void
message(char *fmt, ...)
{
   va_list ap;
   if (status_flag) {
      va_start(ap, fmt);
      vfprintf(message_log, fmt, ap);
      va_end(ap);
      }
}

void
status(char *fmt, ...)
{
   va_list ap;
   if (status_flag) {
      va_start(ap, fmt);
      vfprintf(stderr, fmt, ap);
      va_end(ap);
      }
}

void
warning(char *fmt, ...)
{
   va_list ap;

   if (warnings_flag) {
      fprintf(message_log, "WARNING: ");
      va_start(ap, fmt);
      vfprintf(message_log, fmt, ap);
      va_end(ap);
      if (File_Name_Depth)
         fprintf(message_log, "On or near line %d of file %s\n",
                 yylineno, CurrentFileName);
      }
}

void
error(char *fmt, ...)
{
   va_list ap;

   if (errors_flag) {
      fprintf(message_log, "ERROR: ");
      va_start(ap, fmt);
      vfprintf(message_log, fmt, ap);
      va_end(ap);
      if (File_Name_Depth)
         fprintf(message_log, " on or near line %d of file %s\n",
                 yylineno, CurrentFileName);
      }
   EXIT_NOW(1);
}

void
fatal(char *fmt, ...)
{
   va_list ap;
   if (errors_flag) {
      fprintf(message_log, "FATAL: ");
      va_start(ap, fmt);
      vfprintf(message_log, fmt, ap);
      va_end(ap);
      if (File_Name_Depth)
         fprintf(message_log, "On or near line %d of file %s\n",
                 yylineno, CurrentFileName);
      }
   EXIT_NOW(2);
}
