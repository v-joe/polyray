#if !defined(__POLYRAY_PARSE_DEFS)
#define __POLYRAY_PARSE_DEFS

/* Parser support functions and variables */
extern void ReadSceneFile(char *str);
extern int yylex(void);
extern int yylook(void);
extern int yyinput(void);
extern void yyoutput(int);
extern void yyunput(int);
extern int yyparse(void);
extern char *yyptok(int);
extern Vec White;

#endif /* __POLYRAY_PARSE_DEFS */

