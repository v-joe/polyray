#include <stdio.h>
#define U(x) x
#define NLSTATE yyprevious=YYNEWLINE
#define BEGIN yybgin = yysvec + 1 +
#define INITIAL 0
#define YYLERR yysvec
#define YYSTATE (yysetate-yysvec-1)
#define YYOPTIM 1
#define YYLMAX 200
#define output(c) putc(c, yyout)
#define input() (((yytchar=yysptr>yysbuf?U(*--yysptr):getc(yyin))==10?(yylineno++,yytchar):yytchar)==EOF?0:yytchar)
#define unput(c) {yytchar=(c);if (yytchar=='\n') yylineno--; *yysptr++=yytchar;}
#define yymore() (yymorfg=1)
#define ECHO fprintf(yyout, "%s",yytext)
#define REJECT { nstr = yyreject(); goto yyfussy;}
int yyleng; extern char yytext[];
int yymorfg;
extern char *yysptr, yysbuf[];
int yytchar;
FILE *yyin = { stdin }, *yyout = { stdout };
extern int yylineno;
struct yysvf { 
	struct yywork *yystoff;
	struct yysvf *yyother;
	int *yystops;};
struct yysvf *yyestate;
extern struct yysvf yysvec[], *yybgin;
extern int yywrap();
extern int yylook();
extern int yyback(int *p, int m);
/*

  Copyright (C) 1993, Alexander Enzmann, All rights reserved.
 
  You may not distribute this software, in whole or in part,
  without the express consent of the authors.
 
  There is no warranty or other guarantee of fitness of this software
  for any purpose.  It is provided solely "as is".

*/
#include "defs.h"
#include "ytab.h"
#include "io.h"
#include "memory.h"
#include "symtab.h"
#include "parse.h"

static void skipcom(void);

#define END(v) (v - 1 + sizeof(v) / sizeof(v[0]))
static int lookup(void);
static int copy_yystring(void);
# define YYNEWLINE 10
yylex() {
int nstr; extern int yyprevious;
while ((nstr = yylook()) >= 0)
yyfussy: switch(nstr) {
case 0:
if (yywrap()) return (0); break;
case 1:                        ;
break;
case 2:                           skipcom();
break;
case 3:                        ;
break;
case 4:                             ;
break;
case 5:                   { return copy_yystring(); }
break;
case 6: { yylval.flt = atof(yytext); return NUM; }
break;
case 7: { yylval.flt = atof(yytext); return NUM; }
break;
case 8:                  return lookup();
break;
case 9:                           return AND_SYM;
break;
case 10:                           return OR_SYM;
break;
case 11:                           return LTEQ_SYM;
break;
case 12:                           return EQUAL_SYM;
break;
case 13:                           return GTEQ_SYM;
break;
case 14:                           return NEQUAL_SYM;
break;
case 15:                              return yytext[0];
break;
case -1:
break;
default:
fprintf(yyout,"bad switch yylook %d",nstr);
} return(0); }
/* End of yylex */

/* Screen for reserved words */
struct lookup_table_blk {
   char *name;
   int val;
   } lookup_table[] = {
      {"Blinn",                  BLINN},
      {"Cook",                   COOK},
      {"Ferrari",                FERRARI},
      {"Gaussian",               GAUSSIAN},
      {"Phong",                  PHONG},
      {"Reitz",                  REITZ},
      {"Sturm",                  STURM},
      {"Vieta",                  VIETA},
      {"acceleration",           ACCELERATION},
      {"ambient",                AMBIENT},
      {"angle",                  ANGLE},
      {"antialias",              ANTIALIAS},
      {"antialias_threshold",    ANTIALIAS_THRESHOLD},
      {"aperture",               APERTURE},
      {"aspect",                 ASPECT},
      {"at",                     AT},
      {"avoid",                  AVOID},
      {"background",             BACKGROUND},
      {"bezier",                 BEZIER},
      {"birth",                  BIRTH},
      {"blob",                   BLOB},
      {"bounding_box",           BOUNDING_BOX},
      {"box",                    BOX},
      {"brilliance",             BRILLIANCE},
      {"bump_scale",             BUMP_SCALE},
      {"checker",                CHECKER},
      {"cheight_field",          CHEIGHT_FIELD},
      {"cheight_fn",             CHEIGHT_FN},
      {"color",                  COLOR},
      {"color_map",              COLOR_MAP},
      {"concave",                CONCAVE},
      {"cone",                   CONE},
      {"contour",                CONTOUR},
      {"count",                  COUNT},
      {"cylinder",               CYLINDER},
      {"death",                  DEATH},
      {"define",                 DEFINE},
      {"depth",                  DEPTH},
      {"depthmapped_light",      DEPTHMAPPED_LIGHT},
      {"diffuse",                DIFFUSE},
      {"directional_light",      DIRECTIONAL_LIGHT},
      {"disc",                   DISC},
      {"displace",               DISPLACE},
      {"dither",                 DITHER},
      {"draw",                   DRAW},
      {"else",                   ELSE},
      {"end_frame",              END_FRAME},
      {"file_flush",             FILE_FLUSH},
      {"flare",                  FLARE},
      {"flock",                  FLOCK},
      {"focal_distance",         FOCAL_DISTANCE},
      {"frame",                  FRAME},
      {"frame_time",             FRAME_TIME},
      {"frequency",              FREQUENCY},
      {"from",                   FROM},
      {"function",               FUNCTION},
      {"glyph",                  GLYPH},
      {"gridded",                GRIDDED},
      {"haze",                   HAZE},
      {"height_field",           HEIGHT_FIELD},
      {"height_fn",              HEIGHT_FN},
      {"hexagon",                HEXAGON},
      {"hither",                 HITHER},
      {"hypertexture",           HYPERTEXTURE},
      {"if",                     IF},
      {"image_format",           IMAGE_FORMAT},
      {"image_window",           IMAGE_WINDOW},
      {"include",                INCLUDE},
      {"indexed",                INDEXED},
      {"lathe",                  LATHE},
      {"layered",                LAYERED},
      {"light",                  LIGHT},
      {"lookup_fn",              LOOKUP_FUNCTION},
      {"max_samples",            MAX_SAMPLES},
      {"max_trace_depth",        MAX_TRACE_DEPTH},
      {"microfacet",             MICROFACET},
      {"no_shadow",              NO_SHADOW},
      {"noeval",                 NOEVAL},
      {"noise",                  NOISE},
      {"normal",                 NORMAL},
      {"nurb",                   NURB},
      {"object",                 OBJECT},
      {"octaves",                OCTAVES},
      {"outfile",                OUTFILE},
      {"parabola",               PARABOLA},
      {"parametric",             PARAMETRIC},
      {"particle",               PARTICLE},
      {"patch",                  PATCH},
      {"phase",                  PHASE},
      {"pixel_encoding",         PIXEL_ENCODING},
      {"pixelsize",              PIXELSIZE},
      {"plane",                  PLANE},
      {"point",                  POINT},
      {"polygon",                POLYGON},
      {"polynomial",             POLYNOMIAL},
      {"position",               POSITION},
      {"position_fn",            POSITION_FUNCTION},
      {"position_scale",         POSITION_SCALE},
      {"raw",                    RAW},
      {"reflection",             REFLECTION},
      {"resolution",             RESOLUTION},
      {"root_solver",            ROOT_SOLVER},
      {"rotate",                 ROTATE},
      {"scale",                  SCALE},
      {"seed",                   SEED},
      {"shading_flags",          SHADING_FLAGS},
      {"shear",                  SHEAR},
      {"sheight_field",          SHEIGHT_FIELD},
      {"sheight_fn",             SHEIGHT_FN},
      {"size",                   SIZE},
      {"smooth_cheight_field",   SMOOTH_CHEIGHT_FIELD},
      {"smooth_cheight_fn",      SMOOTH_CHEIGHT_FN},
      {"smooth_height_field",    SMOOTH_HEIGHT_FIELD},
      {"smooth_height_fn",       SMOOTH_HEIGHT_FN},
      {"smooth_sheight_field",   SMOOTH_SHEIGHT_FIELD},
      {"smooth_sheight_fn",      SMOOTH_SHEIGHT_FN},
      {"spacing",                SPACING},
      {"special",                SPECIAL},
      {"specular",               SPECULAR},
      {"sphere",                 SPHERE},
      {"spot_light",             SPOT_LIGHT},
      {"start_frame",            START_FRAME},
      {"static",                 STATIC},
      {"summed",                 SUMMED},
      {"superq",                 SUPERQ},
      {"surface",                SURFACE},
      {"sweep",                  SWEEP},
      {"system",                 SYSTEM},
      {"texture",                TEXTURE},
      {"texture_map",            TEXTURE_MAP},
      {"textured_light",         TEXTURED_LIGHT},
      {"torus",                  TORUS},
      {"total_frames",           TOTAL_FRAMES},
      {"transform",              TRANSFORM},
      {"translate",              TRANSLATE},
      {"transmission",           TRANSMISSION},
      {"turbulence",             TURBULENCE},
      {"u_steps",                U_STEPS},
      {"up",                     UP},
      {"uv",                     UV},
      {"uv_bounds",              UV_BOUNDS},
      {"uv_steps",               UV_STEPS},
      {"v_steps",                V_STEPS},
      {"velocity",               VELOCITY},
      {"viewpoint",              VIEWPOINT},
      {"w_steps",                W_STEPS},
      {"yon",                    YON}
      };
typedef struct lookup_table_blk *tabptr;

int
lookup(void)
{
   tabptr low  = lookup_table;
   tabptr high = END(lookup_table);
   tabptr mid;
   char *strtmp;
   int c, l1;
   void *d;
#ifdef DEBUG_FN_CALLS
message("lookup\n");
#endif

   while (low <= high) {
      mid = low + (high - low) / 2;
      if ((c = strcmp(mid->name, yytext)) == 0) {
         return mid->val;
         }
      else if (c < 0)
         low = mid + 1;
      else
         high = mid - 1;
      }

   /* If the string isn't predefined, then look for it in the symbol table */
   l1 = strlen(yytext);
   strtmp = Lookup_String(yytext);
   if (strtmp == NULL) {
      strtmp = polyray_malloc((l1+1) * sizeof(char));
      if (strtmp == NULL)
         error("Failed to allocate space for a string constant\n");
      memcpy(strtmp, yytext, l1+1);
      yylval.name = strtmp;
      return TOKEN;
      }
   else {
      Lookup_Definition(yytext, &c, &d);
      switch (c) {
      case T_OBJECT:
         yylval.name = strtmp;
         return OBJECT_SYM;
      case T_SURFACE:
         yylval.name = strtmp;
         return SURFACE_SYM;
      case T_TEXTURE:
         yylval.name = strtmp;
         return TEXTURE_SYM;
      case T_TEXTURE_MAP:
         yylval.name = strtmp;
         return TEXTURE_MAP_SYM;
      case T_TRANSFORM:
         yylval.name = strtmp;
         return TRANSFORM_SYM;
      case T_EXPRESSION:
         yylval.name = strtmp;
         return EXPRESSION_SYM;
      case T_PARTICLE:
         yylval.name = strtmp;
         return PARTICLE_SYM;
      default:
         error("Bad type in symbol table: '%s' - %d\n", yytext, c);
	 return 0; /* To remove warning about non-return value */
      } }
}

int
copy_yystring(void)
{
   int len = strlen(yytext);
   char *strtmp;
#ifdef DEBUG_FN_CALLS
message("copy_yystring\n");
#endif
   strtmp = polyray_malloc((len - 1) * sizeof(char));
   if (strtmp == NULL)
      error("Failed to allocate space for a string constant\n");
   memcpy(strtmp, &yytext[1], len-2);
   strtmp[len-2] = '\0';
   yylval.name = strtmp;
   return STRING;
}

void
ReadSceneFile(char *str)
{
#ifdef DEBUG_FN_CALLS
message("ReadSceneFile\n");
#endif
   SetInputFile(str);
   if (yyparse() == 1)
      error("Invalid input specification\n");
   if (yyin != stdin)
      fclose(yyin);
}

/* Skip over comments */
static void
skipcom()
{
   int c;

   c = input();
   while (c != 0) {
      while (c == '*')
         if ((c = input()) == '/')
            return;
      c = input();
      }
   error("EOF inside comment");
}
int yyvstop[] = {
0, 

7, 
0,

7, 
0,

15, 
0,

1, 
15, 
0,

4, 
0,

15, 
0,

15, 
0,

15, 
0,

15, 
0,

15, 
0,

6, 
7, 
8, 
15, 
0,

15, 
0,

15, 
0,

15, 
0,

8, 
15, 
0,

15, 
0,

14, 
0,

5, 
0,

9, 
0,

7, 
0,

2, 
0,

-3,
0,

6, 
0,

6, 
7, 
8, 
0,

8, 
0,

8, 
0,

11, 
0,

12, 
0,

13, 
0,

10, 
0,

3, 
0,

6, 
7, 
0,

6, 
7, 
8, 
0,

7, 
0,

6, 
0,

6, 
7, 
0,
0};
#define YYTYPE char
struct yywork { YYTYPE verify, advance; } yycrank[] ={
0,0,	0,0,	1,3,	0,0,	
0,0,	0,0,	0,0,	0,0,	
7,18,	0,0,	1,4,	1,5,	
18,0,	0,0,	0,0,	0,0,	
7,18,	7,0,	0,0,	0,0,	
0,0,	0,0,	0,0,	0,0,	
0,0,	0,0,	0,0,	0,0,	
0,0,	0,0,	0,0,	0,0,	
0,0,	0,0,	1,6,	1,7,	
0,0,	0,0,	0,0,	1,8,	
2,6,	7,19,	8,20,	0,0,	
1,3,	2,8,	10,22,	1,9,	
1,10,	1,11,	7,18,	10,23,	
0,0,	7,18,	2,10,	7,18,	
40,36,	0,0,	40,36,	0,0,	
0,0,	1,12,	1,13,	1,14,	
6,17,	12,28,	1,15,	2,12,	
2,13,	2,14,	13,29,	14,30,	
7,18,	9,21,	9,21,	9,21,	
9,21,	9,21,	9,21,	9,21,	
9,21,	9,21,	9,21,	27,36,	
0,0,	27,36,	0,0,	0,0,	
27,37,	27,37,	27,37,	27,37,	
27,37,	27,37,	27,37,	27,37,	
27,37,	27,37,	32,38,	0,0,	
32,38,	0,0,	21,32,	32,39,	
32,39,	32,39,	32,39,	32,39,	
32,39,	32,39,	32,39,	32,39,	
32,39,	11,24,	34,40,	11,25,	
11,25,	11,25,	11,25,	11,25,	
11,25,	11,25,	11,25,	11,25,	
11,25,	1,16,	16,31,	0,0,	
0,0,	0,0,	0,0,	2,16,	
11,26,	11,26,	11,26,	11,26,	
11,26,	11,26,	11,26,	11,26,	
11,26,	11,26,	11,26,	11,26,	
11,26,	11,26,	11,26,	11,26,	
11,26,	11,26,	11,26,	11,26,	
11,26,	11,26,	11,26,	11,26,	
11,26,	11,26,	0,0,	0,0,	
0,0,	0,0,	11,26,	0,0,	
11,26,	11,26,	11,26,	11,26,	
11,27,	11,26,	11,26,	11,26,	
11,26,	11,26,	11,26,	11,26,	
11,26,	11,26,	11,26,	11,26,	
11,26,	11,26,	11,26,	11,26,	
11,26,	11,26,	11,26,	11,26,	
11,26,	11,26,	15,26,	15,26,	
15,26,	15,26,	15,26,	15,26,	
15,26,	15,26,	15,26,	15,26,	
0,0,	0,0,	0,0,	0,0,	
0,0,	0,0,	0,0,	15,26,	
15,26,	15,26,	15,26,	15,26,	
15,26,	15,26,	15,26,	15,26,	
15,26,	15,26,	15,26,	15,26,	
15,26,	15,26,	15,26,	15,26,	
15,26,	15,26,	15,26,	15,26,	
15,26,	15,26,	15,26,	15,26,	
15,26,	0,0,	0,0,	0,0,	
0,0,	15,26,	0,0,	15,26,	
15,26,	15,26,	15,26,	15,26,	
15,26,	15,26,	15,26,	15,26,	
15,26,	15,26,	15,26,	15,26,	
15,26,	15,26,	15,26,	15,26,	
15,26,	15,26,	15,26,	15,26,	
15,26,	15,26,	15,26,	15,26,	
15,26,	23,23,	0,0,	0,0,	
0,0,	0,0,	0,0,	0,0,	
0,0,	23,23,	23,33,	24,34,	
24,34,	24,34,	24,34,	24,34,	
24,34,	24,34,	24,34,	24,34,	
24,34,	36,43,	36,43,	36,43,	
36,43,	36,43,	36,43,	36,43,	
36,43,	36,43,	36,43,	0,0,	
0,0,	0,0,	23,23,	0,0,	
0,0,	0,0,	0,0,	0,0,	
0,0,	0,0,	0,0,	23,23,	
0,0,	35,41,	23,23,	35,41,	
23,23,	0,0,	35,42,	35,42,	
35,42,	35,42,	35,42,	35,42,	
35,42,	35,42,	35,42,	35,42,	
0,0,	0,0,	0,0,	0,0,	
24,35,	23,23,	37,37,	37,37,	
37,37,	37,37,	37,37,	37,37,	
37,37,	37,37,	37,37,	37,37,	
38,39,	38,39,	38,39,	38,39,	
38,39,	38,39,	38,39,	38,39,	
38,39,	38,39,	41,42,	41,42,	
41,42,	41,42,	41,42,	41,42,	
41,42,	41,42,	41,42,	41,42,	
0,0};
struct yysvf yysvec[] = {
0,	0,	0,
yycrank+-1,	0,		yyvstop+1,		/* state 0 */
yycrank+-7,	yysvec+1,	yyvstop+3,		/* state 1 */
yycrank+0,	0,		yyvstop+5,		/* state 2 */
yycrank+0,	0,		yyvstop+7,		/* state 3 */
yycrank+0,	0,		yyvstop+10,		/* state 4 */
yycrank+3,	0,		yyvstop+12,		/* state 5 */
yycrank+-7,	0,		yyvstop+14,		/* state 6 */
yycrank+4,	0,		yyvstop+16,		/* state 7 */
yycrank+25,	0,		yyvstop+18,		/* state 8 */
yycrank+4,	0,		yyvstop+20,		/* state 9 */
yycrank+67,	0,		yyvstop+22,		/* state 10 */
yycrank+4,	0,		yyvstop+27,		/* state 11 */
yycrank+9,	0,		yyvstop+29,		/* state 12 */
yycrank+10,	0,		yyvstop+31,		/* state 13 */
yycrank+142,	0,		yyvstop+33,		/* state 14 */
yycrank+2,	0,		yyvstop+36,		/* state 15 */
yycrank+0,	0,		yyvstop+38,		/* state 16 */
yycrank+-2,	yysvec+7,	0,			/* state 17 */
yycrank+0,	0,		yyvstop+40,		/* state 18 */
yycrank+0,	0,		yyvstop+42,		/* state 19 */
yycrank+1,	yysvec+9,	yyvstop+44,		/* state 20 */
yycrank+0,	0,		yyvstop+46,		/* state 21 */
yycrank+-264,	0,		yyvstop+48,		/* state 22 */
yycrank+227,	0,		yyvstop+50,		/* state 23 */
yycrank+0,	yysvec+11,	yyvstop+52,		/* state 24 */
yycrank+0,	yysvec+15,	yyvstop+56,		/* state 25 */
yycrank+40,	yysvec+15,	yyvstop+58,		/* state 26 */
yycrank+0,	0,		yyvstop+60,		/* state 27 */
yycrank+0,	0,		yyvstop+62,		/* state 28 */
yycrank+0,	0,		yyvstop+64,		/* state 29 */
yycrank+0,	0,		yyvstop+66,		/* state 30 */
yycrank+55,	0,		0,			/* state 31 */
yycrank+0,	0,		yyvstop+68,		/* state 32 */
yycrank+13,	yysvec+24,	yyvstop+70,		/* state 33 */
yycrank+266,	0,		0,			/* state 34 */
yycrank+237,	0,		0,			/* state 35 */
yycrank+282,	yysvec+15,	yyvstop+73,		/* state 36 */
yycrank+292,	0,		0,			/* state 37 */
yycrank+0,	yysvec+38,	yyvstop+77,		/* state 38 */
yycrank+13,	yysvec+36,	0,			/* state 39 */
yycrank+302,	0,		0,			/* state 40 */
yycrank+0,	yysvec+41,	yyvstop+79,		/* state 41 */
yycrank+0,	yysvec+36,	yyvstop+81,		/* state 42 */
0,	0,	0};
struct yywork *yytop = yycrank+359;
struct yysvf *yybgin = yysvec+1;
char yymatch[] = {
00  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,011 ,012 ,01  ,01  ,011 ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
011 ,01  ,'"' ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,'+' ,01  ,'+' ,'.' ,01  ,
'0' ,'0' ,'0' ,'0' ,'0' ,'0' ,'0' ,'0' ,
'0' ,'0' ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,
'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,
'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,
'A' ,'A' ,'A' ,01  ,01  ,01  ,01  ,'A' ,
01  ,'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,
'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,
'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,
'A' ,'A' ,'A' ,01  ,01  ,01  ,01  ,01  ,
0};
char yyextra[] = {
0,0,0,1,0,0,0,0,
0,0,0,0,0,0,0,0,
0};
int yylineno = 1;
#define YYU(x) x
#define NLSTATE yyprevious=YYNEWLINE
char yytext[YYLMAX];
struct yysvf *yylstate[YYLMAX], **yylsp, **yyolsp;
char yysbuf[YYLMAX];
char *yysptr = yysbuf;
int *yyfnd;
extern struct yysvf *yyestate;
int yyprevious = YYNEWLINE;
#ifdef LEXDEBUG
int allprint(int c);
#include <ctype.h>
#endif

int
yylook()
{
   register struct yysvf *yystate, **lsp;
   register struct yywork *yyt;
   struct yysvf *yyz;
   int yych, yyfirst;
   struct yywork *yyr;
#ifdef LEXDEBUG
   int debug;
#endif
   char *yylastch;
   /* start off machines */
#ifdef LEXDEBUG
   debug = 0;
#endif
   yyfirst = 1;
   if (!yymorfg)
      yylastch = yytext;
   else {
      yymorfg = 0;
      yylastch = yytext+yyleng;
      }
   for (;;) {
      lsp = yylstate;
      yyestate = yystate = yybgin;
      if (yyprevious == YYNEWLINE)
         yystate++;
      for (;;) {
#ifdef LEXDEBUG
         if (debug) fprintf(yyout, "state %d\n",yystate-yysvec-1);
#endif
         yyt = yystate->yystoff;
         if (yyt == yycrank && !yyfirst) {
            /* May not be any transitions */
            yyz = yystate->yyother;
            if (yyz == 0) {
               break;
               }
            if (yyz->yystoff == yycrank) {
               break;
               }
            }
         *yylastch++ = yych = input();
         yyfirst = 0;
      tryagain:
#ifdef LEXDEBUG
         if (debug) {
            fprintf(yyout, "char ");
            allprint(yych);
            putchar('\n');
            }
#endif
         yyr = yyt;
         if (yyt > yycrank) {
            yyt = yyr + yych;
            if (yyt <= yytop && yyt->verify+yysvec == yystate) {
               if (yyt->advance+yysvec == YYLERR)
                  /* Error transitions */
                  {unput(*--yylastch);break;}
               *lsp++ = yystate = yyt->advance+yysvec;
               goto contin;
               }
            }
#ifdef YYOPTIM
         else if (yyt < yycrank) {
            /* r < yycrank */
            yyt = yyr = yycrank+(yycrank-yyt);
#ifdef LEXDEBUG
            if (debug) fprintf(yyout,"compressed state\n");
#endif
            yyt += yych;
            if (yyt <= yytop && yyt->verify+yysvec == yystate) {
               if (yyt->advance+yysvec == YYLERR)
                  /* Error transitions */
                  {unput(*--yylastch);break;}
               *lsp++ = yystate = yyt->advance+yysvec;
               goto contin;
               }
            yyt = yyr + YYU(yymatch[yych]);
#ifdef LEXDEBUG
            if (debug) {
               fprintf(yyout,"try fall back character ");
               allprint(YYU(yymatch[yych]));
               putchar('\n');
               }
#endif
            if (yyt <= yytop && yyt->verify+yysvec == yystate) {
               if (yyt->advance+yysvec == YYLERR)
                  /* Error transitions */
                  {unput(*--yylastch);break;}
               *lsp++ = yystate = yyt->advance+yysvec;
               goto contin;
               }
            }
         if ((yystate = yystate->yyother) != NULL &&
             (yyt = yystate->yystoff) != yycrank) {
            goto tryagain;
            }
#endif
         else {
            unput(*--yylastch);
            break;
            }

      contin:
#ifdef LEXDEBUG
      if (debug) {
         fprintf(yyout,"stopped at %d with ",*(lsp-1)-yysvec-1);
         allprint(yych);
         putchar('\n');
         }
#endif
         ;
         }
      while (lsp-- > yylstate) {
         *yylastch-- = 0;
         if (*lsp != 0 && (yyfnd = (*lsp)->yystops) != NULL && *yyfnd > 0) {
            yyolsp = lsp;
            if (yyextra[*yyfnd]) {
               /* Must be backup */
               while (yyback((*lsp)->yystops, -*yyfnd) != 1 && lsp > yylstate) {
                  lsp--;
                  unput(*yylastch--);
                  }
               }
            yyprevious = YYU(*yylastch);
            yylsp = lsp;
            yyleng = yylastch-yytext+1;
            yytext[yyleng] = 0;
#ifdef LEXDEBUG
            if (debug) {
               fprintf(yyout,"\nmatch ");
               fprintf(yyout,"%s",yytext);
               fprintf(yyout," action %d\n",*yyfnd);
               }
#endif
            return(*yyfnd++);
            }
         unput(*yylastch);
         }
      if (yytext[0] == 0) {
         yysptr = yysbuf;
         return(0);
         }
      yyprevious = yytext[0] = input();
      if (yyprevious > 0)
         output(yyprevious);
      yylastch = yytext;
#ifdef LEXDEBUG
      if (debug)
         putchar('\n');
#endif
      }
}

int
yyback(int *p, int m)
{
   if (p == 0)
   return (0);
   while (*p) {
      if (*p++ == m)
         return (1);
      }
   return (0);
}

int
yyinput()
{
   return(input());
}

void
yyoutput(int c)
{
   output(c);
}

void
yyunput(int c)
{
   unput(c);
}
#ifdef LEXDEBUG
allprint(c)
   int c;
{
   switch (c) {
      case 014:
         printf("\\f");
	 break;
      case '\n':
         printf("\\n");
	 break;
      case '\t':
         printf("\\t");
	 break;
      case '\b':
        printf("\\b");
	break;
     case ' ':
        printf("\\\bb");
	break;
     default:
        if (!isprint(c)) {
	   printf("\\%-3o",c);
	   }
	else
	   putchar(c);
	break;
	}
  return;
}
#endif
