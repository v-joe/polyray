%{
/*

  Copyright (C) 1993-1996, Alexander Enzmann, All rights reserved.
 
  This software may be used for any private and non-commercial
  use.

  You may not distribute this software, in whole or in part,
  for any commercial purpose, without the express consent of
  the authors.
 
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
%}

%%
[ \t\r]                         ;
"/*"                            skipcom();
"//".*$                         ;
\n                              ;
\"[^\"\n]*\"                    { return copy_yystring(); }
[0-9]+([.]?[0-9]*(e[\-\+]?[0-9]+)?)?  { yylval.flt = atof(yytext); return NUM; }
[0-9]*([.]?[0-9]+(e[\-\+]?[0-9]+)?)?  { yylval.flt = atof(yytext); return NUM; }
[A-Za-z0-9_]+                   return lookup();
"&&"                            return AND_SYM;
"||"                            return OR_SYM;
"<="                            return LTEQ_SYM;
"=="                            return EQUAL_SYM;
">="                            return GTEQ_SYM;
"!="                            return NEQUAL_SYM;
.                               return yytext[0];

%%

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
