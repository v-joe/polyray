/*
   Main driver for ray-tracer

  Copyright (C) 1993-1996, Alexander Enzmann, All rights reserved.

  This software may be used for any private and non-commercial
  use.

  You may not distribute this software, in whole or in part,
  for any commercial purpose, without the express consent of
  the authors.

  There is no warranty or other guarantee of fitness of this software
  for any purpose.  It is provided solely "as is".

*/
#include <time.h>

#include "defs.h"
#include "vector.h"
#include "io.h"
#include "memory.h"
#include "psupport.h"
#include "particle.h"
#include "light.h"
#include "sweep.h"
#include "scan.h"
#include "screen.h"
#include "display.h"
#include "pic.h"
#include "symtab.h"
#include "bound.h"
#include "image.h"
#include "intersec.h"

/* Platform/compiler specific include files */
#if defined( MAC )
#include <console.h>
#include <profile.h>
#endif

static int skip_frame = -1;

static void
show_arguments(void)
{
message("Syntax: polyray filename (options)\n");
message("\n File options:\n");
message("     -o outputfile      [ default is 'out.tga'                  ]\n");
message("     -p pixelsize       [ bits/pixel [8, 16, 24, 32]            ]\n");
message("     -d                 [ Render as a depth map                 ]\n");
message("     -u                 [ Turn off RLE compression of output    ]\n");
message("     -x columns         [ trace 'colunms' pixels per row        ]\n");
message("     -y lines           [ trace 'lines' rows                    ]\n");
message("     -B                 [ flush the output file every line      ]\n");
message("     -b count           [ flush every 'count' pixels            ]\n");
message("     -R                 [ Resume an old trace                   ]\n");
message("     -z y0 y1           [ Render from line y0 to line y1        ]\n");
message("\n Rendering options:\n");
message("     -r method          [ Render: 0=Ray,1=Scan,2=Wire,3=Hidden  ]\n");
message("                        [         4=Gourad,5=Raw,6=uv           ]\n");
message("     -q flags           [ Turn on/off various shading options:  ]\n");
message("                        [ 1=shadow, 2=reflect, 4=transmit, 8=two]\n");
message("                        [ sides, 16=check uv, 32=flip normals   ]\n");
message("                        [ 63=all flags                          ]\n");
message("     -N                 [ Don't generate an image file          ]\n");
message("     -M kbytes          [ Max # of KBytes for image buffer      ]\n");
#ifndef unix
#if defined( MAC )
   message("Press mouse to continue\n");
#else
   message("Press any key to continue\n");
#endif
   while (!kbhit()) ;
#if !defined( MAC )
   (void)getch();
#endif
#endif
message("\n");
message("\nStatus options:\n");
message("     -t freqency        [ Status: 0=none,1=totals,2=line,3=pixel]\n");
message("\n Antialiasing options for ray-tracing:\n");
message("     -a mode            [ AA: 0=none,1=filter,2-4=adaptive      ]\n");
message("     -T threshold       [ Threshold to start oversampling       ]\n");
message("     -S samples         [ Use 'samples' rays/pixel              ]\n");
message("\n Optimization options:\n");
message("     -O optimizer       [ 0 = none, 1 = slabs                   ]\n");
message("\n Display option:\n");
message("     -V mode            [ 0=none, 1-5=8bit, 6-10=15bit,         ]\n");
message("                        [ 11-15=16bit, 16-20=24bit, 21-22=4bit  ]\n");
message("     -P pallette        [ 0=grey, 1=332, 2=666, 3=EGA           ]\n");
message("     -e start           [ Start position in VGA pallette        ]\n");
message("     -W                 [ Wait for key before clearing display  ]\n");
message("     -D flag            [ 0=no video dither, 1 = use dither     ]\n");
message("\n Frame counter option:\n");
message("     -F start_frame     [ Skip frames until start_frame         ]\n");
message("\n Abort option:\n");
message("     -Q abort_option    [ 0 = no abort, 1 = check by pixel      ]\n");
message("                        [ 2 = by line/object                    ]\n");
message("\n");
}

/* get option letter from argument vector */
static int optind = 2, /* index into parent argv vector */
           optopt;     /* character checked for validity */
static char *optarg1,
            *optarg2;  /* argument(s) associated with option */

#define BADCH   (int)'?'
#define EMSG   ""

static int
getopt(int nargc, char **nargv, char *ostr)
{
   static char *place = EMSG; /* option letter processing */
   char *oli;                 /* option letter list index */

   if(!*place) {         /* update scanning pointer */
      if (optind >= nargc || *(place = nargv[optind]) != '-' || !*++place)
         return(EOF);
      if (*place == '-') {   /* found "--" */
         ++optind;
         return(EOF);
         }
      }

   /* option letter okay? */
   if ((optopt = (int)*place++) == (int)':' || 
       (optopt == (int)'+') ||
       !(oli = strchr(ostr,optopt))) {
      if(!*place) ++optind;
      warning(": illegal option -- '%c'", (char)optopt);
      return 0;
      }

   if (*++oli == ':') {
      /* need an argument */
      if (*place)
         optarg1 = place;   /* no white space */
      else if (nargc <= ++optind) { /* no arg */
         place = EMSG;
         warning(": option requires an argument -- '%c'", (char)optopt);
         return 0;
      }
       else optarg1 = nargv[optind];   /* white space */
      place = EMSG;
      ++optind;
      }
   else if (*oli == '+') {
      /* need two arguments */
      if (*place)
         optarg1 = place;   /* no white space */
      else if (nargc <= ++optind) { /* no arg */
         place = EMSG;
         warning(": option requires an argument -- '%c'", (char)optopt);
         return 0;
         }
      else
         optarg1 = nargv[optind];   /* white space */
      if (nargc <= ++optind) { /* no arg */
         place = EMSG;
         warning(": option requires two arguments -- '%c'", (char)optopt);
         return 0;
         }
      else
         optarg2 = nargv[optind];
      place = EMSG;
      ++optind;
      }
   else {
      /* don't need argument */
      optarg1 = NULL;
      if (!*place) ++optind;
   }
   return(optopt);         /* dump back option letter */
}

static void
Polygon_Render(Viewpoint *eye, BinTree *Root)
{
   ostackptr objs;
   Object *obj;
   int i;
#ifdef unix
   char tmp[100];
#endif

   /* Read and render the objects */
   if (setjmp(abort_environ) == 0) {
      for (i=0,objs=Root->members.list;objs!=NULL;i++,objs=objs->next) {
         obj = objs->element;
         if (tickflag > 1)
            status("\r%d ", i);
#ifdef unix
         if (start_frame!=end_frame)
           sprintf(tmp,"F%d/%d, p%d/%d", current_frame-start_frame, 
                end_frame-start_frame, i, Root->members.count);
         else
           sprintf(tmp,"p%d/%d",i,Root->members.count);
         SpecialStatus(tmp);
#endif
         if ((Check_Abort_Flag == 2) && kbhit()) {
#if defined( MAC )
            Abort_Flag = 1;
#else
            Abort_Flag = getch();
#endif
            break;
            }
         render_prim(eye, Root, NULL, obj);
         }
      }
   else {
#if defined( MAC )
      Abort_Flag = 1;
#else
      if (kbhit())
         Abort_Flag = getch();
      else
         Abort_Flag = 1;
#endif
      }
}

static void
StartupDisplay(Viewpoint *eye, char *infilename)
{
   if (Display_Flag != 0)
      /* Initialize the display */
      display_init(eye, infilename);
}

static void
Startup_Statistics(char *instr, char *outstr, Viewpoint *eye)
{
   if (tickflag > 0) {
      status("Data: %s, Targa: %s\n", instr, outstr);
      status("Prims = %lu, Lights = %d, Memory: %lu\n",
             Root.members.count, nLights, nMalloc - nFree);
      status("Res: %dx%d, Aspect: %0.4lg, Shade: %x\n",
              eye->view_xres, eye->view_yres,
              eye->view_aspect, Global_Shade_Flag);
      if (end_frame > start_frame || total_frames != 0)
         status("Frame %d in [%d, %d]\n",
                current_frame, start_frame, end_frame);
      }
}

/* Dump out the statistics for this frame. */
static void
PrintStatistics(time_t ptime, time_t ttime, int wait)
{
   if (Display_Flag && Reset_Display_Flag &&
       ((current_frame == end_frame) ||
        (Check_Abort_Flag != 0) && Abort_Flag))
      display_close(wait);
   if (tickflag > 0) {
      status("preprocess time            %-8ld seconds\n", (long)ptime);
      status("tracing time               %-8ld seconds\n", (long)ttime);
      if (Rendering_Method == RAY_TRACING) {
         status("number of eye rays cast:   %-8lu\n", nRays);
         status("number of shadow rays:     %-8lu\n", nShadows);
         if (nReflected > 0)
            status("number of reflected rays:  %-8lu\n", nReflected);
         if (nRefracted > 0)
            status("number of refracted rays:  %-8lu\n", nRefracted);
         if (nTIR > 0)
            status("number of TIR rays:        %-8lu\n", nTIR);
         if (nJittered > 0)
            status("number of jittered rays:   %-8lu\n", nJittered);
         if (totalQueueResets > 0) {
            status("number of queue resets:    %-8lu\n", totalQueueResets) ;
            status("avg number of queues/ray:  %-8g\n", (Flt)nEnqueued /
                                                        (Flt)totalQueueResets) ;
            status("max queue size:            %-8lu\n", maxQueueSize) ;
            }
         status("number of shadow checks:   %-8lu\n", totalShadows);
         status("number of cache hits:      %-8lu\n", totalShadowCaches);
         status("number of bounds checked:  %-8lu\n", nChecked);
         status("number of bounds queued:   %-8lu\n", nEnqueued);
         status("number of objects checked: %-8lu\n", totalQueues);
         }
      }
}

static void
initialize_scan_data(Viewpoint *eye, Pic *pic,
                     int y_start, int y_end,
                     int Start_Line, int End_Line,
                     int alloc_flag)
{
   /* Set up the view transformation */
   if (eye->WS != NULL)
      polyray_free(eye->WS);
   eye->WS = Normalize_View(eye);
#if 0
printf("Y: %d - %d, Strip: %d, %d\n", y_start, y_end, Start_Line, End_Line);
#endif
   /* This is necessary so that drawing commands will work in
      non file generation modes.  Otherwise we could stuff this
      call into the Allocation routine below. */
   Initialize_Clipping(eye, MAX(y_start,Start_Line), MIN(y_end,End_Line));

   if (alloc_flag)
      Allocate_Scan_Buffers(eye, pic, y_start, y_end);
}

/* Between each generated frame, reclaim any memory that was allocated. */
static void
reset_data(void)
{
   /* Reset the global bean counters */
   Initialize_Bean_Counters();

   /* Allocate memory for the various data structures */
   Initialize_Symtab();

   Initialize_Lights();
}

static void
read_defaults(int argc, char **argv, char *outfilename,
              Viewpoint *eye, int *resume,
              int *wait, int *Start_Line, int *End_Line)
{
   static char *opts = "a:b:D:e:M:o:O:p:P:q:Q:r:S:t:T:V:x:y:z+F:BNRduvW";
   int c;
   int line_buffer = 0;

   /* Set defaults for the output file */
   if (!filebaseflag)
      strcpy(outfilebase, "out");
   strcpy(outfilename, "out.tga");
   optind = 2; /* Reset the index into parent argv vector */

   /* Loop through the command line options */
   while ((c = getopt(argc, argv, opts)) != EOF) {
      switch (c) {
      case 'a':
         antialias = atoi(optarg1);
         if (antialias > 4)
            antialias = 4;
         else if (antialias < 0)
            antialias = 0;
         break;
      case 'b':
         if (line_buffer)
            warning("Line buffering already specified\n");
         else {
            buffer_update = 1;
            buffer_size = atoi(optarg1);
            }
         break;
      case 'B':
         if (buffer_update)
            warning("Pixel buffering already specified\n");
         else {
            line_buffer   = 1;
            buffer_update = 1;
            }
         break;
      case 'd':
         DepthRender = 1;
         break;
      case 'D':
         Dither_Flag = atoi(optarg1);
         if (Dither_Flag < 0 || Dither_Flag > 1)
            error("Dither flag must be either 0 (off) or 1 (on)\n");
         break;
      case 'e':
         Pallette_Start = atoi(optarg1);
         if (Pallette_Start < 0 || Pallette_Start > 240)
            error("First entry of pallette must be between 0 and 240\n");
         break;
      case 'F':
         skip_frame = atoi(optarg1);
         break;
      case 'M':
         MaxBufferRAM = 1024L * atoi(optarg1);
         break;
      case 'N':
         File_Generation_Flag = 0;
         break;
      case 'o':
         strcpy(outfilename, optarg1);
         break;
      case 'O':
         Optimizer = atoi(optarg1);
         if (Optimizer < 0 || Optimizer > 1)
            error("Optimizer must be one of: 0 [none], 1 [slabs]\n");
         break ;
      case 'p':
         pixelsize = atoi(optarg1);
         if (pixelsize != 8 && pixelsize != 16 &&
             pixelsize != 24 && pixelsize != 32)
            error("Pixel size must be 8, 16, 24, or 32 bits\n");
         break;
      case 'P':
         Pallette_Flag = atoi(optarg1);
         if (Pallette_Flag < 0 || Pallette_Flag > 3)
            error("Pallette must be one of: 0 [grey], 1 [8-8-4], 2 [6-6-6], 3 [4 bit]\n");
         break;
      case 'q':
         Global_Shade_Flag = atoi(optarg1);
         if (Global_Shade_Flag > ALL_SHADE_FLAGS)
            error("Too many bits set in shade flag (must be less than %d\n",
                  ALL_SHADE_FLAGS);
         break ;
      case 'Q':
         Check_Abort_Flag = atoi(optarg1);
         if (Check_Abort_Flag < 0 || Check_Abort_Flag > 2) {
            warning("Abort flag must be between 0 and 2");
            Check_Abort_Flag = 1;
            }
         break ;
      case 'r':
         Rendering_Method = atoi(optarg1) + RAY_TRACING;
         if (Rendering_Method < 1 || Rendering_Method > LAST_RENDER_MODE)
            error("Rendering must be one of: 0 [raytrace], 1 [scan], 2 [wire], 3 [hidden line], 4 [gourad], 5 [raw tri], 6 [uv tri]\n");
         if (Rendering_Method == RAW_TRIANGLES ||
             Rendering_Method == UV_TRIANGLES ||
             Rendering_Method == CSG_TRIANGLES)
            File_Generation_Flag = 0;
         break;
      case 'R':
         *resume = 1;
         break;
      case 'S':
         maxsamples = atoi(optarg1);
         if (maxsamples <= 0) {
            warning("Samples must be > 0, reset to 4\n");
            maxsamples = 4;
            }
         break;
      case 'T':
         antialias_threshold = atof(optarg1);
         if (antialias_threshold < 0.0) {
            warning("Antialias threshold must be >= 0.0, reset to 0.02\n");
            antialias_threshold = 0.02;
            }
         antialias_threshold *= antialias_threshold;
         break ;
      case 't':
         tickflag = atoi(optarg1);
         if (tickflag == 0) status_flag = 0;
         break;
      case 'u':
         pixel_encoding = 0;
         break;
      case 'V':
         Display_Flag = atoi(optarg1);
         if (Display_Flag < 0) {
            Display_Flag = -Display_Flag;
            Reset_Display_Flag = 0;
            }
         else if (Display_Flag > FIRST_4BIT_MODE+1)
            error("Unsupported display mode");
         else
            Reset_Display_Flag = 1;
         break;
      case 'W':
         *wait = 1;
         break ;
      case 'x':
         eye->view_xres = atoi(optarg1);
         break;
      case 'y':
         eye->view_yres = atoi(optarg1);
         break;
      case 'z':
         *Start_Line = atoi(optarg1);
         *End_Line = atoi(optarg2);
         break;
      }
   }

   if (line_buffer)
      buffer_size = eye->view_xres;
}

static void
set_line_boundaries(Viewpoint *eye, int *Start_Line, int *End_Line)
{
   if (*Start_Line < 0) {
      eye->view_ystart = 0;
      *Start_Line = 0;
      }
   else
      eye->view_ystart = *Start_Line;
   if (*End_Line < 0) {
      eye->view_yend = eye->view_yres;
      *End_Line = eye->view_yres;
      }
   else
      eye->view_yend = *End_Line;

   if (eye->view_ystart > eye->view_yend) {
      warning("Start line of image %d is after end line %d\n",
            eye->view_ystart, eye->view_yend);
      *Start_Line = -1;
      *End_Line   = -1;
      }
   else if (eye->view_ystart >= eye->view_yres) {
      warning("Start line of image %d is after end of image %d\n",
            eye->view_ystart, eye->view_yres);
      *Start_Line = -1;
      *End_Line   = -1;
      }
   else if (eye->view_yend > eye->view_yres) {
      warning("End line of image %d is too high (reset to %d)\n",
              eye->view_yend, eye->view_yres);
      eye->view_yend = eye->view_yres;
      *End_Line = eye->view_yres;
      }
}

static void
render_scene(int argc, char **argv, Viewpoint *eye, char *infilename,
             char *outfilename, char *outfilebase, int resume,
             time_t *atime, int *wait, int Start_Line, int End_Line)
{
   Pic *pic;
   time_t ptime, ttime;
   int y;
   long ram_needed, max_lines;

   /* At last we get to the main driver loop. */
   for (;;) {
      Start_Line = End_Line = -1;

      /* If we are on the second (or later) frame then reset various
         pieces of data */
      if (Parsed_Flag)
         reset_data();

      /* Parse the data file */
      ReadSceneFile(infilename);

      /* Use command line params to override values that appeared
         in either the ini file or the data file */
      read_defaults(argc, argv, outfilename, eye, &resume,
                    wait, &Start_Line, &End_Line);

      /* If the amount of RAM to use hasn't been specified, then
         assume that the entire image buffer should be allocated. */
      if (MaxBufferRAM < 0)
         MaxBufferRAM = 8L * (eye->view_xres + 1) * (eye->view_yres + 1);

      /* Determine the line boundaries for this image */
      set_line_boundaries(eye, &Start_Line, &End_Line);

      /* If the output of this render is a depth file then make sure
         that the file is neither RLE compressed nor antialiased. */
      if (DepthRender == 1) { pixel_encoding = 0; antialias = 0; }

      /* If this is the first pass, then we need to set up the animation
         counters. */
      if (!Parsed_Flag) {
         Parsed_Flag = 1;
         if (end_frame < 0) end_frame = total_frames;
         current_frame = start_frame;
         if (tickflag > 0)
            status("Polyray v1.8a, Copyright (c) 1994 Alexander Enzmann\n");
         }

      /* Allocate and initialize memory, set the global shading flags. */
      if (Rendering_Method == SCAN_CONVERSION ||
          Rendering_Method == HIDDEN_LINE ||
          Rendering_Method == GOURAD_SHADE ||
          (Rendering_Method == WIRE_FRAME && File_Generation_Flag != 0)) {
         /* By default there is no shadow/reflect/transmit
            performed for polygonal scenes */
         if (Global_Shade_Flag & UNSET_SFLAG)
            Global_Shade_Flag = TWO_SIDED_SURFS | NORMAL_CORRECT;
         }
      else if (Global_Shade_Flag & UNSET_SFLAG)
         /* By default check almost everything when raytracing */
         Global_Shade_Flag = SHADOW_CHECK | REFLECT_CHECK | TRANSMIT_CHECK |
                             UV_CHECK | CAST_SHADOW | NORMAL_CORRECT;

      if (Global_Shade_Flag & SHADOW_CHECK)
         Initialize_Light_Caches();

      if (end_frame > start_frame || total_frames != 0)
         build_outfile_name(outfilebase, outfilename);

      InstantiateParticles(eye);

      if (current_frame >= skip_frame) {
         /* This section performs the rendering of a single
            frame.  An inner loop is required to allow for the
            rendering of a strip of the image at a time. */
         StartupDisplay(eye, infilename);

         AddLightObjects(&Root);

         Startup_Statistics(infilename, outfilename, eye);

         time(&ptime);

         /* Start by opening a file to hold the image.  If it already
            exists, then the image characteristics stored in the file
            will be used for the final image. */
         if (File_Generation_Flag)
            pic = TGAOpen(outfilename, eye, resume, pixel_encoding,
                          pixelsize);
         else
            pic = NULL;

         if ((Optimizer > 0) &&
             ((Rendering_Method == RAY_TRACING) ||
              ((Rendering_Method == GOURAD_SHADE ||
                Rendering_Method == SCAN_CONVERSION) &&
               (Global_Shade_Flag &
                  (SHADOW_CHECK | REFLECT_CHECK | TRANSMIT_CHECK))))) {
            /* We are either ray tracing or we have shading options
               in a scan conversion mode that may lead to a recursive
               call to the ray tracer.  In either case if slab
               optimization is turned on we have to create the slabs. */
            BuildBoundingSlabs(&Root);

            /* Create a list of all objects (prims and composites) that
               contain the eye within their bounds */
            AddEyeObjects(Root.slab_root, &Root.eyeprims);
            }

         /* Figure out how many lines at a time we can render */
         if ((Rendering_Method == WIRE_FRAME && File_Generation_Flag == 0) ||
             Rendering_Method == RAW_TRIANGLES ||
             Rendering_Method == UV_TRIANGLES ||
             Rendering_Method == CSG_TRIANGLES) {
            ram_needed = 0;
            max_lines  = eye->view_yres;
            }
         else {
            ram_needed = 8L * (eye->view_xres + 1) * (eye->view_yres + 1);
            max_lines  = MaxBufferRAM / (8L * (eye->view_xres + 1));
            if (max_lines < 3) {
               warning("Too little RAM allocated for rendering: %ld",
                       MaxBufferRAM);
               max_lines = 3;
               }
            else if (max_lines > eye->view_yres)
               max_lines = eye->view_yres + 1;
            }

         /* Create the image in a series of strips */
         for (y=0;y<eye->view_yres;y+=max_lines) {
            if (ram_needed == 0)
               initialize_scan_data(eye, pic, y, y + max_lines,
                                    Start_Line, End_Line, 0);
            else
               initialize_scan_data(eye, pic, y, y + max_lines,
                                    Start_Line, End_Line, 1);

            /* Now go render each strip (as long as it overlaps the part
               of the image that we want to render) */
            if ((y < Start_Line && y+max_lines >= Start_Line) ||
                (y >= Start_Line && y <= End_Line)) {
               if (Rendering_Method == RAY_TRACING)
                  Screen(eye, MAX(y, Start_Line), MIN(y + max_lines, End_Line));
               else
                  Polygon_Render(eye, &Root);
                }

            /* Perform and point/line/spline drawing commands */
            if (Check_Abort_Flag == 0 || Abort_Flag == 0)
               DoDrawing(eye, Draw_Commands);

            /* Blast the contents of this strip out to the file */
            if (File_Generation_Flag)
               TGADump(eye, pic, y, y+max_lines-1);

            if (ram_needed != 0)
               Destroy_Scan_Buffers(eye);

            /* Check for an abort before continuing with the render */
            if (Check_Abort_Flag && Abort_Flag != 0)
               break;
            }

         /* If we were creating an output file then we close it now */
         if (File_Generation_Flag)
            TGAClose(pic);

         time(&ttime);
         ttime -= ptime;
         ptime -= *atime;
         }
      else if (tickflag > 0)
         status("\rSkipping frame: %d    ", current_frame);

      /* Only deallocate everything if this is the last frame in animation */
      if (current_frame < end_frame &&
          (!Check_Abort_Flag || Abort_Flag == 0))
         Deallocate_Symtab(0);
      else
         Deallocate_Symtab(1);

      Deallocate_Lights();

      /* Deallocate the perspective view transform */
      if (eye->WS != NULL) {
         polyray_free(eye->WS);
         eye->WS = NULL;
         }

      if (current_frame > skip_frame)
         PrintStatistics(ptime, ttime, *wait);

#if defined( MAC )
      if (kbhit() && Check_Abort_Flag > 0)
#else
      if (kbhit() && (Check_Abort_Flag > 0 || (Abort_Flag = getch()) == 27))
#endif
         break;
      time(atime);

#if defined( DEBUG_POINTERS ) && 0
   allocation_status();
#endif

      /* Increment the frame counter and see if we are done */
      current_frame++;
      if ((current_frame > end_frame) ||
          (Check_Abort_Flag != 0 && Abort_Flag != 0))
         break;
      }

   /* Deallocate any memory devoted to particle systems */
   FreeParticles();
#if defined( DEBUG_POINTERS ) && 0
   allocation_status();
#endif
}

int
main(int argc, char **argv)
{
   char *infilename = NULL;
   char outfilename[128];
   int resume = 0;
   int wait = 0;
   int Start_Line = -1, End_Line = -1;
   time_t atime;

#if defined( MAC )
   // InitProfile(1000, 200);
   argc = ccommand(&argv);
#endif

   if (argc < 2) {
      show_arguments();
      exit(1);
      }
   else {
      infilename = argv[1];
      }

   /* Set the warning and status log files to stderr */
   SetMessageLog(NULL);

   /* Here is the actual set of calls that does something */
   time(&atime);
   Initialize_Symtab();
   Read_Initialization_Data();
   read_defaults(argc, argv, outfilename, &Eye, &resume,
                 &wait, &Start_Line, &End_Line);
   render_scene(argc, argv, &Eye, infilename, outfilename, outfilebase,
                resume, &atime, &wait, Start_Line, End_Line);

   /* Close the warning and status log files */
   SetMessageLog(NULL);

#if defined( MAC )
   // freopen("profile.out", "w", stdout);
#elif !defined( _WINDOWS )
   if (kbhit()) getch();
#endif

   exit( Abort_Flag ); /* Return the status */

   return 0; /* Just to keep the compiler happy */
}
