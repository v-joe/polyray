#include <Palettes.h>
#include <console.h>
#include "defs.h"
#include "display.h"

int Xresolution = 512;
int Yresolution = 512;
int Pallette_Start = 0; /* First actual entry in pallette */
int Pallette_Flag = 1;  /* Use 884 pallette */
int Display_Flag = 0;
int Reset_Display_Flag = 1;  /* Turn on the requested video mode */

/* Define the top right, the width, and the length of usuable display area */
int Display_x0 = -1, Display_y0 = -1;
int Display_xl = -1, Display_yl = -1;
int Dither_Flag = 1;

static WindowPtr cwin; // cWindow;
static PaletteHandle PolyPalette;
static unsigned offx = 0;
static unsigned offy = 0;
static int maxx = 320;
static int maxy = 240;
static unsigned screen_maxx = 320;
static unsigned screen_maxy = 240;
static int display_active_flag = 0;
static Flt X_Display_Scale = 1.0;
static Flt Y_Display_Scale = 1.0;

static int EGA_colors[16][3] =
   {{  0,  0,  0}, {  0,  0,255}, {  0,255,  0}, {  0,127,127},
    {255,  0,  0}, {127,  0,127}, {127,127,  0}, { 85, 85, 85},
    {170,170,170}, {127,127,255}, {127,255,127}, {  0,255,255},
    {255,127,127}, {255,  0,255}, {255,255,  0}, {255,255,255}};
static int EGA_remap[16] =
   {0,  6,  5,  4,  3,  2, 1,  7,
    8, 14, 13, 12, 11, 10, 9, 15};
static int Dither_Matrix[8][8] =
   {{ 0, 32,  8, 40,  2, 34, 10, 42},
    {48, 16, 56, 24, 50, 18, 58, 26},
    {12, 44,  4, 36, 14, 46,  6, 38},
    {60, 28, 52, 20, 62, 30, 54, 22},
    { 3, 35, 11, 43,  1, 33,  9, 41},
    {51, 19, 59, 27, 49, 17, 57, 25},
    {15, 47,  7, 39, 13, 45,  5, 37},
    {63, 31, 55, 23, 61, 29, 53, 21}};

static void
ToolBoxInit()
{
   InitGraf(&thePort);
   InitFonts();
   FlushEvents(0xFFFF, 0);
   InitWindows();
   InitMenus();
   TEInit();
   InitDialogs(0L);
   InitCursor();
}

static void
display_clear()
{
   Rect re;

   PmForeColor(0xFF);
   re.top = 0;
   re.left = 0;
   re.right = maxx;
   re.bottom = maxy;
   PaintRect(&re);
}

void
display_init(Viewpoint *eye, char *S)
{
   RGBColor c;
   int index, i, j;
   int r, g, b;
   Rect re;
   int sflag, xres, yres, x1, y1;

   sflag = 0; /* No subscreen window */
   xres = eye->view_xres;
   yres = eye->view_yres;
   if (Display_x0 < 0) Display_x0 = 0; else sflag = 1;
   if (Display_y0 < 0) Display_y0 = 0; else sflag = 1;
   if (Display_xl < 0) Display_xl = xres; else sflag = 1;
   if (Display_yl < 0) Display_yl = yres; else sflag = 1;

   if (display_active_flag) {
      display_clear();
      return;
      }

   display_active_flag = 1;

   ToolBoxInit();

   Xresolution = xres;
   Yresolution = yres;

   re.top = 50;
   re.left = 300;
   re.right = 300 + maxx;
   re.bottom = 50 + maxy;
   cwin = NewCWindow(0L, &re, "\pPolyray Display",
                        FALSE, noGrowDocProc,
                        (WindowPtr)-1L, FALSE, 0L);
   PolyPalette = NewPalette(256, 0L, pmTolerant, 0x0000);
   c.red = c.green = c.blue = 0;
   SetEntryColor(PolyPalette, 0, &c);

   if (Pallette_Flag == 0) {
      for (i=0;i<128;i++) {
         c.red   = i << 9;
         c.green = i << 9;
         c.blue  = i << 9;
         SetEntryColor(PolyPalette, i, &c);
         }
      }
   else if (Pallette_Flag == 1) {
      i = 0;
      for (r=0;r<8;r++)
         for (g=0;g<8;g++)
            for (b=0;b<4;b++) {
            c.red   = r << 13;
            c.green = g << 13;
            c.blue  = b << 14;
            SetEntryColor(PolyPalette, i, &c);
                i++;
            }
      }
   else if (Pallette_Flag == 2) {
      i = 0;
      for (r=0;r<6;r++)
         for (g=0;g<6;g++)
            for (b=0;b<6;b++) {
               c.red   = (r * 51) << 8;
               c.green = (g * 51) << 8;
               c.blue  = (b * 51) << 8;
               SetEntryColor(PolyPalette, i, &c);
                   i++;
               }
      }
   else if (Pallette_Flag == 3) {
      i = 0;
      for (i=0;i<16;i++) {
         c.red   = EGA_colors[i][0] << 8;
         c.green = EGA_colors[i][1] << 8;
         c.blue  = EGA_colors[i][2] << 8;
         SetEntryColor(PolyPalette, i, &c);
         }
      }

   /* Make sure the last entry is true white */
   c.red   = 0xFFFF;
   c.green = 0xFFFF;
   c.blue  = 0xFFFF;
   SetEntryColor(PolyPalette, i-1, &c);
   SetEntryColor(PolyPalette, 0xFF, &c);

   SetPalette(cwin, PolyPalette, TRUE);
   /* Now that the palette has been created, open the window. */
   ShowWindow(cwin);
   SetPort(cwin);

   display_clear();

   /* Do some conditioning on the screen window to ensure it doesn't cause
      pixels to be generated off screen */
   x1 = Display_x0 + Display_xl;
   y1 = Display_y0 + Display_yl;
   if (Display_x0 < 0) { Display_x0 = 0; }
   if (Display_y0 < 0) { Display_y0 = 0; }
   if (x1 >= screen_maxx) { Display_xl = screen_maxx - Display_x0; }
   if (y1 >= screen_maxy) { Display_yl = screen_maxy - Display_y0; }

   /* If there is a rendering subwindow of the screen that we should
      be working in.  If so, then set the scaling values appropriately. */
   if (sflag) {
      offx = Display_x0;
      offy = Display_y0;
      maxx = Display_x0 + Display_xl;
      maxy = Display_y0 + Display_yl;
      }
   else {
      offx = 0;
      offy = 0;
      maxx = screen_maxx;
      maxy = screen_maxy;
      }

#if 0
   if (xres > Display_xl)
      X_Display_Scale = (Flt)Display_xl / (Flt)xres;
   else
      X_Display_Scale = 1.0;
   if (yres > Display_yl)
      Y_Display_Scale = (Flt)Display_yl / (Flt)yres;
   else
      Y_Display_Scale = 1.0;
   if (X_Display_Scale < Y_Display_Scale)
      Y_Display_Scale = X_Display_Scale;
   else if (Y_Display_Scale < X_Display_Scale)
      X_Display_Scale = Y_Display_Scale;
#else
   X_Display_Scale = (Flt)Display_xl / (Flt)xres;
   Y_Display_Scale = (Flt)Display_yl / (Flt)yres;
   if (X_Display_Scale < Y_Display_Scale)
      Y_Display_Scale = X_Display_Scale;
   else if (Y_Display_Scale < X_Display_Scale)
      X_Display_Scale = Y_Display_Scale;
#endif

   return;
}

void
display_close(int wait_flag)
{
   if (wait_flag) {
      while (!kbhit()) ;
      }
}

static int
determine_color_index(int x, int y, Vec color)
{
   float d;
   int i;
   Vec color2;
   unsigned char r, g, b;

   if (Dither_Flag) {
      color2[0] = color[0] + Dither_Matrix[x & 0x07][y & 0x07] / 255.0;
      color2[1] = color[1] + Dither_Matrix[x & 0x07][y & 0x07] / 255.0;
      color2[2] = color[2] + Dither_Matrix[x & 0x07][y & 0x07] / 255.0;
      }
   else {
      VecCopy(color, color2)
      }

   i = 255.0 * color2[2];
   if (i<0) i=0;
   else if (i>=256) i = 255;
   b = (unsigned char)i;

   i = 255.0 * color2[1];
   if (i<0) i=0;
   else if (i>=256) i = 255;
   g = (unsigned char)i;

   i = 255.0 * color2[0];
   if (i<0) i=0;
   else if (i>=256) i = 255;
   r = (unsigned char)i;

   if (Pallette_Flag == 0)
      /* Greyscale */
      return (r > g ? (r > b ? r >> 1 : b >> 1) : (g > b ? g >> 1 : b >> 1));
   else if (Pallette_Flag == 1)
      /* 884 */
      return (r & 0xE0) | ((g & 0xE0) >> 3) | (b >> 6);
   else if (Pallette_Flag == 2)
      /* 666 */
      return 6 * (6 * (r / 51) + (g / 51)) + b / 51;
   else if (Pallette_Flag == 3) {
      /* EGA colors */
      d = 0.5;
      r = 0;
      if (color2[0] + color2[1] - color2[2] > d)
         r |= 0x01;
      if (color2[0] - color2[1] + color2[2] > d)
         r |= 0x02;
      if (-color2[0] + color2[1] + color2[2] > d)
         r |= 0x04;
      d = 1.5;
      g = 0x08;
      if (r == 0) {
         d = 0.5;
         g = 0x07;
         }
      else if (r == 7) {
         d = 2.5;
         r = 0x08;
         g = 0x07;
         }
      if (color2[0] + color2[1] + color2[2] > d)
         r |= g;
      return EGA_remap[r];
      }
   else {
      return 6 * (6 * (r / 51) + (g / 51)) + b / 51;
      }
}

static void
clip_display_point(int *x, int *y)
{
   if (*x < offx)
      *x = offx;
   else if (*x >= maxx)
      *x = maxx - 1;
   if (*y < offy)
      *y = offy;
   else if (*y >= maxy)
      *y = maxy - 1;
}

static void
putbox(int x0, int y0, int x1, int y1, int color)
{
   Rect r;

   r.top = y0;
   r.left = x0;
   r.bottom = y1+1;
   r.right = x1+1;

   PmForeColor(color);
   PaintRect(&r);
   return;
}

void
display_box(int x0, int y0, int x1, int y1, Vec color)
{
   int xt0, yt0, xt1, yt1;

   yt0 = (int)((Flt)offy + Y_Display_Scale * (Flt)y0);
   xt0 = (int)((Flt)offx + X_Display_Scale * (Flt)x0);
   yt1 = (int)((Flt)offy + Y_Display_Scale * (Flt)y1);
   xt1 = (int)((Flt)offx + X_Display_Scale * (Flt)x1);

   clip_display_point(&xt0, &yt0);
   clip_display_point(&xt1, &yt1);

   putbox(xt0, yt0, xt1, yt1, determine_color_index(x0, y0, color));
}

static void
putpixel(x, y, color)
   int x, y, color;
{
   Rect r;

   if (Xresolution > 80 || Yresolution > 60) {
      r.top = y;
      r.left = x;
      r.bottom = y+1;
      r.right = x+1;
      }
   else {
      r.top = 2*y;
      r.left = 2*x;
      r.bottom = 2*(y+1);
      r.right = 2*(x+1);
      }
   PmForeColor(color);
   PaintRect(&r);
   return;
}

void
display_plot(int x, int y, Vec color)
{
   int x0, y0, x1, y1;

   y0 = (int)((Flt)offy + Y_Display_Scale * (Flt)y);
   x0 = (int)((Flt)offx + X_Display_Scale * (Flt)x);
   y1 = (int)((Flt)offy + Y_Display_Scale * (Flt)(y + 1));
   x1 = (int)((Flt)offx + X_Display_Scale * (Flt)(x + 1));

   clip_display_point(&x0, &y0);
   clip_display_point(&x1, &y1);

   if (x0 == x1 && y0 == y1)
      putpixel(x0, y0, determine_color_index(x, y, color));
   else
      putbox(x0, y0, x1, y1, determine_color_index(x0, y0, color));
}

void
display_line(int x0, int y0, int x1, int y1, Vec color)
{
   Flt xt, yt;
   int i, color_index;

   color_index = determine_color_index(x0, y0, color);

   // Scale from image size to actual screen pixel size
   yt = Y_Display_Scale * y0;
   xt = X_Display_Scale * x0;
   if (xt < 0.0)
      x0 = 0;
   else if (xt > maxx) {
      x0 = maxx - 1;
      }
   else x0 = (int)xt;
   if (yt < 0.0)
      y0 = 0;
   else if (yt > maxy) {
      y0 = maxy;
      }
   else
      y0 = (int)yt;

   yt = Y_Display_Scale * y1;
   xt = X_Display_Scale * x1;
   if (xt < 0.0)
      x1 = 0;
   else if (xt > maxx) {
      x1 = maxx - 1;
      }
   else x1 = (int)xt;
   if (yt < 0.0)
      y1 = 0;
   else if (yt > maxy) {
      y1 = maxy;
      }
   else
      y1 = (int)yt;


   PmForeColor(color_index);
   MoveTo(x0, y0);
   LineTo(x1, y1);
}

