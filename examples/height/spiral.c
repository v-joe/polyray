/* Test file to create Targa height fields for POVRay -- Alexander Enzmann */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626
#endif

static float **height_buffer, **height_buffer2;
static int GridHeight, GridWidth;

static FILE *TargaFile;
static int TargaBits = 24;

static void
TargaOpen(char *filename, int x, int y)
{
   unsigned char tgaheader[18];
   if ((TargaFile = fopen(filename, "wb")) == NULL) {
      printf("Failed to open Targa file: %s/n", filename);
      exit(1);
      }
   memset(tgaheader, 0, 18);
   tgaheader[2] = 2;
   tgaheader[12] = (unsigned char)(x & 0xFF);
   tgaheader[13] = (unsigned char)((x >> 8) & 0xFF);
   tgaheader[14] = (unsigned char)(y & 0xFF);
   tgaheader[15] = (unsigned char)((y >> 8) & 0xFF);
   tgaheader[16] = TargaBits;
   tgaheader[17] = 0x20;
   fwrite(tgaheader, 18, 1, TargaFile);
}

static void
TargaWrite(unsigned char r, unsigned char g, unsigned char b)
{
   if (TargaBits == 16) {
      fputc(g, TargaFile);
      fputc(r, TargaFile);
      }
   else if (TargaBits == 24) {
      fputc(b, TargaFile);
      fputc(g, TargaFile);
      fputc(r, TargaFile);
      }
   else {
      fprintf(stderr, "Wrong # of bits/pixel for Targa file: %d/n",
              TargaBits);
      exit(1);
      }
}

static void
TargaClose()
{
   fclose(TargaFile);
}

static void
write_height(float y)
{
   unsigned char r, g, b;
   if (y < -127.0) y = -127.0;
   if (y > 127.0) y = 127.0;
   y += 128.0;
   r = (unsigned char)y;
   y -= (float)r;
   g = (unsigned char)(256.0 * y);
   b = 0;
   TargaWrite(r, g, b);
}

static void
alloc_height_buffer(int width, int height)
{
   int i, j;

   height_buffer = malloc(height * sizeof(float *));
   if (height_buffer == NULL) {
      printf("Failed to allocate height buffer/n");
      exit(1);
      }
   for (i=0;i<height;i++) {
      height_buffer[i] = malloc(width * sizeof(float));
      if (height_buffer[i] == NULL) {
         printf("Failed to allocate element %d of height buffer/n", i);
         exit(1);
         }
      for (j=0;j<width;j++)
         height_buffer[i][j] = 0.0;
      }
}

void
dump_height_buffer(int width, int height)
{
   int i, j;
   for (i=0;i<height;i++)
      for (j=0;j<width;j++)
	 write_height(height_buffer[i][j]);
}

/* Function is: r = c * exp(a * theta) */
static void
make_spiral(int width, int height, int turns, float swidth)
{
   int i, j, k;
   int x, y;
   float fx, fy;
   float theta, r;
   float scale, offset;
   int ang_steps = 8 * 360;
   float delta_theta = 2.0 * M_PI / ((float)ang_steps);
   float a = 1.0;
   float c = 1.0;
   int bump_width = 3;

   /* Figure out how to scale to get the proper # of
      turns into the height field */
   /* r = c * exp(a * 2.0 * M_PI * turns); */
   r = a * 2.0 * M_PI * ((float)turns);
   scale = 0.15 * ((float)width) / (2.0 * r);

   /* Now step through r, plopping height values down as we go */
   for (i=0,theta=0.0;i<turns*ang_steps;i++,theta+=delta_theta) {
      /* r = c * exp(a * theta); */
      r = a * 2.0 * M_PI * theta;
      fx = scale * r * cos(theta);
      fy = scale * r * sin(theta);

      x = (int)fx + width/2;
      y = (int)fy + height/2;

      for (j=-bump_width;j<=bump_width;j++)
	 for (k=-bump_width;k<=bump_width;k++)
	    if (x+j >= 0 && x+j < width &&
	        y+k >= 0 && y+k < height) {
	       offset = 2.0 * sqrt((float)j * (float)j + (float)k * (float)k);
	       height_buffer[x+j][y+k] =
		  (height_buffer[x+j][y+k] + (10.0 - offset)) / 2.0;
	       }
printf("/rT: %g ", theta);
      }
}

void
main(argc, argv)
 int argc;
 char **argv;
{
   int width    = 256;
   int height   = 256;
   int turns    = 5;
   float swidth = 5.0;

   TargaOpen("spirhf.tga", height, width);
   alloc_height_buffer(width, height);
   make_spiral(width, height, turns, swidth);
   dump_height_buffer(width, height);
   TargaClose();
}
