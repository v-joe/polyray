/* Test file to create a Targa height field for POVRay -- Alexander Enzmann */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626
#endif

static float **height_buffer, **height_buffer2;
static int GridHeight, GridWidth;

static FILE *TargaFile;

void
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
   tgaheader[16] = 24;
   tgaheader[17] = 0x20;
   fwrite(tgaheader, 18, 1, TargaFile);
}

void
TargaWrite(unsigned char r, unsigned char g, unsigned char b)
{
   fputc(b, TargaFile);
   fputc(g, TargaFile);
   fputc(r, TargaFile);
}

void
TargaClose()
{
   fclose(TargaFile);
}

void
alloc_height_buffers()
{
   int i, j;
   GridHeight = 160;
   GridWidth = 160;

   height_buffer = malloc(GridHeight * sizeof(float *));
   height_buffer2 = malloc(GridHeight * sizeof(float *));
   if (height_buffer == NULL || height_buffer2 == NULL) {
      printf("Failed to allocate height buffer/n");
      exit(1);
      }
   for (i=0;i<GridHeight;i++) {
      height_buffer[i] = malloc(GridWidth * sizeof(float));
      height_buffer2[i] = malloc(GridWidth * sizeof(float));
      if (height_buffer[i] == NULL || height_buffer2[i] == NULL) {
         printf("Failed to allocate element %d of height buffer/n", i);
         exit(1);
         }
      for (j=0;j<GridWidth;j++)
         height_buffer[i][j] = 0.0;
      }
}

void
clear_temp_buffer()
{
   int i, j;
   for (i=0;i<GridHeight;i++)
      for (j=0;j<GridWidth;j++)
         height_buffer2[i][j] = -10000.0;
}

void
add_buffers()
{
   int i, j;
   for (i=0;i<GridHeight;i++)
      for (j=0;j<GridWidth;j++)
         if (height_buffer2[i][j] != -10000.0)
            height_buffer[i][j] = (height_buffer[i][j]+height_buffer2[i][j])/2;
}

void
evaluate_wake(int x0, int z0)
{
   float wave_phase, crest_angle;
   float y;
   int x, z;
   float wave_count = 4.0;
   int wave_steps  = 100;
   float phase_delta = 2.0 * M_PI / wave_steps;
   float total_phase = 2.0 * M_PI * wave_count;
   int crest_steps = 180;
   float amplitude = 64.0;
   float wave_x_scale = 3;
   float wave_z_scale = 1.4;
   int step_count = 0;
   float maxy = -1000.0;
   float miny = 1000.0;
   float decay_factor = 0.8; /* Amount of decay by wavelength */

   for (wave_phase=0.0;wave_phase<=total_phase;wave_phase+=phase_delta) {
      clear_temp_buffer();
printf("Step %d of %d/r", step_count, (int)(wave_count * wave_steps));
      step_count++;

      /* Determine the amplitude of the wave at this phase angle */
      y = amplitude *                          /* Full wave height */
          cos(fmod(wave_phase, 2.0 * M_PI)) *  /* Change with phase angle */
          pow(decay_factor, wave_phase);       /* Decay with distance */

      if (y < miny) miny = y;
      if (y > maxy) maxy = y;

      /* Steps through the equation for wave crests, generating offsets */
      for (crest_angle=-M_PI/2.0+M_PI/((float)crest_steps);
           crest_angle<M_PI/2.0;
           crest_angle += M_PI / ((float)crest_steps)) {
         x = x0 + wave_x_scale * (float)wave_phase *
             (5 * cos(crest_angle) - cos(3 * crest_angle));
         z = z0 + wave_z_scale * (float)wave_phase *
             (sin(crest_angle) + sin(3 * crest_angle));
         
         /* Add the deflection into the buffer */
         if (x >= 0 && x < GridWidth && z >= -GridWidth/2 && z < GridHeight/2) {
            if (height_buffer2[x][z + GridWidth/2] == -10000.0)
               height_buffer2[x][z + GridWidth/2] = y;
            }
         }
      add_buffers();
      }
printf("Max: %f, Min: %f/n", maxy, miny);
}

void
evaluate_wake1(int x0, int z0)
{
   int x, z;
   float xx0, xx1, y, r, t0, t1, theta, wave_phase;
   float k1, k2;
   float wave_number = 1.0;
   float amplitude = 128.0;
   float ship_velocity = 50.0;
   float wave_velocity = 50.0;
   float miny = 10000.0;
   float maxy = -10000.0;
   float decay_factor = 0.8;
   float dist_scale0 = 0.2;
   float dist_scale1 = 0.1;
   
   clear_temp_buffer();
   for (x=1;x<GridHeight;x++) {
      for (z=-GridWidth/2;z<GridWidth/2;z++) {
         xx0 = x - x0;
         xx1 = z - z0;
         theta = atan(xx1 / xx0);
         k1 = cos(theta);
         k2 = sin(theta);
         t0 = xx0 / ship_velocity;
         t1 = fabs(xx1) / wave_velocity;
         if (t1 <= t0) {
            r = sqrt(xx0 * xx0 + xx1 * xx1);
            wave_phase = (k1 * xx0 + k2 * xx1) * dist_scale0;
            /* Determine the amplitude of the wave at this phase angle */
            y = amplitude *        /* Full wave height */
                cos(wave_phase) *  /* Change with phase angle */
                pow(decay_factor, r * dist_scale1);  /* Decay with distance */
            if (y < miny) miny = y;
            if (y > maxy) maxy = y;
            height_buffer2[x][z + GridWidth/2] = y;
            }
         }
      }
   add_buffers();
}

void
dump_height_buffer()
{
   unsigned char r, g, b;
   float height;
   int y, i, j;
   for (i=0;i<GridHeight;i++) {
      for (j=0;j<GridWidth;j++) {
         height = height_buffer[i][j];
         if (height < -128.0) height = -128.0;
         if (height > 127.0) height = 127.0;
         height += 128.0;
         r = height;
         height -= (float)r;
         g = 0; /* (unsigned char)(256.0 * height); */
         b = 0;
         TargaWrite(r, g, b);
         }
      }
}

void
main()
{
   alloc_height_buffers();
   evaluate_wake(0, 0);
   TargaOpen("wave.tga", GridWidth, GridHeight);
   dump_height_buffer();
   TargaClose();
}
