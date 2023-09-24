/* Test file to create Targa height fields for POVRay -- Alexander Enzmann */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

/* Noise generator variables */
#define LARGE 100000.0
#define HASH_SIZE 78737
#define spline(x) ((3 - 2 * (x)) * (x) * (x))
static int mult_table_size = 5;
static unsigned long mult_table[] = { 99137, 113713, 55237, 994472, 441973 };

static float **height_buffer;
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

static float
hash3d(unsigned long x, unsigned long y, unsigned long z)
{
   int i;
   unsigned long K, Kt, P;
   float result;
   K = ((x & 0x03ffL) << 20) |
       ((y & 0x03ffL) << 10) |
        (z & 0x03ffL);
   Kt = 0;
   for (i=0;i<mult_table_size;i++)
      Kt = ((Kt + K) * mult_table[i]) % HASH_SIZE;
   result = (float)Kt / (float)(HASH_SIZE - 1);
   return result;
}

static float
flt_noise(float x, float y, float z)
{
   unsigned long ix, iy, iz, jx, jy, jz;
   float sx0, sy0, sz0, sx1, sy1, sz1;
   float result;

   x += LARGE;   y += LARGE;  z += LARGE;
   ix = x;      iy = y;      iz = z;
   jx = ix + 1; jy = iy + 1; jz = iz + 1;
   
   /* Compute cubic interpolated influences of surrounding
      lattice points */
   sx0 = spline(x - ix);
   sy0 = spline(y - iy);
   sz0 = spline(z - iz);
   sx1 = 1.0 - sx0; sy1 = 1.0 - sy0; sz1 = 1.0 - sz0;

   result  = hash3d(ix, iy, iz) * sx1 * sy1 * sz1;
   result += hash3d(ix, iy, jz) * sx1 * sy1 * sz0;
   result += hash3d(ix, jy, iz) * sx1 * sy0 * sz1;
   result += hash3d(ix, jy, jz) * sx1 * sy0 * sz0;
   result += hash3d(jx, iy, iz) * sx0 * sy1 * sz1;
   result += hash3d(jx, iy, jz) * sx0 * sy1 * sz0;
   result += hash3d(jx, jy, iz) * sx0 * sy0 * sz1;
   result += hash3d(jx, jy, jz) * sx0 * sy0 * sz0;
   return result;
}

static float
noise3d(float x, float y, float z,
        float recursion_scale, int count)
{
   float x1, y1, z1;
   float result = 0.0;
   float scale  = 1.0;
   float iscale = 1.0 / recursion_scale;
   float magnitude = 0.0;
   int i;

   x1 = x;
   y1 = y;
   z1 = z;
   for (i=0;i<count;i++) {
     result += scale * flt_noise(x1, y1, z1);
     magnitude += scale;
     if (i < count-1) {
        x1 *= 2.0;
        y1 *= 2.0;
        z1 *= 2.0;
        scale *= recursion_scale;
        }
     }
   result /= magnitude;
   return result;
}

static write_height(float y)
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
make_noise_field(int width, int height,
                 float xlow, float xhigh,
                 float zlow, float zhigh,
                 float seed,
                 float recursion_scale, int octaves)
{
   float x, xdelta;
   float y;
   float z, zdelta;
   int i, j;

   xdelta = (xhigh - xlow) / ((float)(width-1));
   zdelta = (zhigh - zlow) / ((float)(height-1));

   /* Loop through all possible values of z and x, generating
      altitudes as we go. */
   for (i=0,z=zlow;i<height;i++,z+=zdelta) {
      printf("%d  /r", i);
      for (j=0,x=xlow;j<width;j++,x+=xdelta) {
         /* Calculate the altitude */
         y = noise3d(x, seed, z, recursion_scale, octaves);
         
         /* Scale it to fit into -128 -> 128 */
         y = 256.0 * (y - 0.5);
         
         /* Write out the altitude */
         write_height(y);
         }
      }
}

void
main(argc, argv)
 int argc;
 char **argv;
{
   int h = 256;       /* height of resulting height field */
   int w = 256;       /* width of resulting height field */
   int oct = 3;       /* # of octaves of noise to use */
   float x0 = -8.0;   /* Low x to interpolate noise over */
   float x1 =  8.0;   /* High x to interpolate noise over */
   float z0 = -8.0;   /* Low z to interpolate noise over */
   float z1 =  8.0;   /* High z to interpolate noise over */
   float seed = 1.11; /* Seed to change resulting pattern */
   float sc =  0.25;  /* Contribution of each octave of noise */

   TargaOpen("fracthf.tga", h, w);
   make_noise_field(h, w, x0, x1, z0, z1, seed, sc, oct);
   TargaClose();
}
