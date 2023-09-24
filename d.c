typedef unsigned char byte;

#include <assert.h>
#include <stdio.h>
#include <time.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#ifdef MIT_SHM
#include <sys/ipc.h>
#include <sys/shm.h>
#include <X11/extensions/XShm.h>
#endif
#include "d.xbm"

Display		*theDisplay = NULL;
int		theScreen;
int		theDepth, bits_per_RGB;
unsigned long	theBlackPixel;
unsigned long	theWhitePixel;
XEvent		theEvent;
Window		theWindow, openWindow();
GC		theGC;

Visual		*theVisual;
XStandardColormap theStdcmap;
XImage		*theImage;
unsigned long	theCmap;
int		theCells;

Atom		wm_delw;
XColor		*theColors;
int		visible=0,truecolor=0;

byte		*bitmap,*bitmapPos;

int		rowlen, planelen;

unsigned int	rdiv[256],rmod[256],rmul[256],
		gdiv[256],gmod[256],gmul[256],
		bdiv[256],bmod[256],bmul[256];

void bestColor(Display *dpy)
{
  XVisualInfo template,*trues,*pseudos,*besttrue,*bestpseudo,*vl;
  int ok,npseudos,ntrues;
  unsigned int i;
  
  /* get all truecolor and pseudocolor visuals */
  template.class = TrueColor;
  trues = XGetVisualInfo(dpy, VisualClassMask, &template, &ntrues);

  template.class = PseudoColor;
  pseudos = XGetVisualInfo(dpy, VisualClassMask, &template, &npseudos);

  /* find the deepest one */
  besttrue = NULL;
  for(i=0;i<ntrues;i++) {
    if(!(besttrue)) besttrue = &trues[i];
    if(trues[i].depth > besttrue->depth) besttrue = &trues[i];
  }

  bestpseudo = NULL;
  for(i=0;i<npseudos;i++) {
    if(!(bestpseudo)) bestpseudo = &pseudos[i];
    if(pseudos[i].depth > bestpseudo->depth) bestpseudo = &pseudos[i];
  }

  if(besttrue) {
    if(bestpseudo) {
     /* If the deepest TrueColor is deeper than the deepest PseudoColor
	we imply that this is a native TrueColor server.
	Even native TrueColor servers sometimes state to be Pseudo for
	compatibility reasons. */
      vl = ((bestpseudo->depth >= besttrue->depth) ? bestpseudo : besttrue);
    } else { vl = besttrue; }
  } else {
    vl = bestpseudo;
    if(!(vl)) {
      fprintf(stderr,"Neither a TrueColor nor a PseudoColor visual found");
      exit(1);
    }
  }

  /* get the values for later use */
  theVisual = vl->visual;
  theScreen = vl->screen;
  theDepth  = vl->depth;
  bits_per_RGB = (theDepth==24?32:theDepth);

  if(vl->class == TrueColor) {
    unsigned long msk,mul;

    truecolor = 1;

    fprintf(stderr,"Using %d bit TrueColor visual\n",vl->depth);

    theStdcmap.colormap = XCreateColormap(dpy,RootWindow(dpy,vl->screen),
					  vl->visual,AllocNone);

    /* convert the masks to max/mult couples for pixel generation */
    msk = vl->red_mask;
    mul = 1;
    while(!(msk & 1)) {
      mul *= 2;
      msk >>= 1;
    }
    theStdcmap.red_max = msk;
    theStdcmap.red_mult = mul;

    msk = vl->green_mask; mul = 1;
    while(!(msk & 1)) { mul *= 2; msk >>= 1; }
    theStdcmap.green_max = msk;
    theStdcmap.green_mult = mul;

    msk = vl->blue_mask; mul = 1;
    while(!(msk & 1)) { mul *= 2; msk >>= 1; }
    theStdcmap.blue_max = msk;
    theStdcmap.blue_mult = mul;
    theStdcmap.base_pixel = 0;
    theStdcmap.visualid = vl->visualid;
    theStdcmap.killid = 0;
  } else {
    truecolor = 0;

    fprintf(stderr, "Using on depth %d PseudoColor visual (",vl->depth);

    /* We NEED a stdcmap on PseudoColor screens */
    ok=XGetStandardColormap(dpy, RootWindow(dpy,vl->screen),
			    &theStdcmap, XA_RGB_BEST_MAP);
    if(ok) {
      fprintf(stderr,"RGB_BEST_MAP).\n");
    } else {
      fprintf(stderr,"custom colormap).\n");
      theStdcmap.colormap = 0;
      /* we should allocate the colormap here instead of only marking
	 this job to be done later. But this weird hack of code doesn't
	 have the xcolors vector allocated yet... */
    }
  }
  /* all info is retrieved -- free the lists */
  XFree(trues);
  XFree(pseudos);
  for (i=0;i<256;i++) {
    rdiv[i]=(i*theStdcmap.red_max)>>8;
    rmod[i]=((i*theStdcmap.red_max)>>2) & 63;
    rmul[i]=i*theStdcmap.red_mult;
    gdiv[i]=(i*theStdcmap.green_max)>>8;
    gmod[i]=((i*theStdcmap.green_max)>>2) & 63;
    gmul[i]=i*theStdcmap.green_mult;
    bdiv[i]=(i*theStdcmap.blue_max)>>8;
    bmod[i]=((i*theStdcmap.blue_max)>>2) & 63;
    bmul[i]=i*theStdcmap.blue_mult;
  }
}

unsigned long rgb_to_pixel (byte Red, byte Green, byte Blue)
{
  return theStdcmap.base_pixel+rmul[rdiv[Red]]+gmul[gdiv[Green]]+
			       bmul[bdiv[Blue]];

 /* Compute the Pixel from the information about the pixel anatomy */
 /* that was either built by hand or got from RGB_BEST_MAP         */

/*
register unsigned long color,red,green,blue;
 red = (unsigned long)Red;
 green = (unsigned long)Green;
 blue = (unsigned long)Blue;
 color = theStdcmap.base_pixel
  +((red   *(theStdcmap.red_max+1)  )>>8) * theStdcmap.red_mult
  +((green *(theStdcmap.green_max+1))>>8) * theStdcmap.green_mult
  +((blue  *(theStdcmap.blue_max+1) )>>8) * theStdcmap.blue_mult;
 return(color);
*/
}

/* Init all basic X related stuff */
void initX()
{
  theDisplay = XOpenDisplay(NULL);
  if (theDisplay == NULL) {
    fprintf(stderr,"ERROR: Cannot establish a connection to the X server %s\n",
            XDisplayName(NULL));
    exit(1);
  }
  bestColor(theDisplay);
  theCmap = theStdcmap.colormap;
  if(!truecolor) {
    theCells = DisplayCells(theDisplay,theScreen);
    if(theCells!=256) fprintf(stderr,
			     "%d color cells reported\n",theCells);
    theWhitePixel = WhitePixel(theDisplay, theScreen);
    theBlackPixel = BlackPixel(theDisplay, theScreen);
  } else {
    theWhitePixel = rgb_to_pixel(255,255,255);
    theBlackPixel = rgb_to_pixel(0,0,0);
  }
}

#ifdef MIT_SHM
XShmSegmentInfo shminfo;
char Shm;
#endif

int ytab[4096];
int xshift=0;

/* special color related initializations */
void initColor(int width,int height)
{
  int i;

#ifdef MIT_SHM
  if (XShmQueryExtension(theDisplay)) {
    theImage = XShmCreateImage(theDisplay,theVisual,
		theDepth,ZPixmap,0,&shminfo,width,height);
    shminfo.shmid = shmget( IPC_PRIVATE, theImage->bytes_per_line *
		theImage->height, IPC_CREAT|0777);
    shminfo.shmaddr = theImage->data = shmat( shminfo.shmid,0,0);
    shminfo.readOnly = False;
    if (XShmAttach(theDisplay,&shminfo)) {
      fprintf(stderr,"Shared Memory enabled\n");
      Shm=1;
    } else {
      fprintf(stderr,"Shared Memory disabled\n");
      XDestroyImage(theImage);
      shmdt(shminfo.shmaddr);
      shmctl(shminfo.shmid,IPC_RMID,0);
      goto NoShm;
    }
  } else {
NoShm:
    Shm=0;
    theImage = XCreateImage(theDisplay,theVisual,
		theDepth,ZPixmap,0,bitmap,width,height,
		32, 0);
    theImage->data = (unsigned char *) malloc(theImage->bytes_per_line*height );
  }
#else
  theImage = XCreateImage(theDisplay,theVisual,
              theDepth,ZPixmap,0,bitmap,width,height,
              32, 0);
  theImage->data = (unsigned char *) malloc(theImage->bytes_per_line*height );
#endif
  ytab[0]=0;
  for (i=1;i<theImage->height;i++)
    ytab[i]=ytab[i-1]+theImage->bytes_per_line;
  switch (theImage->bits_per_pixel) {
    case 32: xshift++;
    case 15: case 16: xshift++;
  }

  if (visible) XFillRectangle(theDisplay, theWindow, theGC, 0,0,width,height);
}

time_t last_flush=0;

Window openWindow(int x,int y,int width,int height,int flag,GC *theNewGC,
		  char *s)
{
  XSetWindowAttributes	theWindowAttributes;
  XSizeHints		theSizeHints;
  unsigned long		theWindowMask;
  Window		theNewWindow;
  Pixmap		theIconPixmap;
  XWMHints		theWMHints;
  static XClassHint	theClassHint  = { "render", "Render" };

  theWindowAttributes.border_pixel	= theWhitePixel;
  theWindowAttributes.background_pixel	= theBlackPixel;
  theWindowAttributes.colormap		= theCmap;
  theWindowAttributes.event_mask	= 
    ButtonPressMask | KeyPressMask | ExposureMask | VisibilityChangeMask;
  theWindowMask = CWBackPixel | CWBorderPixel | CWEventMask | CWColormap;

  theNewWindow = XCreateWindow( theDisplay,
				RootWindow(theDisplay,theScreen),
				x,y,
				width, height,
				2,
				theDepth,
				InputOutput,
				theVisual,
				theWindowMask,
				&theWindowAttributes);

  theIconPixmap = XCreateBitmapFromData(theDisplay,
                                        theNewWindow,
                                        d_bits,
                                        d_width,
                                        d_height);

  theWMHints.icon_pixmap        = theIconPixmap;
  theWMHints.initial_state      = NormalState;
  theWMHints.flags              = IconPixmapHint | StateHint;

  XSetWMHints(theDisplay,theNewWindow,&theWMHints);

  theSizeHints.flags            = PPosition | PSize | PMinSize | PMaxSize;
  theSizeHints.x                = x;
  theSizeHints.y                = y;
  theSizeHints.width            = width;
  theSizeHints.height           = height;
  theSizeHints.min_width        = width;
  theSizeHints.min_height       = height;
  theSizeHints.max_width        = width;
  theSizeHints.max_height       = height;

  XSetNormalHints(theDisplay,theNewWindow,&theSizeHints);

  XSetClassHint(theDisplay,theNewWindow,&theClassHint);

  if (createGC(theNewWindow, theNewGC) == 0) {
    XDestroyWindow(theDisplay, theNewWindow);
    return((Window) 0);
  }

  XStoreName(theDisplay, theNewWindow, s);
  XSetIconName(theDisplay, theNewWindow, s);

  /* Its somewhat complicated to do it by hand. But the Xt Philosophy
     doesn't match the PoV previewer interface in any way */

  wm_delw = XInternAtom(theDisplay,"WM_DELETE_WINDOW",False);
  XSetWMProtocols(theDisplay,theNewWindow,&wm_delw,1);

  XMapWindow(theDisplay,theNewWindow);
  last_flush=time(NULL);
  XFlush(theDisplay);
  return(theNewWindow);
}

int createGC(Window theNewWindow, GC *theNewGC)
{
  XGCValues theGCValues;
  *theNewGC = XCreateGC(theDisplay, theNewWindow, (unsigned long) 0, 
			&theGCValues);

  if (*theNewGC == 0) {
    return(0);
  } else {
    XSetForeground(theDisplay, *theNewGC, theWhitePixel);
    XSetBackground(theDisplay, *theNewGC, theBlackPixel);
    return(1);
  }
}

int error_handler(Display *dpy,XErrorEvent *ev)
{
  fprintf(stderr,"X error received\n"); /* do no quitting */
  return(0);
}

void InitRegionTable();

char inited=0;

void InitDisplay(int width,int height,char *s)
{
  if (inited) return; inited=1;
  initX();
  InitRegionTable();
  if (!truecolor) {
    theColors = (XColor *) malloc(sizeof(XColor) * theCells);
    if (theStdcmap.colormap==0) {
      int rgbbits=theDepth/3,rgbmod=theDepth%3,rbits,gbits,bbits;
      int cnt,mul1=1,mul2=1,i,j,k,l;

      rbits = rgbbits; if(rgbmod) {rbits++; rgbmod--;}
      gbits = rgbbits; if(rgbmod) gbits++;
      bbits = rgbbits;
      theStdcmap.red_mult = 1;
      theStdcmap.base_pixel = 0;
      cnt = rbits; mul1 = mul2 = 1;
      while(cnt) {mul1 *= 2; mul2 *= 2; cnt--;}
      theStdcmap.red_max = mul1-1;
      theStdcmap.green_mult = mul2;
      cnt = gbits; mul1 = 1;
      while(cnt) {mul1 *=2; mul2 *= 2; cnt--;}
      theStdcmap.green_max = mul1-1;
      theStdcmap.blue_mult = mul2;
      cnt = bbits; mul1 = 1;
      while(cnt) {mul1 *=2; cnt--;}
      theStdcmap.blue_max = mul1-1;
      l = theStdcmap.base_pixel;
      for(i=0;i<=theStdcmap.blue_max;i++) {
        for(j=0;j<=theStdcmap.green_max;j++) {
	  for(k=0;k<=theStdcmap.red_max;k++) {
	   theColors[l].red   = 0xffff * k / theStdcmap.red_max;
	   theColors[l].green = 0xffff * j / theStdcmap.green_max;
	   theColors[l].blue  = 0xffff * i / theStdcmap.blue_max;
	   theColors[l].flags = DoRed|DoGreen|DoBlue;
	   theColors[l].pixel = l;
	   l++;
	  }
        }
      }
      theStdcmap.colormap = XCreateColormap(theDisplay,
          RootWindow(theDisplay,theScreen), theVisual,AllocAll);
      XStoreColors(theDisplay, theStdcmap.colormap, theColors, theCells);
      theCmap = theStdcmap.colormap;
      for (i=0;i<256;i++) {
	rdiv[i]=(i*theStdcmap.red_max)>>8;
	rmod[i]=((i*theStdcmap.red_max)>>2) & 63;
	rmul[i]=i*theStdcmap.red_mult;
	gdiv[i]=(i*theStdcmap.green_max)>>8;
	gmod[i]=((i*theStdcmap.green_max)>>2) & 63;
	gmul[i]=i*theStdcmap.green_mult;
	bdiv[i]=(i*theStdcmap.blue_max)>>8;
	bmod[i]=((i*theStdcmap.blue_max)>>2) & 63;
	bmul[i]=i*theStdcmap.blue_mult;
      }
    }
  } else {
   theColors = NULL;
  }
  /* alloc X Image store */
  /* bitmap = (unsigned char *) malloc( theImage->bytes_per_line *
  height); */
  /* bitmap = (unsigned char *) malloc( ((bits_per_RGB *width+32)/8) *
  height ); */

  /* bitmap = (unsigned char *) malloc( ((theDepth*width+32)/8) * height); */
  /* ((theDepth+7)/8) * width * height ); */

  /* alloc 24bit image store */

  rowlen = 2 * ((width * 8 + 15) / 16);
  planelen = rowlen * height;

  /* if (((theColors == NULL)&& !truecolor)
     ||(bitmap == NULL)) {
    fprintf(stderr,"ERROR: Can not allocate the buffers..\n");
    exit(1);
  }*/

  /* this avoids from crashing on minor error conditions. otherwise it may
     force core crashs... */
  XSetErrorHandler(error_handler);

  theWindow = openWindow(0,0,width,height,0,&theGC,s);
  initColor(width,height);
  XFlush(theDisplay);
}

void DeinitDisplay(void)
{
  if(!(theDisplay)) return; /* there's nothing to close */
  XFlush(theDisplay);
#ifdef MIT_SHM
  if (Shm) {
    XShmDetach(theDisplay, &shminfo);
    XDestroyImage(theImage);
    shmdt(shminfo.shmaddr);
    shmctl(shminfo.shmid,IPC_RMID,0);
  }
#endif
  XCloseDisplay(theDisplay);
  theDisplay = NULL;
  if(theColors) free(theColors);
  /* free(bitmap); */
}

void endpreview ()
{
  DeinitDisplay();
  fprintf(stderr,"Preview canceled. Stop.\n");
  exit(0);
}

void refreshWindow(Window theExposedWindow)
{
  if(!visible) return;  /* why do it, nobody sees me. */
  switch(theEvent.type) {
   case Expose : /* only refresh exposed region */
#ifdef MIT_SHM
    XShmPutImage(theDisplay,theWindow,theGC,theImage,
		theEvent.xexpose.x,theEvent.xexpose.y,
		theEvent.xexpose.x,theEvent.xexpose.y,
		theEvent.xexpose.width,theEvent.xexpose.height,0);
#else
    XPutImage(theDisplay,theWindow,theGC,theImage,
              theEvent.xexpose.x,theEvent.xexpose.y,
              theEvent.xexpose.x,theEvent.xexpose.y,
              theEvent.xexpose.width,theEvent.xexpose.height);
#endif
    break;
   default :     /* refresh all */
    XPutImage(theDisplay,theWindow,theGC,theImage,0,0,0,0,
              theImage->width,theImage->height);
  }
  XFlush(theDisplay);
  last_flush=time(NULL);
}


struct Region { int x,y,c; struct Region* n; };
typedef struct Region region;

#define rtsize 1024	/* Size of the cache. Reduce to have more frequent
			   screen updates */
#define maxlines 4096	/* Maximum number of rows. */
/* Sorry, but the cache scheme used here needs a max. number of lines to be 
   defined. It cannot easily made dynamic without losing speed. This
   number does not mean I waste much memory, since it's the size of two
   tables of pointers, so 4096 means 4096*pointer size (4 on 32 bit
   systems) or 16k for each table. Compared to the size of the whole
   bitmap this is "peanuts". */

region rt[rtsize];
int rtfree;
region *firstfree, *yt[maxlines],*yl[maxlines];

/* The cache works this way : I have a table with a linked list of
   regions for each display line (a region is a range of pixels on a
   line, given by the start colunn (x) and the width of the region (c)).
   The table is for speed, a single linked list for all lines would
   work, too, but it would be much slower, since for each pixel you
   would have to search through a potentially big list. This code
   minimizes the length of the list to search with the table lookup -
   you only have to look through the list of pixels for the current
   line. Looking through a list is tedious in any case, that's why I
   have regions. I don't have a list element for each pixel, only for
   regions. If you watch Polyray's screen output, the common case is
   several pixels in the same row on adjacent columns. This code makes
   use of the "adjacent" property by keeping not only a table with the
   start of the list for each row, but also a table with a pointer to
   the most recently used region for each row. So if the adjacent pixel
   is set, this code does not have to search the list, it can directly
   modify the right element. There are several cases to watch, for
   example joining two regions if the pixel between them is set and so
   on. This makes the code look more complicated than it is.

   Another speed gain has been achieved by not allocating the list
   elements with malloc() or calloc(). Instead, I have a global table
   with regions, and the allocation and deallocation of an element takes
   a constant time, no matter how much of the "dynamic store" is in
   use.  */

void InitRegionTable()
{
  int i;
  for (i=0;i<rtsize;i++) {
    rt[i].x=-1; rt[i].y=-1; rt[i].c=1; rt[i].n=&rt[i+1];
  }
  rt[rtsize-1].n=0;
  firstfree=&rt[0];
  rtfree=rtsize;
  for (i=0;i<maxlines;i++) { yt[i]=0; yl[i]=0; }
}

void FlushRegionTable()
{
  int i;
  region* r;

/* This redraws the screen roughly in the same order the PutPixels were
   issued. Looks a little less cached : */
  for (i=0;i<rtsize;i++) {
    r=&rt[i];
    if (r->x!=-1)
      XPutImage(theDisplay,theWindow,theGC,theImage,r->x,r->y,r->x,r->y,
		r->c,1);
  }
  /*   This redraws the screen from the first to the last row :
  for (i=0;i<maxlines;i++) {
    r=yt[i];
    while (r) {
      XPutImage(theDisplay,theWindow,theGC,theImage,r->x,i,r->x,i,
		r->c,1);
      r=r->n;
    }
  } */
  InitRegionTable();
}

region *ralloc()
{
  region *r;
  /* if (rtfree==0) FlushRegionTable(); */
  rtfree--;
  r=firstfree;
  firstfree=(rtfree?firstfree->n:0);
  assert(r->x==-1);
  return r;
}

void rfree(region* r)
{
  rtfree++;
  assert(r->x!=-1);
  r->x=-1;
  r->n=firstfree;
  firstfree=r;
}

void CachedPutPixel(int x,int y)
{
  region *r,*r1;

  if (rtfree==0) FlushRegionTable();

  if (yl[y]) {
    if (yl[y]->x+yl[y]->c==x) {
#ifdef DEBUG
      printf("Cache hit on line %d\n",y);
#endif
      yl[y]->c++;
      if (yl[y]->n) { /* Look if we can join with next region */
	r=yl[y]->n;
	if (r->x==x+1) {
#ifdef DEBUG
	  printf("Join on line %d\n",y);
#endif
	  yl[y]->n=r->n;
	  yl[y]->c+=r->c;
	  rfree(r);
	}
      }
      return;
    } else if (yl[y]->x<=x && yl[y]->x+yl[y]->c>x) {
#ifdef DEBUG
      printf("100%% Cache hit on line %d\n",y);
#endif
      return;
    }
  }
  /* OK, cache miss. Now check the rest of the regions on that line */
  if (yt[y]) {
    r=yt[y];
    if (r->x>x) {	/* Before the first region */
      if (r->x==x+1) {  /* ... directly ? */
	r->x--; r->c++; /* Yep. Join */
      } else {		/* Nope. Insert region */
	r=ralloc(); r->x=x; r->y=y; r->c=1; r->n=yt[y]; yt[y]=r;
      }
      yl[y]=r;
      return;
    }
    for (;;) {
      if (r->x<=x && r->x+r->c>x) {
#ifdef DEBUG
	printf("100%% 2nd level cache hit on line %d\n",y);
#endif
	yl[y]=r;
	return;
      }
      if (r->x+r->c==x) {
#ifdef DEBUG
	printf("2nd level cache hit on line %d\n",y);
#endif
	r->c++;
	if (r->n && r->n->x==r->x+r->c) {
#ifdef DEBUG
	  printf("2nd level cache join on line %d\n",y);
#endif
	  yl[y]=r; r=r->n; yl[y]->n=r->n; yl[y]->c+=r->c; rfree(r);
	  return;
	}
	yl[y]=r;
	return;
      }
      if (r->n) {
	if (r->n->x==x+1) {	/* Join with next region */
	  r=r->n;
	  r->x--; r->c++;
	  yl[y]=r;
	  return;
	}
	if (r->n->x>=x) { r=r->n; continue; }
#ifdef DEBUG
	printf("Insert region in line %d\n",y);
#endif
	r1=ralloc(); r1->n=r->n; r->n=r1;
	r1->x=x; r1->y=y; r1->c=1; yl[y]=r1;
	return;
      } else {
#ifdef DEBUG
	printf("Append region in line %d\n",y);
#endif
	r->n=ralloc(); r=r->n;
	r->x=x; r->y=y; r->c=1; r->n=0; yl[y]=r;
	return;
      }
    };
  } else {
#ifdef DEBUG
    printf("Create first region on line %d\n",y);
#endif
    r=ralloc();
    r->x=x; r->y=y; r->c=1; r->n=0; yt[y]=r; yl[y]=r;
  }
}


int dm[] = {     /* Dither Matrix */
       0,  48,  12,  60,   3,  51,  15,  63,
      32,  16,  44,  28,  35,  19,  47,  31,
       8,  56,   4,  52,  11,  59,   7,  55,
      40,  24,  36,  20,  43,  27,  39,  23,
       2,  50,  14,  62,   1,  49,  13,  61,
      34,  18,  46,  30,  33,  17,  45,  29,
      10,  58,   6,  54,   9,  57,   5,  53,
      42,  26,  38,  22,  41,  25,  37,  21
};

int rmax=0,gmax=0,bmax=0;


int PollEvents() {
  int numEvents;
  int i=1;

  numEvents = XPending(theDisplay);     /* Hmm, polling. but theres no way */
  if (numEvents > 0) {                  /* now deal with the events.. */
    XNextEvent(theDisplay,&theEvent);

    switch (theEvent.type) {
    case Expose:
      visible = True;                   /* yep. expose implies this */
      refreshWindow(theEvent.xany.window);
      break;

    case VisibilityNotify:
      if(theEvent.xvisibility.state==VisibilityFullyObscured) visible=False;
      break;

    case ClientMessage:
      /* we should also check for message_type beeing XA_PROTOCOLS */
      /* we act on WM_DELETE_WINDOW even if no +x was given */
      if(theEvent.xclient.data.l[0] == wm_delw) {endpreview(); return; };
      break;

    case KeyPress:
    case ButtonPress:
      i=0;
      break;
    }
  }
  return i;
}


int pixels=0;

void X11_PutPixel(int x,int y,byte R,byte G,byte B)
{
  unsigned long         color,tmp;
  unsigned char         *rgbpos;
  time_t		t;

  if(!(theDisplay)) return;
  PollEvents();

  color=theStdcmap.base_pixel;
  tmp=rdiv[R];
  if (tmp<theStdcmap.red_max & rmod[R]>dm[(x&7)+((y&7)<<3)]) tmp++;
  color+=rmul[tmp];
  tmp=gdiv[G];
  if (tmp<theStdcmap.green_max & gmod[G]>dm[(x&7)+((y&7)<<3)]) tmp++;
  color+=gmul[tmp];
  tmp=bdiv[B];
  if (tmp<theStdcmap.blue_max & bmod[B]>dm[(x&7)+((y&7)<<3)]) tmp++;
  color+=bmul[tmp];

/* color = rgb_to_pixel(R,G,B);*/

  /* put the pixel to the bitmap */
  /* XPutPixel(theImage, x, y, color ); */
  rgbpos=theImage->data+ytab[y]+(x<<xshift);
  switch (theImage->bits_per_pixel) {
    case 32: *((unsigned long int *) rgbpos)=color; break;
    case 15: case 16: *((unsigned short int *) rgbpos)=color & 65535; break;
    case 8: *rgbpos = color; break;
  }

  CachedPutPixel(x,y);
  if (++pixels==100) {
    t=time(NULL);
    if (t-last_flush>1) {
      FlushRegionTable();
      XFlush(theDisplay);
      last_flush=t;
    }
    pixels=0;
  }
}

void FlushDisplay()
{
  if (theDisplay) {
    FlushRegionTable();
    XFlush(theDisplay);
    last_flush=time(NULL);
  }
}

void SpecialStatus(const char *s)
{
  if (theDisplay) {
    XSetIconName(theDisplay, theWindow, s);
    XFlush(theDisplay);
  }
}
