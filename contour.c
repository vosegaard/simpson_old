#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tcl.h>
#include <math.h>
#include "config.h"
#include "iodata.h"
#include <stdarg.h>
#ifdef WITH_GD
#include "gd.h"
#include "gdfontmb.h"
#endif

extern FD **fd;

#define DECLARE_CMD(func) \
int func _ANSI_ARGS_((ClientData clientData, \
		   Tcl_Interp *interp, int argc, char *CONST argv[]))

#define DEFINE_CMD(func) \
int func(clientData, interp, argc, argv) \
    ClientData clientData; \
    Tcl_Interp *interp;	\
    int argc; \
    char *CONST argv[];

#define ADD_CMD(cmdName, cmdProc) \
Tcl_CreateCommand(interp, #cmdName, cmdProc, NULL, NULL)

#define OPEN   0
#define FILLED 1

int* point2(int x1, int y1, int x2, int y2)
{
  static int p[4];
  p[0] = x1; p[1] = y1;
  p[2] = x2; p[3] = y2;
  return p;
}
int* point3(int x1, int y1, int x2, int y2, int x3, int y3)
{
  static int p[6];
  p[0] = x1; p[1] = y1;
  p[2] = x2; p[3] = y2;
  p[4] = x3; p[5] = y3;
  return p;
}
int* point4(int x1, int y1, int x2, int y2, int x3, int y3,
            int x4, int y4)
{
  static int p[8];
  p[0] = x1; p[1] = y1;
  p[2] = x2; p[3] = y2;
  p[4] = x3; p[5] = y3;
  p[6] = x4; p[7] = y4;
  return p;
}
int* point5(int x1, int y1, int x2, int y2, int x3, int y3,
            int x4, int y4, int x5, int y5)
{
  static int p[10];
  p[0] = x1; p[1] = y1;
  p[2] = x2; p[3] = y2;
  p[4] = x3; p[5] = y3;
  p[6] = x4; p[7] = y4;
  p[8] = x5; p[9] = y5;
  return p;
}

#ifdef WITH_GD
struct tmp {
  int bytes;
  int alloc;
  char *content;
};

static int dump(void *context, const char *buffer, int len)
{
  struct tmp *c = (struct tmp*)context;
  if (c->bytes+len > c->alloc) {
    fprintf(stderr,"Allocation error in dump\n");
    exit(1);
  }
  memcpy(c->content+c->bytes, buffer, len);
  c->bytes += len;
  return 1;
}
#endif

int XFIG_AllocateColor(void *ptr, int r, int g, int b)
{
  static int color = 31;
  FILE *file = (FILE*)ptr;
  color++;
  if (file != NULL)
    fprintf(file, "0 %d #%02x%02x%02x\n", color, r, g, b);
  return color;
}

int GD_AllocateColor(void *ptr, int r, int g, int b)
{
#ifdef WITH_GD
  return gdImageColorAllocate((gdImagePtr)ptr,r,g,b);
#else
  return 0;
#endif
}


int XFIG_DrawPolygon(void *ptr, int color, int depth, 
  int npt, int *p, int fill)
{
  FILE *file = (FILE*)ptr;
  int i;

  fprintf(file, "2 3 0 1 %d %d %d 0 %d 0.000 1 0 -1 0 0 %d\n\t",
    color,color,depth,fill==FILLED ? 20:-1, npt+1);
  for (i=0; i < npt*2; i+=2) {
    fprintf(file,"%d %d ", p[i], p[i+1]);
    if (i % 6 == 4) fprintf(file,"\n\t");
  }
  fprintf(file,"%d %d\n", p[0],p[1]);
  return 1;
}

int GD_DrawPolygon(void *ptr, int color, int depth, 
  int npt, int *p, int fill)
{
#ifdef WITH_GD
  gdImagePtr img = (gdImagePtr)ptr;

  if (fill == FILLED) {
    gdImageFilledPolygon(img,(gdPointPtr)p,npt,color);
  } else {
    gdImagePolygon(img,(gdPointPtr)p,npt,color);
  }
#endif
  return 1;
}

int GD_DrawLines(void *ptr, int color, int depth, 
  int npt, int *p)
{
#ifdef WITH_GD
  int i;
  gdImagePtr img = (gdImagePtr)ptr;
  for (i=0; i<(npt-1)*2;i+=2) {
    gdImageLine(img, p[i], p[i+1],p[i+2],p[i+3], color);
  }
#endif
  return 1;
}

int XFIG_DrawLines(void *ptr, int color, int depth,
  int npt, int *p)
{
  FILE *file = (FILE*)ptr;
  int i;

  fprintf(file, "2 1 0 1 %d %d %d 0 -1 0.000 0 0 -1 0 0 %d\n\t",
    color,color,depth,npt);
  for (i=0; i < npt*2; i+=2) {
    fprintf(file,"%d %d ", p[i], p[i+1]);
    if (i%6 == 4 && i != npt*2-2) fprintf(file,"\n\t");
  }
  fprintf(file,"\n");
  return 1;
}








int (*AllocateColor)(void *,int,int,int);
int (*DrawPolygon)(void*,int,int,int,int*,int);
int (*DrawLines)(void*,int,int,int,int*);



typedef struct _Point {
  double x, y, f;
} Point;

Point* min(Point* p1, Point* p2, Point* p3)
{
  if (p1->f <= p2->f && p1->f <= p3->f) return p1;
  if (p2->f <= p1->f && p2->f <= p3->f) return p2;
  return p3;
}
Point* max(Point* p1, Point* p2, Point* p3)
{
  if (p3->f >= p1->f && p3->f >= p2->f) return p3;
  if (p1->f >= p2->f && p1->f >= p3->f) return p1;
  return p2;
}
Point* mid(Point* p1, Point* p2, Point* p3)
{
  Point *a = max(p1,p2,p3);
  Point *b = min(p1,p2,p3);
  if (p1 != a && p1 != b) return p1;
  if (p2 != a && p2 != b) return p2;
  return p3;
}


DECLARE_CMD(FContour);

int AddContourCommand(Tcl_Interp *interp)
{
/* Syntax
  ADD_CMD(Tcl_Command_name, Subroutine name);
*/
  ADD_CMD(fcontour, FContour);
  return TCL_OK;
}


#define PLOT_POS (1<<0)
#define PLOT_NEG (1<<1)

#define X(ii) (int)((ii)+0.5)
#define Y(ii) (int)((ii)+0.5)
   

#define INT(i,j) (re ? fdp->data[((j)-1+jfrom)*fdp->np+(i)+ifrom].re:fdp->data[((j)-1+jfrom)*fdp->np+(i)+ifrom].im)

#define XFIG 1
#define GD   2
#define NOXFIGINIT   0
#define ONLYXFIGINIT 1
#define XFIGALL      2

DEFINE_CMD(FContour)
{
  int i, j, autolevels=1, nlevels, skiplevels, pos, itmp;
  int desc, w, h, re, fill, l, ac, yup=0, verbose = 0;
  int ifrom, ito, jfrom, jto, nx, ny, ni, np,add=0;
  int xreverse=0, black,white, wset=0, hset=0,output=XFIGALL;
  double lowcolor[3], highcolor[3], xmin, xmax, ymin, ymax;
  double *levels, x1, y1, x2, y2, cl, r1, r2, fmax;
  double *xpos, *ypos, **z, iratio, jratio, nlowcolor[3], nhighcolor[3];
  char outputfile[256], **av;
  FILE *file;
  FD *fdp;
  Point *p, **pp;
  int tri, c, **t, type = XFIG;
  int *colors;
  void *ptr;
#ifdef WITH_GD
  gdImagePtr img;
#endif
    
  pp = (Point**)malloc(sizeof(Point*)*6);
  p = (Point*)malloc(sizeof(Point)*5);
  t = (int**)malloc(sizeof(int*)*4);
  t[0] = (int*)malloc(sizeof(int)*13);
  for (i=1; i<4; i++) {
    t[i] = t[i-1]+3;
  }
  t[0][0] = 0; t[0][1] = 1; t[0][2] = 4;
  t[1][0] = 1; t[1][1] = 3; t[1][2] = 4;
  t[2][0] = 3; t[2][1] = 2; t[2][2] = 4;
  t[3][0] = 2; t[3][1] = 0; t[3][2] = 4;
    
  if (argc < 2) {
    fprintf(stderr,"Usage: %s desc ?options?\n\n", argv[0]);
    fprintf(stderr,"Options:\n"
    "         -im                 Takes the imaginary part of the data\n"
    "         -levels {list}      Specify a list of contour levels\n"
    "                             If not specified, autolevels will be used\n"
    "         -nlevels num        Number of autolevels\n"
    "         -skiplevels num     Skip lowest levels\n"
    "         -pos                Display only positive contours (Default)\n"
    "         -neg                Display only negative contours\n"
    "         -posneg             Display both positive and negative contours\n"
    "                             In this case the same levels will be used for\n"
    "                             positive and negative levels\n"
    "         -lowcolor {r g b}   Color of the lowest level. r, g, and b must\n"
    "                             be specified between 0 and 1 (Default 0 0 1)\n"
    "         -highcolor {r g b}  Color of the highest level (Default 1 0 0)\n"
    "         -lowcolor2 {r g b}  Color of the lowest level of negative contours in posneg mode(Default 1 0 0)\n"
    "         -highcolor2 {r g b} Color of the highest level of negative contours in posneg mode(Default 1 0 0)\n"
    "         -fill               Fill contours. This causes the contour lines\n"
    "                             to appear in black and the filled levels\n"
    "                             appear with the specified color\n"
    "         -xmin val           first point in the x dimension\n"
    "         -xmax val           last point in x dimension\n"
    "         -ymin val           first point in the y dimension\n"
    "         -ymax val           last point in y dimension\n"
    "         -ydown              y axis points downwards (default)\n"
    "         -yup                y axis points upwards\n"
    "         -width num          width of the plot window (5000)\n"
    "         -height num         height of the plot window (5000)\n"
    "         -type type          output type xfig/gif\n"
    "         -noxfiginit         don't output colors to xfig device\n"
    "         -onlyxfiginit       only output colors - no spectrum drawn\n"
    "         -xreverse           reverse x axis\n"
    "         -v                  verbose output\n");
    exit(0);
  }

  if (Tcl_GetInt(interp,argv[1],&desc) != TCL_OK) {
    fprintf(stderr, "Error: %s: specify spectrum descriptor\n",
      argv[1]);
    exit(0);
  }
  fdp = fd[desc];

  /* Setting default values */
  sprintf(outputfile, "%s/.xfig", getenv("HOME"));
  autolevels = 1;
  nlevels = 10;
  skiplevels = 2;
  pos = PLOT_POS;
  re = 1;
  fill = 0;
  w = 5000;
  h = 5000;
  lowcolor[1]  = lowcolor[2]  = 0; lowcolor[0]  = 1;
  highcolor[0] = highcolor[1] = 0; highcolor[2] = 1;
  nlowcolor[1]  = nlowcolor[2]  = 1; nlowcolor[0]  = 0;
  nhighcolor[0] = nhighcolor[1] = 1; nhighcolor[2] = 0;
  xmin = ymin = 1;
  xmax = fdp->np;
  ymax = fdp->ni;
  fmax = -1;
  np = fdp->np;
  ni = fdp->ni;

  for (i=2; i<argc; i++) {
    if (!strcmp(argv[i], "-im")) {
      re = 0;
    } else if (!strcmp(argv[i], "-pos")) {
      pos = PLOT_POS;
    } else if (!strcmp(argv[i], "-neg")) {
      pos = PLOT_NEG;
    } else if (!strcmp(argv[i], "-posneg")) {
      pos = PLOT_POS|PLOT_NEG;
    } else if (!strcmp(argv[i], "-fill")) {
      fill = 1;
    } else if (!strcmp(argv[i], "-xreverse")) {
      xreverse = 1;
    } else if (!strcmp(argv[i], "-ydown")) {
      yup = 0;
    } else if (!strcmp(argv[i], "-noxfiginit")) {
      output = NOXFIGINIT;
    } else if (!strcmp(argv[i], "-onlyxfiginit")) {
      output = ONLYXFIGINIT;
    } else if (!strcmp(argv[i], "-yup")) {
      yup = 1;
    } else if (!strcmp(argv[i], "-v")) {
      verbose = 1;
    } else if (!strcmp(argv[i], "-levels") && i < argc-1) {
      autolevels = 0;
      i++;
      if (Tcl_SplitList(interp,argv[i],
                        &ac, &av) != TCL_OK) {
        fprintf(stderr,
	  "Warning: %s: unable to read list\n",
	  argv[0]);
	autolevels=1;
      } else {
        nlevels = ac;
        levels = (double*)malloc(nlevels*sizeof(double));
	for (j=0; j<nlevels; j++) {
	  if (Tcl_GetDouble(interp,av[j],&levels[j]) != TCL_OK) {
	    fprintf(stderr,
	      "Warning: %s: element %d of level list is not a double\n",
	      argv[0], j);
	    autolevels=1;
	    free((char*)levels);
	  }
	}
	Tcl_Free((char*)av);
      }
    } else if (!strcmp(argv[i], "-nlevels") && i < argc-1) {
      if (Tcl_GetInt(interp,argv[i+1],&itmp) != TCL_OK) {
        fprintf(stderr, 
	  "Warning: %s: unable to read integer '%s' for argument -nlevels",
	  argv[0], argv[i+1]);
      }
      i++;
      nlevels = itmp;
    } else if (!strcmp(argv[i], "-skiplevels") && i < argc-1) {
      if (Tcl_GetInt(interp,argv[i+1],&itmp) != TCL_OK) {
        fprintf(stderr, 
	  "Warning: %s: unable to read integer '%s' for argument -skiplevels",
	  argv[0], argv[i+1]);
      }
      i++;
      skiplevels = itmp;
    } else if (!strcmp(argv[i], "-lowcolor") && i < argc-1) {
      if (Tcl_SplitList(interp,argv[i+1],&ac,&av) != TCL_OK) {
        fprintf(stderr,
	  "Warning: %s: unable to read list '%s' for argument -lowcolor\n",
	  argv[0], argv[i+1]);
      }
      if (ac == 3) {
        for (j=0; j<3; j++) {
	  Tcl_GetDouble(interp,av[j],&lowcolor[j]);
	}
      } else {
        fprintf(stderr,
	  "Warning: %s: list '%s' for argument -lowcolor does not have three elements {r g b}\n",
	  argv[0], argv[i+1]);
      }
      Tcl_Free((char*)av);
      i++;
    } else if (!strcmp(argv[i], "-highcolor") && i < argc-1) {
      if (Tcl_SplitList(interp,argv[i+1],&ac,&av) != TCL_OK) {
        fprintf(stderr,
	  "Warning: %s: unable to read list '%s' for argument -highcolor\n",
	  argv[0], argv[i+1]);
      }
      if (ac == 3) {
        for (j=0; j<3; j++) {
	  Tcl_GetDouble(interp,av[j],&highcolor[j]);
	}
      } else {
        fprintf(stderr,
	  "Warning: %s: list '%s' for argument -highcolor does not have three elements {r g b}\n",
	  argv[0], argv[i+1]);
      }
      Tcl_Free((char*)av);
      i++;
    } else if (!strcmp(argv[i], "-lowcolor2") && i < argc-1) {
      if (Tcl_SplitList(interp,argv[i+1],&ac,&av) != TCL_OK) {
        fprintf(stderr,
	  "Warning: %s: unable to read list '%s' for argument -lowcolor2\n",
	  argv[0], argv[i+1]);
      }
      if (ac == 3) {
        for (j=0; j<3; j++) {
	  Tcl_GetDouble(interp,av[j],&nlowcolor[j]);
	}
      } else {
        fprintf(stderr,
	  "Warning: %s: list '%s' for argument -lowcolor2 does not have three elements {r g b}\n",
	  argv[0], argv[i+1]);
      }
      Tcl_Free((char*)av);
      i++;
    } else if (!strcmp(argv[i], "-highcolor2") && i < argc-1) {
      if (Tcl_SplitList(interp,argv[i+1],&ac,&av) != TCL_OK) {
        fprintf(stderr,
	  "Warning: %s: unable to read list '%s' for argument -highcolor2\n",
	  argv[0], argv[i+1]);
      }
      if (ac == 3) {
        for (j=0; j<3; j++) {
	  Tcl_GetDouble(interp,av[j],&nhighcolor[j]);
	}
      } else {
        fprintf(stderr,
	  "Warning: %s: list '%s' for argument -highcolor does not have three elements {r g b}\n",
	  argv[0], argv[i+1]);
      }
      Tcl_Free((char*)av);
      i++;
    } else if (!strcmp(argv[i], "-xmin") && i < argc-1) {
      if (Tcl_GetDouble(interp,argv[i+1],&xmin) != TCL_OK) {
        fprintf(stderr, 
	  "Warning: %s: unable to read double '%s' for argument -xmin",
	  argv[0], argv[i+1]);
      }
      i++;
    } else if (!strcmp(argv[i], "-xmax") && i < argc-1) {
      if (Tcl_GetDouble(interp,argv[i+1],&xmax) != TCL_OK) {
        fprintf(stderr, 
	  "Warning: %s: unable to read double '%s' for argument -xmax",
	  argv[0], argv[i+1]);
      }
      i++;
    } else if (!strcmp(argv[i], "-ymin") && i < argc-1) {
      if (Tcl_GetDouble(interp,argv[i+1],&ymin) != TCL_OK) {
        fprintf(stderr, 
	  "Warning: %s: unable to read integer '%s' for argument -ymin",
	  argv[0], argv[i+1]);
      }
      i++;
    } else if (!strcmp(argv[i], "-ymax") && i < argc-1) {
      if (Tcl_GetDouble(interp,argv[i+1],&ymax) != TCL_OK) {
        fprintf(stderr, 
	  "Warning: %s: unable to read integer '%s' for argument -ymax",
	  argv[0], argv[i+1]);
      }
      i++;
    } else if (!strcmp(argv[i], "-type") && i < argc-1) {
      if (!strcasecmp(argv[i+1],"gif")) {
        type = GD;
        if (!wset) w /= 10;
	if (!hset) h /= 10;
      }
      i++;
    } else if (!strcmp(argv[i], "-width") && i < argc-1) {
      if (Tcl_GetInt(interp,argv[i+1],&itmp) != TCL_OK) {
        fprintf(stderr, 
	  "Warning: %s: unable to read integer '%s' for argument -width",
	  argv[0], argv[i+1]);
      }
      i++;
      w = itmp;
      wset = 1;
    } else if (!strcmp(argv[i], "-height") && i < argc-1) {
      if (Tcl_GetInt(interp,argv[i+1],&itmp) != TCL_OK) {
        fprintf(stderr, 
	  "Warning: %s: unable to read integer '%s' for argument -height",
	  argv[0], argv[i+1]);
      }
      i++;
      h = itmp;
      hset = 1;
    } else if (!strcmp(argv[i], "-output") && i < argc-1) {
    } else {
      fprintf(stderr, "Warning: fcontour: unknown option '%s'\n", argv[i]);
    }
  }

/* Now all input parameters are read and we're ready to proceed */

  if (xmin >= xmax) {
    fprintf(stderr, "Error: %s: xmin must be smaller than xmax\n", argv[0]);
    return TCL_ERROR;
  }
  if (ymin >= ymax) {
    fprintf(stderr, "Error: %s: ymin must be smaller than ymax\n", argv[0]);
    return TCL_ERROR;
  }
  ifrom = xmin < 1 ? 1:(int)xmin;
  if (ifrom > np) ifrom = np;
  ito   = xmax > (double)np ? np:(int)ceil(xmax);
  if (ito < 1) ito = 1;
  jfrom = ymin < 1 ? 1:(int)ymin;
  if (jfrom > ni) jfrom = ni;
  jto   = ymax > (double)ni ? ni:(int)ceil(ymax);
  if (jto < 1) jto = 1;
  nx = ito - ifrom + 1;
  ny = jto - jfrom + 1;
  xpos = (double*)malloc(sizeof(double)*nx);
  ypos = (double*)malloc(sizeof(double)*ny);
  z = (double**)malloc(sizeof(double*)*ny);
  z[0] = (double*)malloc(sizeof(double)*ny*nx);
  for (j=0; j<ny; j++) {
    if (j > 0) z[j] = z[j-1]+nx; /* setup matrix */
    jratio = 1;
    ypos[j] = (double)(j+jfrom);
    if (j == 0 && ymin > 1) {
      jratio = 1.0-(ymin-(double)jfrom);
      ypos[j] = ymin;
    } else if (j == ny-1 && ymax < (double)ni) {
      jratio = (double)j+jfrom-ymax-1.0;
      ypos[j] = ymax;
    }
    if (yup == 1) {
      ypos[j] = (1.0-(ypos[j]-ymin)/(ymax-ymin))*(double)h;
    } else {
      ypos[j] = (ypos[j]-ymin)/(ymax-ymin)*(double)h;
    }
    for (i=0; i<nx; i++) {
      iratio = 1;
      xpos[i] = (double)i+ifrom;
      if (i == 0 && xmin > 1) {
        iratio = 1.0-(xmin-(double)ifrom);
        xpos[i] = xmin;
      } else if (i == nx-1 && xmax < (double)np) {
        iratio = (double)i+ifrom-xmax-1.0;
        xpos[i] = xmax;
      }
      if (xreverse == 1) {
        xpos[i] = (xpos[i]-xmin)/(xmax-xmin)*(double)w;
      } else {
        xpos[i] = (1.0-(xpos[i]-xmin)/(xmax-xmin))*(double)w;
      }
      if (fabs(iratio) < 1e-6) {
	if (fabs(jratio) < 1e-6) {
	  z[j][i] = INT(i,j);
	} else if (jratio < 0) {
	  z[j][i] = -jratio*INT(i,j) + (jratio+1.0)*INT(i,j-1);
	} else {
	  z[j][i] = jratio*INT(i,j) + (1.0-jratio)*INT(i,j+1);
	}
      } else if (iratio < 0) {
	if (fabs(jratio) < 1e-6) {
	  z[j][i] = -iratio*INT(i,j) + (iratio+1.0)*INT(i-1,j);
	} else if (jratio < 0) {
	  z[j][i] = -iratio*(-jratio*INT(i,j) + (jratio+1.0)*INT(i,j-1)) +
	    (iratio+1.0)*(-jratio*INT(i-1,j) + (jratio+1.0)*INT(i-1,j-1));
	} else {
	  z[j][i] = -iratio*(jratio*INT(i,j) + (1.0-jratio)*INT(i,j+1)) +
	    (iratio+1.0)*(jratio*INT(i-1,j) + (1.0-jratio)*INT(i-1,j+1));
	}
      } else {
	if (fabs(jratio) < 1e-6) {
	  z[j][i] = iratio*INT(i,j) + (1.0-iratio)*INT(i-1,j);
	} else if (jratio < 0) {
	  z[j][i] = iratio*(-jratio*INT(i,j) + (jratio+1.0)*INT(i,j-1)) +
	    (1.0-iratio)*(-jratio*INT(i+1,j) + (jratio+1.0)*INT(i+1,j-1));
	} else {
	  z[j][i] = iratio*(jratio*INT(i,j) + (1.0-jratio)*INT(i,j+1)) +
	    (1.0-iratio)*(jratio*INT(i+1,j) + (1.0-jratio)*INT(i+1,j+1));
	}
      }
      if (pos & PLOT_NEG) {
        if (-z[j][i] > fmax) fmax = -z[j][i];
      } 
      if (pos & PLOT_POS) {
        if (z[j][i] > fmax) fmax = z[j][i];
      }
    }
  }
  if (pos & PLOT_POS && pos & PLOT_POS) {
    colors = (int*)malloc(2*nlevels*sizeof(int));
  } else {
    colors = (int*)malloc(nlevels*sizeof(int));
  }
  if (type == GD) {
#ifdef WITH_GD
    img = gdImageCreate(w,h);
    ptr = (void*)img;
    AllocateColor = GD_AllocateColor;
    DrawPolygon   = GD_DrawPolygon;
    DrawLines     = GD_DrawLines;
    white = AllocateColor(ptr,255,255,255);
    gdImageColorTransparent(img,white);
    black = AllocateColor(ptr,0,0,0);
#else
    fprintf(stderr, "Error: %s: You are trying to use GD but it is not\n"
	    "compiled into this version of SIMPSON\n", argv[0]);
    return TCL_ERROR;
#endif
  } else {
    if ((file = fopen(outputfile,(output==NOXFIGINIT) ? "a":"w")) == NULL) {
      fprintf(stderr, "Warning: %s: Unable to write file '%s'. Using stdout\n",
      argv[0], outputfile);
      file = stdout;
    }
    if (verbose == 1) printf("Using file: %s\n", outputfile);
    if (output != NOXFIGINIT) 
      fprintf(file, "#FIG 3.2\nPortrait\nCenter\nMetric\nA4\n100.00\nSingle\n-2\n1200 2\n");
    ptr = (void*)file;
    AllocateColor = XFIG_AllocateColor;
    DrawPolygon   = XFIG_DrawPolygon;
    DrawLines     = XFIG_DrawLines;
    black = AllocateColor((output == NOXFIGINIT) ? NULL:ptr,0,0,0);
  }
  
  if (autolevels == 1) {
    levels = (double*)malloc(nlevels*sizeof(double));
    for (i=0; i<nlevels; i++) {
      levels[i] = fmax/(nlevels+skiplevels+1)*(i+skiplevels+1);
      if (verbose)
	printf("level %d = %f (%.0f %%)\n", i+1,levels[i],100*levels[i]/fmax);
    }
  }

  for (i=0; i<nlevels; i++) {
    cl = (double)i/(double)(nlevels-1);
    colors[i] = AllocateColor(
        (output == NOXFIGINIT && type == XFIG) ? NULL:ptr,
        (int)(255.0*(cl*highcolor[0] + (1.0-cl)*lowcolor[0])),
        (int)(255.0*(cl*highcolor[1] + (1.0-cl)*lowcolor[1])),
        (int)(255.0*(cl*highcolor[2] + (1.0-cl)*lowcolor[2])));
  }
  if (pos & PLOT_POS && pos & PLOT_NEG) {
    for (i=0; i<nlevels; i++) {
      cl = (double)i/(double)(nlevels-1);
      colors[i+nlevels] = AllocateColor(
          (output == NOXFIGINIT && type == XFIG) ? NULL:ptr,
          (int)(255.0*(cl*nhighcolor[0] + (1.0-cl)*nlowcolor[0])),
          (int)(255.0*(cl*nhighcolor[1] + (1.0-cl)*nlowcolor[1])),
          (int)(255.0*(cl*nhighcolor[2] + (1.0-cl)*nlowcolor[2])));
    }
  }
  
  if (type == XFIG && output == ONLYXFIGINIT) goto DONE;
  
  if (type == GD) {
    DrawPolygon(ptr,black,50,4,point4(0,0,w-1,0,w-1,h-1,0,h-1),OPEN);
  } else {
    DrawPolygon(ptr,black,50,4,point4(0,0,w,0,w,h,0,h),OPEN);
  }
  if (pos & PLOT_POS) {
    for (j = 0; j < ny-1; j++) {
      for (i = 0; i < nx-1; i++) {
        if (z[j][i] < levels[0] && z[j+1][i] < levels[0] &&
  	    z[j][i+1] < levels[0] && z[j+1][i+1] < levels[0]) continue;
        p[0].x = xpos[i];     p[0].y = ypos[j];     p[0].f = z[j][i];
        p[1].x = xpos[i+1];   p[1].y = ypos[j];     p[1].f = z[j][i+1];
        p[2].x = xpos[i];     p[2].y = ypos[j+1];   p[2].f = z[j+1][i];
        p[3].x = xpos[i+1];   p[3].y = ypos[j+1];   p[3].f = z[j+1][i+1];
        p[4].x = (xpos[i]+xpos[i+1])/2.0; p[4].y = (ypos[j]+ypos[j+1])/2.0;
        p[4].f = (z[j][i] + z[j+1][i] + z[j][i+1] + z[j+1][i+1])/4.0;
        for (tri = 0; tri < 4; tri++) {
          for (c = 0; c < 3; c++) {
 	    pp[c+3] = &p[t[tri][c]];
	  }

          pp[0] = min(pp[3],pp[4],pp[5]);
          pp[1] = mid(pp[3],pp[4],pp[5]);
          pp[2] = max(pp[3],pp[4],pp[5]);

/* Plot contours */

          if (fill) {
            c = -1;
  	    for (l=0; l<nlevels; l++) {
              if (levels[l] - pp[0]->f < 1e-5) c = l;
	    }
	    if (c > -1) {
	      DrawPolygon(ptr,colors[c],50+2*nlevels-2*c+1,
	        3, point3(X(pp[0]->x), Y(pp[0]->y), 
	        X(pp[1]->x), Y(pp[1]->y),
	        X(pp[2]->x), Y(pp[2]->y)),FILLED);
            }
          }
	
 	  for (l=0; l<nlevels; l++) {
	    cl = levels[l];
	    if (cl < pp[0]->f || cl > pp[2]->f) {
	      continue;
	    }

 	    if (fabs(cl-pp[0]->f) < 1e-5 && fabs(cl-pp[1]->f) < 1e-5) {
	      x1 = pp[0]->x;
	      y1 = pp[0]->y;
	      x2 = pp[1]->x;
	      y2 = pp[1]->y;
	    } else if (fabs(cl-pp[1]->f) < 1e-5 &&  fabs(cl-pp[2]->f) < 1e-5) {
	      x1 = pp[1]->x;
	      y1 = pp[1]->y;
	      x2 = pp[2]->x;
	      y2 = pp[2]->y;
	    } else {
	      r1 = (cl-pp[0]->f)/(pp[2]->f - pp[0]->f);
	      x1 = pp[2]->x*r1 + pp[0]->x*(1.0-r1);
	      y1 = pp[2]->y*r1 + pp[0]->y*(1.0-r1);
              if (cl < pp[1]->f) {
	        r2 = (cl-pp[0]->f)/(pp[1]->f - pp[0]->f);
	        x2 = pp[1]->x*r2 + pp[0]->x*(1.0-r2);
	        y2 = pp[1]->y*r2 + pp[0]->y*(1.0-r2);
	      } else {
	        r2 = (cl-pp[1]->f)/(pp[2]->f - pp[1]->f);
	        x2 = pp[2]->x*r2 + pp[1]->x*(1.0-r2);
	        y2 = pp[2]->y*r2 + pp[1]->y*(1.0-r2);
	      }
	    }
	    if (fill) {
	      if (cl < pp[1]->f) {
	        DrawPolygon(ptr,colors[l],50+2*nlevels-2*l+1,4,
	          point4(X(x1), Y(y1), X(x2), Y(y2), 
  	  	  X(pp[1]->x), Y(pp[1]->y),
		  X(pp[2]->x), Y(pp[2]->y)),FILLED);
	      } else {
                DrawPolygon(ptr,colors[l],50+2*nlevels-2*l+1,3,
	          point3(X(x1), Y(y1), X(x2), Y(y2), 
  		  X(pp[2]->x), Y(pp[2]->y)), FILLED);
	      }
  	      DrawLines(ptr,black,50+2*nlevels-2*l, 2,
	        point2(X(x1), Y(y1), X(x2), Y(y2)));
            } else {
  	      DrawLines(ptr,colors[l],50+2*nlevels-2*l, 2,
	        point2(X(x1), Y(y1), X(x2), Y(y2)));
	    }
 	  }
        }
      }
    }
  }
  if (pos & PLOT_NEG) {
    add = (pos & PLOT_POS) ? nlevels:0;
    for (j = 0; j < ny-1; j++) {
      for (i = 0; i < nx-1; i++) {
        if (-z[j][i] < levels[0] && -z[j+1][i] < levels[0] &&
  	    -z[j][i+1] < levels[0] && -z[j+1][i+1] < levels[0]) continue;
        p[0].x = xpos[i];     p[0].y = ypos[j];     p[0].f = -z[j][i];
        p[1].x = xpos[i+1];   p[1].y = ypos[j];     p[1].f = -z[j][i+1];
        p[2].x = xpos[i];     p[2].y = ypos[j+1];   p[2].f = -z[j+1][i];
        p[3].x = xpos[i+1];   p[3].y = ypos[j+1];   p[3].f = -z[j+1][i+1];
        p[4].x = (xpos[i]+xpos[i+1])/2.0; p[4].y = (ypos[j]+ypos[j+1])/2.0;
        p[4].f = -(z[j][i] + z[j+1][i] + z[j][i+1] + z[j+1][i+1])/4.0;
        for (tri = 0; tri < 4; tri++) {
          for (c = 0; c < 3; c++) {
 	    pp[c+3] = &p[t[tri][c]];
	  }

          pp[0] = min(pp[3],pp[4],pp[5]);
          pp[1] = mid(pp[3],pp[4],pp[5]);
          pp[2] = max(pp[3],pp[4],pp[5]);

/* Plot contours */

          if (fill) {
            c = -1;
  	    for (l=0; l<nlevels; l++) {
              if (levels[l] - pp[0]->f < 1e-5) c = l;
	    }
	    if (c > -1) {
	      DrawPolygon(ptr,colors[c+add],50+2*nlevels-2*(c+add)+1,
	        3, point3(X(pp[0]->x), Y(pp[0]->y), 
	        X(pp[1]->x), Y(pp[1]->y),
	        X(pp[2]->x), Y(pp[2]->y)),FILLED);
            }
          }
	
 	  for (l=0; l<nlevels; l++) {
	    cl = levels[l];
	    if (cl < pp[0]->f || cl > pp[2]->f) {
	      continue;
	    }

 	    if (fabs(cl-pp[0]->f) < 1e-5 && fabs(cl-pp[1]->f) < 1e-5) {
	      x1 = pp[0]->x;
	      y1 = pp[0]->y;
	      x2 = pp[1]->x;
	      y2 = pp[1]->y;
	    } else if (fabs(cl-pp[1]->f) < 1e-5 &&  fabs(cl-pp[2]->f) < 1e-5) {
	      x1 = pp[1]->x;
	      y1 = pp[1]->y;
	      x2 = pp[2]->x;
	      y2 = pp[2]->y;
	    } else {
	      r1 = (cl-pp[0]->f)/(pp[2]->f - pp[0]->f);
	      x1 = pp[2]->x*r1 + pp[0]->x*(1.0-r1);
	      y1 = pp[2]->y*r1 + pp[0]->y*(1.0-r1);
              if (cl < pp[1]->f) {
	        r2 = (cl-pp[0]->f)/(pp[1]->f - pp[0]->f);
	        x2 = pp[1]->x*r2 + pp[0]->x*(1.0-r2);
	        y2 = pp[1]->y*r2 + pp[0]->y*(1.0-r2);
	      } else {
	        r2 = (cl-pp[1]->f)/(pp[2]->f - pp[1]->f);
	        x2 = pp[2]->x*r2 + pp[1]->x*(1.0-r2);
	        y2 = pp[2]->y*r2 + pp[1]->y*(1.0-r2);
	      }
	    }
	    if (fill) {
	      if (cl < pp[1]->f) {
	        DrawPolygon(ptr,colors[l+add],50+2*nlevels-2*(l+add)+1,4,
	          point4(X(x1), Y(y1), X(x2), Y(y2), 
  	  	  X(pp[1]->x), Y(pp[1]->y),
		  X(pp[2]->x), Y(pp[2]->y)),FILLED);
	      } else {
                DrawPolygon(ptr,colors[l+add],50+2*nlevels-2*(l+add)+1,3,
	          point3(X(x1), Y(y1), X(x2), Y(y2), 
  		  X(pp[2]->x), Y(pp[2]->y)), FILLED);
	      }
  	      DrawLines(ptr,black,50+2*nlevels-2*(l+add), 2,
	        point2(X(x1), Y(y1), X(x2), Y(y2)));
            } else {
  	      DrawLines(ptr,colors[l+add],50+2*nlevels-2*(l+add), 2,
	        point2(X(x1), Y(y1), X(x2), Y(y2)));
	    }
 	  }
        }
      }
    }
  }
  if (type == XFIG) {
    fclose(file);
  } else {
#ifdef WITH_GD
    Tcl_Obj *tmp;
    gdSink sink;
    struct tmp context;
    context.bytes = 0;
    context.alloc = 10*w*h;
    context.content = malloc(context.alloc);    
    sink.context = (void*)&context;
    sink.sink = dump;
    gdImageGifToSink(img, &sink);
    
    tmp = Tcl_NewByteArrayObj(context.content, context.bytes);
    Tcl_SetObjResult(interp,tmp);
    free(context.content);
    gdImageDestroy(img);
#endif  
  }
DONE:
  free((char*)levels);
  free((char*)pp);
  free((char*)p);
  free((char*)t[0]);
  free((char*)t);
  free((char*)xpos);
  free((char*)ypos);
  free((char*)z[0]);
  free((char*)z);
  return TCL_OK;
}
