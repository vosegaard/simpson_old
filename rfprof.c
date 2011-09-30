/*
    Reads rf inhomogeneity files
    Copyright (C) 2007 Zdenek Tosner

    This file is part of the SIMPSON General NMR Simulation Package

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
    
    The rf inhomogeneity file must have the following data format:

       <N> <M>
       <scalefactor1 chan1> <scalefactor1 chan2> ... <weight1>
       <scalefactor2 chan1> <scalefactor2 chan2> ... <weight2>
       <scalefactor3 chan1> <scalefactor3 chan2> ... <weight3>
       ...
       <scalefactorN chan1> <scalefactorN chan2> ... <weightN>

    where <N> is the number of rf isochromats, <M> is the number of rf channels (is
    M is set to 1 the same profile is assumed on all channels), scalefactor is the 
    rf scale factor typically ranging from 0 to 2 (nominal field = 1). 
    Weight is relative occurance of the isochromat; sum of all weights doesn't need 
    to sum up to 1, program takes care of that...
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "matrix_new.h"
#include "defs.h"

/****
 * ZT: reads in rf profile from a file
 ****/
void read_rfproffile(const char* name,int nchan, mv_double** rfscalefactor,mv_double** rfweight,double* sumweight)
{
  FILE* fp;
  char fname[256];
  int N,i,j,M;
  int ver;
  mv_double *Xrfscalefactor,*Xrfweight, *Xdummatrix;
  double dum;
  
  ver= verbose & VERBOSE_RFPROF;
  strcpy(fname,name);
  if (!strcmp(fname,"")) {
    Xdummatrix=double_matrix_alloc(1,nchan);
    Xrfweight=double_vector_alloc(1);
    for (i=0;i<nchan;i++) {
      Xdummatrix->data[i]=1.0;
    }
    Xrfweight->data[0]=1.0;
    *sumweight=1.0;
    *rfscalefactor=Xdummatrix;
    *rfweight=Xrfweight;
    return;
  }
  
  
  *sumweight = 0.0;

#ifdef UNIX
  if (name[0] == '~') {
    char* p=getenv("HOME");
    if (p != NULL) {
      strcpy(fname,p);
      strcat(fname,&name[1]);
    }
  }
#endif
  if (!(fp=fopen(fname,"r"))) {
    strcat(fname,".rf");
    fp=fopen(fname,"r");
  }
  
  if (!fp) {
    fprintf(stderr,"error: unable to open file '%s'\n",fname);
    fprintf(stderr,"\n");
    exit(1);
  } 
  
  if (ver) printf("loading external rfprofile file '%s'\n",fname);

  if (fscanf(fp,"%d%d",&N,&M) != 2) {
    fprintf(stderr,"unable to read rfprofile from file '%s'\n",fname);
    exit(1);
  }
  Xrfscalefactor=double_matrix_alloc(N,M);
  Xrfweight=double_vector_alloc(N);
  Xdummatrix=double_matrix_alloc(N,nchan);

  for (i=1;i<=N;i++) {
    for (j=1;j<=M;j++) {
       if (fscanf(fp,"%lg",&dum) != 1) {
          fprintf(stderr,"error: unable to read line %d, column %d in file %s\n",i,j,fname);
       }
       double_put_elem(Xrfscalefactor,i,j,dum);
    }
    if (fscanf(fp,"%lg",&dum) != 1) {
          fprintf(stderr,"error: unable to read weight on line %d in file %s\n",i,fname);
    }
    if (dum == 0.0) {
         fprintf(stderr,"error: rf isochromat number %d in file '%s' has zero weight\n",i,name);
         exit(1);
    }
    double_put_elem(Xrfweight,i,1,dum);
    *sumweight += dum;
  }

  fclose(fp);

  for (i=1;i<=nchan;i++) {
    for (j=1;j<=N;j++) {
      if (i>M) {
	double_put_elem(Xdummatrix,j,i,1.0);
      } else {
	double_put_elem(Xdummatrix,j,i,double_get_elem(Xrfscalefactor,j,i));
      }
    }
  }
  *rfscalefactor=Xdummatrix;
  *rfweight=Xrfweight;
  double_matrix_free(Xrfscalefactor);
  
 if (ver) printf("Reading rfprofile done!'\n");
}

 
