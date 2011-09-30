/*
    Reads powder averaging/crystallite files
    Copyright (C) 1999 Mads Bak

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
    
    
    The crystallite file must have the following data format:

       <N>
       <alpha1> <beta1> <weight1>
       <alpha2> <beta2> <weight2>
       <alpha3> <beta3> <weight3>
       ...
       <alphaN> <betaN> <weightN>

    where <N> is the number of crystallites, alpha is the 
    alpha Euler angle ranging from 0 to 360 degrees, beta is
    the beta Euler angle ranging from 0 to 180 degrees.
    Weight is the solid angle covered by the crystallite
    and all weights must sum up to 1. Be avare that the number
    of significant digits must large, otherwise errors can be
    introduced.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "matrix_new.h"
#include "cryst.h"
#include "defs.h"
/* AB::beg */
/* Has to be included since the flag "use_3_angle_set" is used. */
#include "sim.h"
/* AB::end */     

/* Pass the simulation as parameter (instead of just the filename)   */
/* because we need to access the flag "s->use_3_angle_set".          */
/* Pass the struct omega containing the Euler angles instead of just */
/* alpha and beta.                                                   */

mv_double * read_crystfile(Sim* s)
{
  FILE* fp;
  char fname[256];
  CRYSTALLITE* c;
  int N,i,j;
  int n=0,ver;
/* Use struct for Euler angles */
  mv_double *Xomega;
  int file_status;
  double da, db, dg, dw;
  
  ver= verbose & VERBOSE_POWDER;
  strcpy(fname,s->crystfile);
  strcat(fname,"_cryst");
  while (strlen(cryst_names[n])) {
    if (!strcmp(cryst_names[n],fname)) {
       if (ver) {
         printf("found internal crystallite file '%s'\n",s->crystfile);
         printf("to overwrite with external file, specify './%s' instead\n",s->crystfile);
       }
       N=cryst_numbers[n];
       if (ver) printf("crystallites: %d\n",N);
       Xomega = double_matrix_alloc(N,4);
       c=cryst_pointers[n];
       for (i=0;i<N;i++) {
         j=i+1;
/* Set the powder angles, gamma is set to zero for the internal powder files */
	 double_put_elem(Xomega,j,1,c[i].alpha);
	 double_put_elem(Xomega,j,2,c[i].beta);
	 double_put_elem(Xomega,j,3,0.0);
	 double_put_elem(Xomega,j,4,c[i].weight);
         if (c[i].weight == 0.0) {
           fprintf(stderr,"error: crystallite number %d in file '%s' has zero weight\n",i+1,s->crystfile);
           exit(1);
         }
         if (ver) 
      	  printf("%5d %15g %15g %15g %15g\n",i+1, c[i].alpha,c[i].beta,0.0,c[i].weight);
       }
       return Xomega;
    }
    n++;
  }

  strcpy(fname,s->crystfile);

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
    strcat(fname,".cry");
    fp=fopen(fname,"r");
  }
  
  if (!fp) {
    int n=0;
    int l,tl=0;

    fprintf(stderr,"error: unable to open file '%s'\n",fname);
    fprintf(stderr,"internal crystallite files are: \n");
    while (strlen(cryst_names[n])) {
       char nam[256];
       strcpy(nam,cryst_names[n++]);
       l=strlen(nam)-6;
       nam[l]=0;
       tl += l;
       fprintf(stderr,"%s ",nam);
       if (tl > 32) {
         fprintf(stderr,"\n");
         tl=0;
       }
    }
    fprintf(stderr,"\n");
    exit(1);
  }
  
  if (ver) printf("loading external crystallite file '%s'\n",fname);

  if (fscanf(fp,"%d",&N) != 1) {
    fprintf(stderr,"unable to read crystallite from file '%s'\n",fname);
    exit(1);
  }
  Xomega = double_matrix_alloc(N,4);
  for (i=1;i<=N;i++) {
/* Read in 3 or 2 angle set depending on the flag */ 
    if( s->use_3_angle_set == 1 ) {
      file_status = (fscanf(fp,"%lg%lg%lg%lg",&da,&db,&dg,&dw) != 4);
    } else {
      dg = 0.0;
      file_status = (fscanf(fp,"%lg%lg%lg",&da,&db,&dw) != 3);
    }
    if( file_status ) {
       fprintf(stderr,"error: unable to read line %d in file %s\n",i,fname);
    }
    if (dw == 0.0) {
      fprintf(stderr,"error: crystallite number %d in file '%s' has zero weight\n",i,s->crystfile);
      exit(1);
    }
    double_put_elem(Xomega,i,1,da);
    double_put_elem(Xomega,i,2,db);
    double_put_elem(Xomega,i,3,dg);
    double_put_elem(Xomega,i,4,dw);

    if (ver) 
      printf("%5d %15g %15g %15g %15g\n",i, da,db,dg,dw);
  }  

  fclose(fp);
  return Xomega;
}
     

mv_double * read_crystfile_byname(char* crystname)
{
  FILE* fp;
  char fname[256];
  CRYSTALLITE* c;
  int N,i,j;
  int n=0,ver;
/* Use struct for Euler angles */
  mv_double *Xomega;
  int file_status;
  double da, db, dg, dw;
  
  ver= verbose & VERBOSE_POWDER;
  strcpy(fname,crystname);
  strcat(fname,"_cryst");
  while (strlen(cryst_names[n])) {
    if (!strcmp(cryst_names[n],fname)) {
       if (ver) {
         printf("found internal crystallite file '%s'\n",crystname);
         printf("to overwrite with external file, specify './%s' instead\n",crystname);
       }
       N=cryst_numbers[n];
       if (ver) printf("crystallites: %d\n",N);
       Xomega = double_matrix_alloc(N,4);
       c=cryst_pointers[n];
       for (i=0;i<N;i++) {
         j=i+1;
/* Set the powder angles, gamma is set to zero for the internal powder files */
	 double_put_elem(Xomega,j,1,c[i].alpha);
	 double_put_elem(Xomega,j,2,c[i].beta);
	 double_put_elem(Xomega,j,3,0.0);
	 double_put_elem(Xomega,j,4,c[i].weight);
         if (c[i].weight == 0.0) {
           fprintf(stderr,"error: crystallite number %d in file '%s' has zero weight\n",i+1,crystname);
           exit(1);
         }
         if (ver) 
      	  printf("%5d %15g %15g %15g %15g\n",i+1, c[i].alpha,c[i].beta,0.0,c[i].weight);
       }
       return Xomega;
    }
    n++;
  }

  strcpy(fname,crystname);

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
    strcat(fname,".cry");
    fp=fopen(fname,"r");
  }
  
  if (!fp) {
    int n=0;
    int l,tl=0;

    fprintf(stderr,"error: unable to open file '%s'\n",fname);
    fprintf(stderr,"internal crystallite files are: \n");
    while (strlen(cryst_names[n])) {
       char nam[256];
       strcpy(nam,cryst_names[n++]);
       l=strlen(nam)-6;
       nam[l]=0;
       tl += l;
       fprintf(stderr,"%s ",nam);
       if (tl > 32) {
         fprintf(stderr,"\n");
         tl=0;
       }
    }
    fprintf(stderr,"\n");
    exit(1);
  }
  
  if (ver) printf("loading external crystallite file '%s'\n",fname);

  if (fscanf(fp,"%d",&N) != 1) {
    fprintf(stderr,"unable to read crystallite from file '%s'\n",fname);
    exit(1);
  }
  Xomega = double_matrix_alloc(N,4);
  for (i=0;i<N;i++) {
    dg = 0.0;
    file_status = (fscanf(fp,"%lg%lg%lg",&da,&db,&dw) != 3);
    if( file_status ) {
       fprintf(stderr,"error: unable to read line %d in file %s\n",i,fname);
    }
    if (dw == 0.0) {
      fprintf(stderr,"error: crystallite number %d in file '%s' has zero weight\n",i,crystname);
      exit(1);
    }
    double_put_elem(Xomega,i,1,da);
    double_put_elem(Xomega,i,2,db);
    double_put_elem(Xomega,i,3,dg);
    double_put_elem(Xomega,i,4,dw);

    if (ver) 
      printf("%5d %15g %15g %15g %15g\n",i+1, da,db,dg,dw);
  }  

  fclose(fp);
  return Xomega;
}


