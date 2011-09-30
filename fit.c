/*
    Minimization routines, wrappers and drivers
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
  

    Makes available the 'fit' Tcl command that can be used
    for generic minimization. Available methods are simplex
    and hookejeeves algorithms.      
*/

#include <stdio.h>
#include <memory.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <tcl.h>
#include "tclutil.h"

#define PARMAX 256

int ndim; /* number used parameters */
int nb[PARMAX]; /* translation between all, and used vectors */
double (*origfunc)(double[]); /* the function to be minimized */
double xvalues[PARMAX];

double mdriver(double x[])  /* wrapper function between used and all*/
{ 
  int i;
  double r;
  for (i=1;i<=ndim;i++) {
    xvalues[nb[i]]=x[i]; /* translate selected parameters to all parameters */
  }
  r=(*origfunc)(xvalues);
  return r;
}

double simplex(double (*func)(double[]),int npar,int *usewho,
               double *start,double *scale)
{
  int simplx(double (*func)(double* v),int n,double *x,
            double *scale, double *fx);
  int i;
  double xi[PARMAX+1],y[PARMAX+1];
  double res;
  
  origfunc=func;  
  ndim=0;

  for (i=1;i<=npar;i++) {
    if (usewho[i]){
      ndim++;    /* count number of parameters to minimize on */
      nb[ndim]=i; /* this parameter will be next real parameter (at place i)*/
    }
    xvalues[i]=start[i]; /* save all start parameters */
  }
  for (i=1;i<=ndim;i++) {
    xi[i]=scale[nb[i]];
    y[i]=start[nb[i]];
  }
  simplx(mdriver,ndim,y,xi,&res);
  for (i=1;i<=ndim;i++) {
    start[nb[i]]=y[i];
  }
  return 0;
}

int iter;
int npar;
char function[256];
Tcl_Interp* intrp;
char name[512][64];
char array[256];

double func (double x[])
{
  int i;
  double value;
  char buf2[256];
  char buf[2048];

  sprintf(buf,"%d",++iter);
  if (NULL == Tcl_SetVar2(intrp,array,"iter",buf,
       TCL_GLOBAL_ONLY | TCL_LEAVE_ERR_MSG)) {
    fprintf(stderr,"Error: '%s'\n",intrp->result);
    exit(1);
  }
  
  strcpy(buf,function);
  strcat(buf," {");
  for (i=1;i<=npar;i++) {      
    sprintf(buf2," { %s %g }",name[i],x[i]);
    strcat(buf,buf2);
  }
  strcat(buf," }");

  if (Tcl_Eval(intrp,buf) != TCL_OK) {
    fprintf(stderr,"Error: '%s'\n",intrp->result);
    exit(1);
  }
  
  if (Tcl_GetDouble(intrp,intrp->result,&value) != TCL_OK) {
    fprintf(stderr,"Error: '%s'\n",intrp->result);
    exit(1);
  }
  return value;
}


int tclFit(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int i;
  int npar2;
  double val[512];
  double scal[512];
  int   used[512];
  char method[256];
  char *p,**pp,**pp2;

  if (argc != 2) {
    fprintf(stderr,"%s: specify array with parameters\n",argv[0]);
    return 1;
  }
  iter=0;
  strcpy(array,argv[1]);
  intrp = interp;

  TclGetString(interp,function,array,"function",1,"");

  TclGetString(interp,method,array,"fitmethod",0,"simplex");  
  
  p=Tcl_GetVar2(interp,array,"values",0);
  if (p == NULL)
    TclError(interp,"array must contain a variable 'values'");

  if (Tcl_SplitList(interp,p,&npar,&pp) != TCL_OK) {
    fprintf(stderr,"%s\n",interp->result);
    exit(1);
  }
  for (i=0;i<npar;i++) {
    if (Tcl_SplitList(interp,pp[i],&npar2,&pp2) != TCL_OK) {
      fprintf(stderr,"%s\n",interp->result);
      exit(1);
    }
    if (npar2 != 4)
       TclError(interp,"invalid number of parameters in 'values'"); 
    strcpy(name[i+1],pp2[0]);
    if (Tcl_GetDouble(interp,pp2[1],&val[i+1]) != TCL_OK)
       TclError(interp,"getdouble(1)");
    if (Tcl_GetDouble(interp,pp2[2],&scal[i+1]) != TCL_OK)
       TclError(interp,"getdouble(2)");    
    if (Tcl_GetInt(interp,pp2[3],&used[i+1]) != TCL_OK)
       TclError(interp,"getint(3)");
    free(pp2);
  }
  free(pp);
  if (!strcmp(method,"simplex")) {
    simplex(func,npar,used,val,scal);
  } else {
    return TclError(interp,"fit: unknown method '%s' (known method: simplex)",method);
  }
  return TCL_OK;
}


void tclcmd_fit(Tcl_Interp* interp)
{
  Tcl_CreateCommand(interp,"fit",tclFit,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
}

