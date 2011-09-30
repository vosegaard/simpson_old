/*
    Alderman Tent averaging for static CSA calculations.
    Copyright (C) 1999 Jimmy T. Rasmussen, Mads Bak

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
    
*/
#include <tcl.h>
#include "iodata.h"


void csastat_triangle(FD* desc,double f1,double f2,double f3,double h)
{
   int i1,i2,i3,i;
   double d,f;
   double2* dat;

   dat=(double2*)desc->data;

   i1=FD_INDEX(desc,f1);
   i2=FD_INDEX(desc,f2);
   i3=FD_INDEX(desc,f3);

   d=h/(f2-f1);
   for (i=i1;i<i2;i++) {
      f=FD_FREQ(desc,i);
      if (f > f2) break;
      dat[i].re += (f-f1)*d;
   }
   dat[i].re += h;
   
   d=h/(f3-f2);
   for (i=i2+1;i<=i3;i++) {
      f=FD_FREQ(desc,i);
      if (f > f3) break;
      dat[i].re += (f3-f)*d;
   }
}

extern FD** fd;
extern int nfd;

int tclStent(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  FD* f;
  int fidN;
  double2* dat;

  if (argc != 2) {
    interp->result="usage: stent <desc>";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR) {
    interp->result="stent: argument 1 must be integer <descriptor>";
    return TCL_ERROR;
  }
  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL) {
    sprintf(interp->result,"stent: descriptor %d was not previously loaded\n",fidN);
    return TCL_ERROR;
  }
  f=fd[fidN];
  dat=(double2*)f->data;
  
  csastat_triangle(f,-200,-110,150,50);
  return TCL_OK;
}

void tclcmd_csastat(Tcl_Interp* interp)
{
  Tcl_CreateCommand(interp,"stent",tclStent,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
}
