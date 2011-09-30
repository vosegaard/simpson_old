/*
    Tcl initialization and main program
    Copyright (C) 1999 Mads Bak, 2000-2005 Thomas Vosegaard

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
    
    This is the main program that initates the Tcl commands declared
    in various source code files via the 'tclcmd' commands.
    Evaluates the 'main.tcl' Tcl code that loads other statically linked
    Tcl code and evaluates the input file. Next step in program
    flow is main.tcl.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tcl.h>
#include <unistd.h>
#include "tclutil.h"
#include "tclcmd.h"
#include "defs.h"
#include "config.h"
#include "cm_new.h"

#include "timing.h"
#ifdef TIMING
struct timezone timing_tz;
struct timeval  timing_tv1;
struct timeval  timing_tv2;
#endif
extern int AddContourCommand(Tcl_Interp* interp);
int verbose=0;
int various=0;

/* #ifdef LIBSIMPSON */
int Simpson_Init(Tcl_Interp *interp)
{
  DECLARE_TCL_COMMANDS;
  AddContourCommand(interp);
  return 1;
}
/* #else */

int mainargc=0;

/* 
  NOT_STATIC is without the TCL files complied into the binary 
  Usually not defined, but the code is leaved here to facilitate
  debugging of the tcl-code.
*/
#ifdef NOT_STATIC
char* findfile(char* file)
{
  static char path[256];
  char* libpath;

  libpath=getenv("SIMPSON_LIB");

  if (!libpath) {
    fprintf(stderr,"error: SIMPSON_LIB environment variable must be set\n");
    exit(1);    
  }
  strcpy(path,libpath);
  strcat(path,"/");
  strcat(path,file);
  return path;
}
#else /* NOT_STATIC */
extern char main_tcl[];
extern char init_tcl[];
#endif /* NOT_STATIC */

/* Required to link with f2c'ed fortran routines */
int MAIN__ () {return 0;}



int TclAppInit(Tcl_Interp* interp)
{
  char buf[256];

  if (Tcl_Init(interp) == TCL_ERROR) return TCL_ERROR;

  TclSetSignalHandler(interp,"signalhandler");

  DECLARE_TCL_COMMANDS;
  AddContourCommand(interp);

  if (NULL == Tcl_SetVar(interp,"simpson_version",VERSION,
      TCL_GLOBAL_ONLY | TCL_LEAVE_ERR_MSG)) {
    fprintf(stderr,"error: %s\n",interp->result);
    return TCL_ERROR;
  }
  sprintf(buf,"%d",mainargc);
  if (NULL == Tcl_SetVar(interp,"mainargc",buf,
      TCL_GLOBAL_ONLY | TCL_LEAVE_ERR_MSG)) {
    fprintf(stderr,"error: %s\n",interp->result);
    return TCL_ERROR;
  }
#ifdef NOT_STATIC
  if (Tcl_EvalFile(interp,findfile("main.tcl")) == TCL_ERROR) {
    fprintf(stderr,"error: %s\n",interp->result);
    return TCL_ERROR;
  }
#else
  if (Tcl_Eval(interp,main_tcl) != TCL_OK) {
    fprintf(stderr,"error: %s\n",interp->result);
    return TCL_ERROR;
  }
#endif

  /* ZT: memory clean-up */
  free_wsp(); free_mv_static();
  exit(0);
  return TCL_OK;
}        

FILE* _IO_stderr_;
FILE* _IO_stdout_;
FILE* _IO_stdin_;
#include <errno.h>
int main (int argc,char *argv[])
{
   _IO_stderr_=stderr;
   _IO_stdout_=stdout;
   _IO_stdin_=stdin;

#ifdef MPI
   /* RA: vars for MPI env. */
   //char hostname[80];
   int a;
   char **b;
   MPI_Init(&a, &b);
   errno = 0;
   //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   //MPI_Comm_size(MPI_COMM_WORLD, &np);
   //gethostname(hostname, 80);
   //printf("Welcome to simpson, I am process %d of %d. Im running on %s\n", rank, np, hostname);
#endif

   mainargc=argc;
   Tcl_Main(argc, argv, TclAppInit);

   return 0;
}

/* #endif */
