/*
    Tcl/C utillity routines
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
    
    Makes available som Tcl function to 
      - evaluate statically linkes Tcl code
      - get parameters from the 'par' array.
      - signal handling so Ctrl-C is handled properly
*/


#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <signal.h>
#include <stdarg.h>
#include "matrix_new.h"
#include "tclutil.h"
#include "defs.h"


int tclEvalInternalCode(ClientData data,Tcl_Interp* interp,
      int argc, char *argv[])
{
  int found=0;
  TCLCODE *p;
  
  if (argc != 2)
    return TclError(interp,"Usage: evalinternalcode <variable name>");
  p=tclcode_pointers;
  while (*(p->name)) {
     if (*(p->name) == 0) {
       fprintf(stderr,"error: evalinternalcode: linker error\n");
       exit(-1);
     }
     if (!strcmp(argv[1],p->name)) {
       if (Tcl_Eval(interp,p->code) != TCL_OK) {
         fprintf(stderr,"error: %s\n",interp->result);
         return TCL_ERROR;
       }     
       found=1;
       break;
     }
     p++;
  }
  if (!found) {
    fprintf(stderr,"evalinternalcode: error, cannot find data for file '%s.tcl'\n",argv[1]);
    exit(-1);
  }
  return TCL_OK;
}


int TclAppendResult(Tcl_Interp* interp,const char* format, ...)
{
   char buffer [512];
   va_list argptr;
   va_start(argptr, format);
   vsprintf(buffer, format, argptr);
   va_end(argptr);
   Tcl_AppendElement(interp,buffer);
   return TCL_OK;
}

int TclSetResult(Tcl_Interp* interp,const char* format, ...)
{
   char buffer [512];
   va_list argptr;
   va_start(argptr, format);
   vsprintf(buffer, format, argptr);
   va_end(argptr);
   Tcl_SetResult(interp,buffer,TCL_VOLATILE);
   return TCL_OK;
}

int TclError(Tcl_Interp* interp,const char* format, ...)
{
   char buffer [512];
   va_list argptr;
   va_start(argptr, format);
   vsprintf(buffer, format, argptr);
   va_end(argptr);
   Tcl_SetResult(interp,buffer,TCL_VOLATILE);   
   return TCL_ERROR;
}


int TclGetInt(Tcl_Interp* interp,char *aryname,char* varname,
                     int mustexist,int defval)
{
  int  val;
  char* src;

  if ((src=Tcl_GetVar2(interp,aryname,varname,0)) == NULL) {
     if (mustexist) {
       fprintf(stderr,"error: could not read integer variable %s(%s)\n",aryname,varname);
       exit(1);
     }
     if (verbose & VERBOSE_PAR)
       printf("integer variable %s in array %s is set to default value %d\n",varname,aryname,defval);

     return defval;
  }

  if (Tcl_GetInt(interp,src,&val) != TCL_OK) TclError(interp,"GetInt(2)");
 
  if (verbose & VERBOSE_PAR)
    printf("integer variable %s in array %s is set to %d\n",varname,aryname,val);
  return val;
}

double TclGetDouble(Tcl_Interp* interp,char *aryname,char* varname,
                     int mustexist,double defval)
{
  double val;
  char* src;
  
  if ((src=Tcl_GetVar2(interp,aryname,varname,0)) == NULL) {
     if (mustexist) {
       fprintf(stderr,"error: could not read double variable %s(%s)\n",aryname,varname);
       exit(-1);
     }
     if (verbose & VERBOSE_PAR)
       printf("double variable %s in array %s is set to default value %g\n",varname,aryname,defval);
     return defval;
  }
  if (Tcl_GetDouble(interp,src,&val) != TCL_OK) TclError(interp,"GetInt(2)");
 
  if (verbose & VERBOSE_PAR)
    printf("double variable %s in array %s is set to %g\n",varname,aryname,val);
  return val;
}

char* TclGetString(Tcl_Interp* interp,char *dst,char* aryname,char* varname,
                     int mustexist,char* defval)
{
  char* src;

  if ((src=Tcl_GetVar2(interp,aryname,varname,0)) == NULL) {
     if (mustexist) {
       fprintf(stderr,"error: could not read string variable %s(%s)\n",aryname,varname);
       exit(-1);
     }
     if (verbose & VERBOSE_PAR)
       printf("string variable %s in array %s is set to default value %s\n",varname,aryname,defval);
     strcpy(dst,defval);
     return dst;
  }
  strcpy(dst,src);

  if (verbose & VERBOSE_PAR) {
    printf("string variable %s in array %s is set to %s\n",varname,aryname,dst);
  }
  return dst;
}

double* TclGetVector(Tcl_Interp* interp,char* aryname,char* varname,
                     int mustexist,double* defval) 
{
  int i,argc;
  char** argv;
  char* list;
  double* v;
  
  list=Tcl_GetVar2(interp,aryname,varname,0);
  if (!list) {
     if (mustexist) {
       fprintf(stderr,"error: could not read vector variable %s(%s)\n",aryname,varname);
       exit(-1);
     }
     if (verbose & VERBOSE_PAR) {
       printf("vector variable %s in array %s is set to default value ",varname,aryname);
       if (defval != NULL) {
         for (i=1;i<=LEN(defval);i++)
           printf("%f ",defval[i]);
       } else {
           printf("<null>\n");
       }
     }
     return defval;
  }
  if (Tcl_SplitList(interp,list,&argc,&argv) != TCL_OK) TclError(interp,"GetVector(2)");
  if (!argc) return NULL;
  v = double_vector(argc);
  for (i=0;i<argc;i++) {
    if (Tcl_GetDouble(interp,argv[i],&v[i+1]) != TCL_OK) TclError(interp,"GetVector(3)");
  }
  free(argv);
  if (verbose & VERBOSE_PAR) {
    printf("vector variable %s in array %s is set to value ",varname,aryname);
    for (i=1;i<=LEN(v);i++)
       printf("%f ",v[i]);
  }
  return v;
}


char* sigtxt[] = {
"-", /* 0 */
"-", /* 1 */
"Interrupt from keyboard", /* 2 */
"Quit from keyboard", /* 3 */
"-", /* 4 */
"-", /* 5 */
"Abort", /* 6 */
"-", /* 7 */
"-", /* 8 */
"Kill signal", /* 9 */
"-", /* 10 */
"-", /* 11 */
"-", /* 12 */
"-", /* 13 */
"-", /* 14 */
"Termination signal", /* 15 */
};

int lastsig;

int TclSignalHandler(ClientData clientData,Tcl_Interp *interp,int code)
{
   char buf[256];
   if (interp == NULL) {
     fprintf(stderr,"signalhandler was called with code %d "
                    "when no Tcl interpreter was available\n",code);
     return TCL_OK;
   }
   sprintf(buf,"signalhandler %d \"%s\"",lastsig,sigtxt[lastsig]);
   return Tcl_Eval(interp,buf);
}           

Tcl_AsyncHandler asynchandler;

void signal_handler(int sig)
{
  lastsig=sig;
  Tcl_AsyncMark(asynchandler);
}

void TclSetSignalHandler(Tcl_Interp* interp,char* function)
{
  asynchandler=Tcl_AsyncCreate(TclSignalHandler,(ClientData)NULL);

#ifdef SIGINT
  signal(SIGINT,signal_handler);
#else
  signal(SIGBREAK,signal_handler);
#endif
}
#define BUFLEN 1024

int TclAppendMatrix(Tcl_Interp* interp,mv_complx * m)
{
   int i,j,r,c;
   char num[64];
   char buf[BUFLEN];
   int nbuf;
   complx z;
   
   r = m->row;   
   c = m->col;
   for (i=0;i<r;i++) {
     nbuf=0;
     buf[0]=0;
     for (j=0;j<c;j++) {
        z = m->data[i+j*r];
        sprintf(num,"{%g %g}",z.re,z.im);
        nbuf += strlen(num);
        if (nbuf >= BUFLEN) {
          Tcl_SetResult(interp,"getmatrix: internal buffer overflow\n",TCL_STATIC);
          return TCL_ERROR;          
        }
        if (j != 0) strcat(buf," ");
        strcat(buf,num);
     } 
     Tcl_AppendElement(interp,buf);
  }
  return TCL_OK;
}

int TclAppendRealMatrix(Tcl_Interp* interp,mv_double * m)
{
   int i,j,r,c;
   char num[64];
   char buf[BUFLEN];
   int nbuf;
   
   r = m->row;   
   c = m->col;
   for (i=1;i<=r;i++) {
     nbuf=0;
     buf[0]=0;
     for (j=1;j<=c;j++) {
        sprintf(num,"%g",m->data[i+j*r]);
        nbuf += strlen(num);
        if (nbuf >= BUFLEN) {
          Tcl_SetResult(interp,"getmatrix: internal buffer overflow\n",TCL_STATIC);
          return TCL_ERROR;          
        }
        if (j != 1) strcat(buf," ");
        strcat(buf,num);
     } 
     Tcl_AppendElement(interp,buf);
  }
  return TCL_OK;
}



void tclcmd_tclutil(Tcl_Interp* interp)
{
  Tcl_CreateCommand(interp,"evalinternalcode",tclEvalInternalCode,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
}
