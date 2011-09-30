/*
    Flatcoil restrictions: simulation and plotting
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
#include "complx.h"
#include "matrix_new.h"
#include "cm_new.h"
#include "tclutil.h"
#include "ham.h"

#define STEP 2
#define BDIM (180*STEP)
#define ADIM (360*STEP)
/* postscript linewidth should be 1.0/STEP */
#define LINEWIDTH "0.5"

#define D2R (M_PI/180.0)

#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))


void Rtensor(complx* R,double waniso,double eta)
{
  R[0]=R[4]= Complx( -waniso*eta/sqrt(6.0) ,0.0);
  R[1]=R[3]= Complx(0.0,0.0);
  R[2]=      Complx(  waniso,0.0);
}

void getdouble2(Tcl_Interp* interp, char* ary, char* str,
                double* var1,double* var2)
{
  double* vec;
  
  vec=TclGetVector(interp,ary,str,1,0);
  if (LEN(vec) != 2) {
    fprintf(stderr,"error: variable name '%s' must have two arguments\n",str);
    exit(1);
  }
  *var1= vec[1];
  *var2= vec[2];
  free_double_vector(vec);
}

double fread_val(FILE* fp)
{
   float f;   
   if (fread(&f,sizeof(float),1,fp) != 1) {
     fprintf(stderr,"error: unable to read float from file\n");
     exit(-1);
   }
   return (double)f;
}

void fwrite_val(FILE* fp,double d)
{
   float f=(float)d;   
   if (fwrite(&f,sizeof(float),1,fp) != 1) {
     fprintf(stderr,"error: unable to write float to file\n");
     exit(-1);
   }
}

int tclResShift(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  char outname[256],ary[128];
  int A,B,i,s_step,alpha_ML_step,beta_ML_step,i_xx,i_yy,i_zz,ia;
  double s_xx,s_yy,s_zz,s_obs,d_xx,d_yy,d_zz,d_obs,astep,da;
  double tilt_L,costilt_L,d00_L,*omega_PM;
  double scal_max[ADIM+1][BDIM+1],scal_min[ADIM+1][BDIM+1];
  complx *d20_ML[ADIM+1][BDIM+1],*d2_PM,*d2_a,*R_P,*R_M,*R_M1;
  const int N=5;
  
  FILE *of;

  if (argc != 3) {
    interp->result = "Usage: resshift <array> <output file>";
    return TCL_ERROR;
  }
  strcpy(ary,argv[1]);
  strcpy(outname,argv[2]);

  getdouble2(interp,ary,"shift_xx",&s_xx,&d_xx);
  getdouble2(interp,ary,"shift_yy",&s_yy,&d_yy);
  getdouble2(interp,ary,"shift_zz",&s_zz,&d_zz);
  getdouble2(interp,ary,"shift_obs",&s_obs,&d_obs);
  s_step =TclGetInt(interp,ary,"shift_step",0,2);
  omega_PM=TclGetVector(interp,ary,"omega_PM",1,0);

  astep =TclGetInt(interp,ary,"gamma_PM_step",0,1);
  da =TclGetDouble(interp,ary,"gamma_PM_delta",0,0);
  if (da > 0.0 && astep == 1) astep=2;
  
  tilt_L =TclGetDouble(interp,ary,"tilt_L",0,0);
  alpha_ML_step=TclGetInt(interp,ary,"alpha_ML_step",1,0);
  beta_ML_step =TclGetInt(interp,ary,"beta_ML_step",1,0);

  costilt_L=cos(tilt_L*D2R);
  d00_L=0.5*(3.0*costilt_L*costilt_L-1.0);
  d2_PM = (complx*)malloc(N*N*sizeof(complx));  
  wigner2(d2_PM,omega_PM[1],omega_PM[2],omega_PM[3]);

  d_xx *= 2.0;
  d_yy *= 2.0;
  d_zz *= 2.0;

  for (A=0; A<ADIM; A += alpha_ML_step) {
    for (B=0; B<BDIM; B += beta_ML_step) {
       d20_ML[A][B] = (complx*)malloc(N*sizeof(complx));       
       wigner20(d20_ML[A][B],(double)A/(double)STEP,(double)B/(double)STEP);
       scal_max[A][B]= -1.0e10;
       scal_min[A][B]=  1.0e10;
    }
  }

  d2_a = (complx*)malloc(N*N*sizeof(complx));
  R_P = (complx*)malloc(N*sizeof(complx));
  R_M = (complx*)malloc(N*sizeof(complx));
  R_M1 = (complx*)malloc(N*sizeof(complx));
  
  
  for (i_xx=0;i_xx<s_step;i_xx++) { 
    double S_xx= s_xx+((double)i_xx/(double)(s_step-1)-0.5)*d_xx;
    for (i_yy=0;i_yy<s_step;i_yy++) {
      double S_yy= s_yy+((double)i_yy/(double)(s_step-1)-0.5)*d_yy;
      for (i_zz=0;i_zz<s_step;i_zz++) {
        double S_zz= s_zz+((double)i_zz/(double)(s_step-1)-0.5)*d_zz;
        double wiso  = (S_xx+S_yy+S_zz)/3.0;
        double waniso= S_zz-wiso;
        double eta   = (S_yy-S_xx)/waniso;

        Rtensor(R_P,waniso,eta);
        /* m_mulvm(R_M,R_P,d2_PM); */
	zgemv_("T",&N,&N,&CPLX1,d2_PM,&N,R_P,&INTONE,&CPLX0,R_M,&INTONE);

        for (ia=0;ia<astep;ia++) {
          if (astep == 1)
            /* m_copyv(R_M1,R_M); */
	    memcpy(R_M1,R_M,N*sizeof(complx));
          else {
            double a=((double)ia/(double)(astep-1)-0.5)*da;
            wigner2(d2_a,a,0.0,0.0);          
            /* m_mulvm(R_M1,R_M,d2_a); */
  	    zgemv_("T",&N,&N,&CPLX1,d2_a,&N,R_M,&INTONE,&CPLX0,R_M1,&INTONE);
          }

          for (A=0; A<ADIM; A += alpha_ML_step) {
            for (B=0; B<BDIM; B += beta_ML_step) {
              complx *a=R_M1;
              complx *b=d20_ML[A][B];
              double w_l,scali,w_o = 0;
              /* Real part of the scalar product */
              for (i=0;i<N;i++) {
                 w_o += a[i].re * b[i].re - a[i].im * b[i].im;
              }
              w_l = w_o * d00_L;
              scali= (w_l+wiso)-s_obs;
              if (scali > scal_max[A][B]) scal_max[A][B]=scali;
              if (scali < scal_min[A][B]) scal_min[A][B]=scali;
            }
          }
        }
      }
    }
  }

  for (A=0; A<=ADIM; A += alpha_ML_step) {
    scal_max[A][BDIM]=scal_max[A][0];
    scal_min[A][BDIM]=scal_min[A][0];
  }
  for (B=0; B<=BDIM; B += beta_ML_step){
    scal_max[ADIM][B]=scal_max[0][B];
    scal_min[ADIM][B]=scal_min[0][B];
  }

  of=fopen(outname,"w");
  if (!of)
    return TclError(interp,"error: unable to create file '%s'",outname);

  fwrite_val(of,alpha_ML_step);
  fwrite_val(of,beta_ML_step);

  for (A=0; A<=ADIM; A += alpha_ML_step) {        
    for (B=0; B<=BDIM; B += beta_ML_step) {
       double mx=scal_max[A][B] + d_obs;
       double mn=scal_min[A][B] - d_obs;
       fwrite_val(of,min(mx,-mn));
    }
  }
  fclose(of);
  for (A=0; A<ADIM; A += alpha_ML_step) {
    for (B=0; B<BDIM; B += beta_ML_step) {
       free((char*)d20_ML[A][B]);
    }
  }
  free((char*)R_P);
  free((char*)R_M);
  free((char*)R_M1);
  free((char*)d2_PM);
  free((char*)d2_a);
  free_double_vector(omega_PM);
  return TCL_OK;
}

int tclResDipole(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  char outname[256],ary[128];
  int A,B,i,alpha_ML_step,beta_ML_step;
  double split_obs,delta_obs,split_max,delta_max;
  double tilt_L,costilt_L,d00_L,*omega_PM;
  double scal_max[ADIM+1][BDIM+1],scal_min[ADIM+1][BDIM+1];
  complx *d20_ML,*d2_PM;
  FILE *of;
  const int N=5;

  if (argc != 3) {
    interp->result = "Usage: resdipole <array> <output file>";
    return TCL_ERROR;
  }
  strcpy(ary,argv[1]);
  strcpy(outname,argv[2]);

  getdouble2(interp,ary,"dipole_obs",&split_obs,&delta_obs);
  getdouble2(interp,ary,"dipole_max",&split_max,&delta_max);
  omega_PM=TclGetVector(interp,ary,"omega_PM",1,0);

  tilt_L =TclGetDouble(interp,ary,"tilt_L",0,0);
  alpha_ML_step=TclGetInt(interp,ary,"alpha_ML_step",1,0);
  beta_ML_step =TclGetInt(interp,ary,"beta_ML_step",1,0);

  costilt_L=cos(tilt_L*D2R);
  d00_L=0.5*(3.0*costilt_L*costilt_L-1.0);

  d2_PM = (complx*)malloc(N*N*sizeof(complx));  
  wigner2(d2_PM,omega_PM[1],omega_PM[2],omega_PM[3]);

  d20_ML = (complx*)malloc(N*sizeof(complx));       
  for (A=0; A<ADIM; A += alpha_ML_step) {
    for (B=0; B<BDIM; B += beta_ML_step) {
      double w;
      complx c;
      wigner20(d20_ML,(double)A/(double)STEP,(double)B/(double)STEP);

      w=0.0;
      for (i=0;i<5;i++) {         
         c = d2_PM[2+i*5];
         w += c.re * d20_ML[i].re - c.im * d20_ML[i].im;
      }
      scal_max[A][B] = fabs(w * d00_L * (split_max+delta_max)) - split_obs + delta_obs;      
      scal_min[A][B] = fabs(w * d00_L * (split_max-delta_max)) - split_obs - delta_obs;
    }
  }

  for (A=0; A<=ADIM; A += alpha_ML_step) {
    scal_max[A][BDIM]=scal_max[A][0];
    scal_min[A][BDIM]=scal_min[A][0];
  }
  for (B=0; B<=BDIM; B += beta_ML_step){
    scal_max[ADIM][B]=scal_max[0][B];
    scal_min[ADIM][B]=scal_min[0][B];
  }

  of=fopen(outname,"w");
  if (!of)
    return TclError(interp,"error: unable to create file '%s'",outname);

  fwrite_val(of,alpha_ML_step);
  fwrite_val(of,beta_ML_step);
  for (A=0; A<=ADIM; A += alpha_ML_step) {        
    for (B=0; B<=BDIM; B += beta_ML_step) {
       double mx=scal_max[A][B];
       double mn=scal_min[A][B];
       fwrite_val(of,min(mx,-mn));
    }
  }
  fclose(of);

  free((char*)d2_PM);
  free((char*)d20_ML);
  return TCL_OK;
}


const char* header[] = {
"%!PS-Adobe-2.0",
"%%BoundingBox: 187 _yy5 393 591",
"147 199 translate 0.5 dup scale",
"/D { exch 360 exch sub exch  moveto 0 1 rlineto stroke } bind def",
"/L { /y2 exch def /y exch def /x exch 360 exch sub def",
"     x y moveto x y2 1 add lineto stroke } bind def",
"/Rshow {dup stringwidth pop neg 0 rmoveto show} def",
"/FT0 { /Times-Roman findfont 20 scalefont setfont} def",
"/FT1 { /Times-Roman findfont 25 scalefont setfont} def",
"/FS {/Symbol findfont 25 scalefont setfont} def",
"/FT {/Times-Roman findfont 18 scalefont setfont} def",
"/Sc { dup stringwidth pop -2 div 0 rmoveto show } def",
"150 750 translate -90 rotate 2 setlinewidth 1.7 dup scale",
"gsave 366 -10 translate FT 90 rotate",
"0 360 moveto (360) Rshow 0 0 moveto (0) Rshow FS",
"-10 177 moveto (a) Rshow  grestore",
"/axis { gsave translate",
"newpath 0 0 moveto 360 0 rlineto 0 180 rlineto",
"-360 0 rlineto closepath stroke /d 5 def",
"0 0 moveto 0 d neg rlineto stroke",
"90 0 moveto 0 d neg rlineto stroke",
"180 0 moveto 0 d neg rlineto stroke",
"270 0 moveto 0 d neg rlineto stroke",
"360 0 moveto 0 d neg rlineto stroke",
"360 0 moveto d 0 rlineto stroke",
"360 90 moveto d 0 rlineto stroke",
"360 180 moveto d 0 rlineto stroke",
"383 -4 translate 90 rotate FT",
"167 0 moveto (180) show 0 0 moveto (0) show FS",
"88 -15 moveto (b) show grestore } def",
"/SUP { 0 5 rmoveto show 0 -5 rmoveto } def",
"/f -15 def /g  40 def FT1",
"0 0 axis",
"/cliprect { newpath 0 0 moveto 360 0 rlineto 0 180 rlineto -360 0 rlineto",
"closepath clip newpath } def",
"0 setlinecap " LINEWIDTH " setlinewidth",
""};

const char* header4[] = {
"%!PS-Adobe-2.0",
"%%BoundingBox: 0 0 596 842",
"%%DocumentPaperSizes: A4",
"/D { exch 360 exch sub exch  moveto 0 1 rlineto stroke } bind def",
"/L { /y2 exch def /y exch def /x exch 360 exch sub def",
"     x y moveto x y2 1 add lineto stroke } bind def",
"/Rshow {dup stringwidth pop neg 0 rmoveto show} def",
"/FT0 { /Times-Roman findfont 20 scalefont setfont} def",
"/FT1 { /Times-Roman findfont 25 scalefont setfont} def",
"/FS {/Symbol findfont 25 scalefont setfont} def",
"/FT {/Times-Roman findfont 18 scalefont setfont} def",
"/Sc { dup stringwidth pop -2 div 0 rmoveto show } def",
"140 65 translate 2 setlinewidth 0.9 dup scale",
"gsave 366 -10 translate FT 90 rotate",
"0 360 moveto (360) Rshow 0 0 moveto (0) Rshow FS",
"-10 177 moveto (a) Rshow  grestore",
"/axis { gsave translate",
"newpath 0 0 moveto 360 0 rlineto 0 180 rlineto",
"-360 0 rlineto closepath stroke /d 5 def",
"0 0 moveto 0 d neg rlineto stroke",
"90 0 moveto 0 d neg rlineto stroke",
"180 0 moveto 0 d neg rlineto stroke",
"270 0 moveto 0 d neg rlineto stroke",
"360 0 moveto 0 d neg rlineto stroke",
"360 0 moveto d 0 rlineto stroke",
"360 90 moveto d 0 rlineto stroke",
"360 180 moveto d 0 rlineto stroke",
"383 -4 translate 90 rotate FT",
"167 0 moveto (180) show 0 0 moveto (0) show FS",
"88 -15 moveto (b) show grestore } def",
"/SUP { 0 5 rmoveto show 0 -5 rmoveto } def",
"/f -15 def /g  40 def FT1",
"gsave f g translate 0 0 moveto",
"90 rotate FT0 () SUP FT1 () show grestore",
"gsave f 210 g add translate 0 0 moveto 90 rotate",
"FT0 () SUP FT1 () show grestore",
"gsave f 420 g add 35 sub translate",
"0 0 moveto 90 rotate  FT0 () SUP FT1 () show FT0 () SUP",
"FT1 () show grestore",
"gsave f 630 g add translate 0 0 moveto 90 rotate",
"FT1 () show grestore",
"0 0   axis 0 210 axis 0 420 axis 0 630 axis",
"/cliprect { newpath 1 1 moveto 358 0 rlineto 0 178 rlineto -358 0 rlineto",
"closepath clip newpath } def",
"0 setlinecap " LINEWIDTH " setlinewidth",
""};

#define DATASIZE ((ADIM+2)*(BDIM+2))

double acc_i,acc_j,acc_j2;

void accum_beg(FILE *fo)
{
  acc_i=-1;
  acc_j=-1;
  acc_j2=-1;
}

void accum_end(FILE *fo)
{
  if (acc_i != -1) {  /* there was'nt a previous point */
    if (acc_j == acc_j2) { /* a point but not a line */
      fprintf(fo,"%g %g D\n",acc_i,acc_j);
    } else {
      fprintf(fo,"%g %g %g L\n",acc_i,acc_j,acc_j2);
    }
  }
  acc_i=-1;
  acc_j=-1;
  acc_j2=-1;  
}

void accum(FILE *fo,double i,double j)
{

   if (acc_i == -1) { /* this is first time */
     acc_i=i;
     acc_j=j;
     acc_j2=j;
   } else if (acc_i == i) { /* next step on a line */
     acc_j2=j;
   } else {  /* acc_i != i:   new i == new line */
     accum_end(fo);
     acc_i=i;
     acc_j=j;
     acc_j2=j;
   }
}

int tclResMakePs4(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int i,f,A,B;
  double d,alpha_ML_step,beta_ML_step;
  FILE *fo,*fi;

  if (argc != 6) {
    interp->result = "Usage: resmakeps4 <.out 1H shift> <.out 15N shift> <.out 1H-15N coup> <.out combined> <new .ps file>";
    return TCL_ERROR;
  }

  fo=fopen(argv[5],"w");
  if (!fo) {
    sprintf(interp->result,"error: unable to create file '%s'",argv[5]);
    return TCL_ERROR;
  }
  i=0;
  while (*header4[i] != 0) {
    fprintf(fo,"%s\n",header4[i++]);
  }
  for (f=1;f<=4;f++) {
    if (!strcmp(argv[f],"NONE")) continue;
    fi=fopen(argv[f],"r");
    if (!fi)
      return TclError(interp,"error: unable to open file '%s'",argv[f]);

    alpha_ML_step=fread_val(fi);
    beta_ML_step=fread_val(fi);

    fprintf(fo,"gsave 0 %d translate cliprect\n",(f-1)*210);
    accum_beg(fo);

    for (A=0; A<=ADIM; A += alpha_ML_step) {        
      for (B=0; B<=BDIM; B += beta_ML_step) {
         d=fread_val(fi);
         if (d >= 0.0)
           accum(fo,A/(double)STEP,B/(double)STEP);
         else
           accum_end(fo);
      }
    }
    accum_end(fo);
    fprintf(fo,"grestore\n");
    fclose(fi);
  }
  fprintf(fo,"showpage\n");
  fclose(fo);
  return TCL_OK;
}


int tclResMakePs(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int i,A,B;
  double d,alpha_ML_step,beta_ML_step;
  FILE *fo,*fi;

  if (argc != 3) {
    interp->result = "Usage: resmakeps <inputfile> <postscript-output>";
    return TCL_ERROR;
  }

  fi=fopen(argv[1],"r");
  if (!fi) {
    sprintf(interp->result,"error: unable to open file '%s'",argv[1]);
    return TCL_ERROR;
  }
  alpha_ML_step=fread_val(fi);
  beta_ML_step=fread_val(fi);

  fo=fopen(argv[2],"w");
  if (!fo) {
    sprintf(interp->result,"error: unable to create file '%s'",argv[2]);
    return TCL_ERROR;
  }

  i=0;
  while (*header[i] != 0) {
    fprintf(fo,"%s\n",header[i++]);
  }
  fprintf(fo,"gsave cliprect\n");
  accum_beg(fo);

  for (A=0; A<=ADIM; A += alpha_ML_step) {        
    for (B=0; B<=BDIM; B += beta_ML_step) {
       d=fread_val(fi);
       if (d >= 0.0)
         accum(fo,A/(double)STEP,B/(double)STEP);
       else
         accum_end(fo);
    }
  }
  accum_end(fo);
  fclose(fi);

  fprintf(fo,"grestore\nshowpage\n");
  fclose(fo);
  return TCL_OK;
}


#define OPER_AND 1
#define OPER_XOR 2
#define OPER_OR 3

void oper(char* fout,char* fin1,char* fin2,int oper)
{
  FILE *fi1,*fi2,*fo;
  int A,B,di1,dj1,di2,dj2;
  double d1,d2;
 
  fi1=fopen(fin1,"r");
  if (!fi1) {
    fprintf(stderr,"error: unable to open file '%s'\n",fin1);
    exit(1);
  }

  fi2=fopen(fin2,"r");
  if (!fi2) {
    fprintf(stderr,"error: unable to open file '%s'\n",fin2);
    exit(1);
  }

  fo=fopen(fout,"w");
  if (!fout) {
    fprintf(stderr,"error: unable to open file '%s'\n",fin2);
    exit(1);
  }
  
  di1=(int)fread_val(fi1);
  dj1=(int)fread_val(fi1);
  di2=(int)fread_val(fi2);
  dj2=(int)fread_val(fi2);


  if (di1 != di2 || dj1 != dj2) {
    fprintf(stderr,"error: files '%s' and '%s' are of different size\n",fin1,fin2);
    exit(-1);
  }

  fwrite_val(fo,(double)di2);
  fwrite_val(fo,(double)di1);

  if (oper == OPER_AND) {
    for (A=0; A<=ADIM; A += di1) {        
      for (B=0; B<=BDIM; B += di2) {
         d1 = fread_val(fi1);
         d2 = fread_val(fi2);       
         fwrite_val(fo,min(d1,d2));
      }
    }
  } else if (oper == OPER_XOR) {
    for (A=0; A<=ADIM; A += di1) {        
      for (B=0; B<=BDIM; B += di2) {
         d1 = fread_val(fi1);
         d2 = fread_val(fi2);       
         if ((d1 > 0.0 && d2 < 0.0) || (d1 < 0.0 && d2 > 0.0))
           fwrite_val(fo,1);
         else
           fwrite_val(fo,-1);
      }
    }
  } else if (oper == OPER_OR) {
    for (A=0; A<=ADIM; A += di1) {        
      for (B=0; B<=BDIM; B += di2) {
         d1 = fread_val(fi1);
         d2 = fread_val(fi2);       
         fwrite_val(fo,max(d1,d2));
      }
    }
  } else {
    fprintf(stderr,"error: wrong argument to oper()\n");
    exit(-1);
  }
  fclose(fi1);
  fclose(fi2);
  fclose(fo);
}

int tclResAND(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  if (argc != 4) {
    interp->result = "Usage: resand <.out file> <.out file> <new .out file>";
    return TCL_ERROR;
  }
  oper(argv[3],argv[1],argv[2],OPER_AND);
  return TCL_OK;
}

int tclResXOR(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  if (argc != 4) {
    interp->result = "Usage: resxor <.out file> <.out file> <new .out file>";
    return TCL_ERROR;
  }
  oper(argv[3],argv[1],argv[2],OPER_XOR);
  return TCL_OK;
}


int tclResOR(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  if (argc != 4) {
    interp->result = "Usage: resor <.out file> <.out file> <new .out file>";
    return TCL_ERROR;
  }
  oper(argv[3],argv[1],argv[2],OPER_OR);
  return TCL_OK;
}

/*
        z   baxis
        |  / 
        | / 
        |/_____y
       /|
      /a|
     /  |
    x   haxis

   The peptide plane fixed coordinate system is x,y,z with x along
   the N-H bond, and N in the origin. The helix axis is in the x-y plane
   ie. the peptide plane. The a value (alpha in the C-code) is the
   angle in degrees between the x and the helix axis. The B-field
   axis is described by the angles A and B ie. the polar angles.
   The angle between the haxis and the baxis is the one we want to
   describe by the lines in the output file. If it hits (less away
   than the deviation in degrees) one of the predefined angles we 
   draw a black dot.
   
   The angle beta is the second spherical angle for the helix-axis
*/

int tclResTiltLines(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  FILE *of;
  double alpha,beta,adev,bdev,alpha2,beta2;
  int i,A,B,alpha_ML_step,beta_ML_step,hit;
  char* outname;
  double  bitvalue[ADIM+1][BDIM+1];
  double ang[6],dot,haxis[6][4],baxis[4];
  double angles[256],a,b;
  double cosa,sina,sinb,minang,maxang,tmp;
  int nangles,argn;
  char **par;

  
  if (argc != 9) 
    return TclError(interp,"Usage: restiltlines <outfile> <alpha> <beta> <adev> <bdev> <astep> <bstep>  <list of angles 0-90 degrees>");

  argn=0;
  outname=argv[++argn];
  
  if (Tcl_GetDouble(interp,argv[++argn],&alpha) == TCL_ERROR) {
     sprintf(interp->result,"restiltlines: unable to make double out of '%s'\n",argv[argn]);
     return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp,argv[++argn],&beta) == TCL_ERROR) {
     sprintf(interp->result,"restiltlines: unable to make double out of '%s'\n",argv[argn]);
     return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp,argv[++argn],&adev) == TCL_ERROR) {
     sprintf(interp->result,"restiltlines: unable to make double out of '%s'\n",argv[argn]);
     return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp,argv[++argn],&bdev) == TCL_ERROR) {
     sprintf(interp->result,"restiltlines: unable to make double out of '%s'\n",argv[argn]);
     return TCL_ERROR;
  }
  if (Tcl_GetInt(interp,argv[++argn],&alpha_ML_step) == TCL_ERROR) {
     sprintf(interp->result,"restiltlines: unable to make integer out of '%s'\n",argv[argn]);
     return TCL_ERROR;
  }
  if (Tcl_GetInt(interp,argv[++argn],&beta_ML_step) == TCL_ERROR) {
     sprintf(interp->result,"restiltlines: unable to make integer out of '%s'\n",argv[argn]);
     return TCL_ERROR;
  }

  if (Tcl_SplitList(interp,argv[++argn],&nangles,&par) != TCL_OK) {
     sprintf(interp->result,"restiltlines: list is not formed correctly\n");
     return TCL_ERROR;
  }
  for (i=0;i<nangles;i++) {
    if (Tcl_GetDouble(interp,par[i],&angles[i]) != TCL_OK) {
      sprintf(interp->result,"restiltlines: unable to convert '%s' to a value in list element number %d\n",par[i],i+1);
      return TCL_ERROR;
    }
  }
  free(par);


  sinb=sin((beta-bdev)*M_PI/180.0);
  haxis[1][1]=cos((alpha-adev)*M_PI/180.0)*sinb;
  haxis[1][2]=sin((alpha-adev)*M_PI/180.0)*sinb;
  haxis[1][3]=cos((beta-bdev)*M_PI/180.0);

  sinb=sin((beta+bdev)*M_PI/180.0);
  haxis[2][1]=cos((alpha-adev)*M_PI/180.0)*sinb;
  haxis[2][2]=sin((alpha-adev)*M_PI/180.0)*sinb;
  haxis[2][3]=cos((beta+bdev)*M_PI/180.0);

  sinb=sin((beta+bdev)*M_PI/180.0);
  haxis[3][1]=cos((alpha+adev)*M_PI/180.0)*sinb;
  haxis[3][2]=sin((alpha+adev)*M_PI/180.0)*sinb;
  haxis[3][3]=cos((beta+bdev)*M_PI/180.0);

  sinb=sin((beta-bdev)*M_PI/180.0);
  haxis[4][1]=cos((alpha+adev)*M_PI/180.0)*sinb;
  haxis[4][2]=sin((alpha+adev)*M_PI/180.0)*sinb;
  haxis[4][3]=cos((beta-bdev)*M_PI/180.0);

  sinb=sin((beta)*M_PI/180.0);
  haxis[5][1]=cos((alpha)*M_PI/180.0)*sinb;
  haxis[5][2]=sin((alpha)*M_PI/180.0)*sinb;
  haxis[5][3]=cos((beta)*M_PI/180.0);

  alpha2= 180 + alpha;
  if (alpha2 > 360.0) alpha2 -= 360.0;
  beta2= 180 - beta;
  if (beta2 < 0.0) beta2 += 180.0;
  
  for (A=0; A<ADIM; A += alpha_ML_step) {
    a=A/(double)STEP;
    cosa=cos(a*M_PI/180.0);
    sina=sin(a*M_PI/180.0);
    for (B=0; B<BDIM; B += beta_ML_step) {
       maxang= -1.0e6;
       minang= 1.0e6;
       b=B/(double)STEP;
       sinb=sin(b*M_PI/180.0);
       baxis[1]=sinb*cosa;
       baxis[2]=sinb*sina;
       baxis[3]=cos(b*M_PI/180.0);

       for (i=1;i<=5;i++) {
         dot= haxis[i][1]*baxis[1]+haxis[i][2]*baxis[2]+haxis[i][3]*baxis[3];
         ang[i]=acos(dot)*180.0/M_PI;         

         if (ang[i] < minang) minang=ang[i];
         if (ang[i] > maxang) maxang=ang[i];
       } 
       
       if (maxang >= 90 && minang >= 90) {
         tmp=minang;
         minang=180-maxang;
         maxang=180-tmp;
         
       } else if (minang < 90 && maxang > 90) {
         maxang=180-maxang;
         if (maxang < minang) minang=maxang;
         maxang=90;
       }

/*       k=0;      
       if (a < alpha + adev && a > alpha - adev) {
          tmp= ang[5]-bdev;
          if (tmp < minang) minang=tmp;
          k++;
       }
       if (b < beta  + bdev && b > beta - bdev) {
          tmp= ang[5]-adev;
          if (tmp < minang) minang=tmp;
          k++;
       }
       if (k == 2) 
         minang=0;

       k=0;      
       if ( a < alpha2 + adev && a > alpha2 - adev) {
          tmp= ang[5]-bdev;
          if (tmp < minang) minang=tmp;
          k++;
       }
       if (b < beta2  + bdev && b > beta2 - bdev) {
          tmp= ang[5]-adev;
          if (tmp < minang) minang=tmp;
          k++;
       }
       if (k == 2) 
         minang=0;
*/


       hit=-1;

       for (i=0;i<nangles;i++) {
         if (angles[i] == 0) {
           if ( (a < alpha  + adev && a > alpha  - adev && 
                 b < beta   + bdev && b > beta   - bdev) ||
                (a < alpha2 + adev && a > alpha2 - adev && 
                 b < beta2  + bdev && b > beta2  - bdev) ) {
             hit=1;
             break;
           }
         } else
         if (angles[i] <= maxang && angles[i] >= minang) {
           hit=1; 
           break; 
         }
       }
       bitvalue[A][B]=hit;
    }
  }
  for (A=0; A<=ADIM; A += alpha_ML_step) {
    bitvalue[A][BDIM]=bitvalue[A][0];
  }
  for (B=0; B<=BDIM; B += beta_ML_step){
    bitvalue[ADIM][B]=bitvalue[0][B];
  }
  
  of=fopen(outname,"w");
  if (!of) {
    sprintf(interp->result,"error: unable to create file '%s'",outname);
    return TCL_ERROR;
  }
  fwrite_val(of,alpha_ML_step);
  fwrite_val(of,beta_ML_step);
  
  for (A=0; A<=ADIM; A += alpha_ML_step) {        
    for (B=0; B<=BDIM; B += beta_ML_step) {
       fwrite_val(of,bitvalue[A][B]);
    }
  }
  fclose(of);
  return TCL_OK;
}

double dd(double a,double b,double a2,double b2)
{
  double x,y;
  x=(a-a2);
  y=(b-b2);
  return sqrt(x*x+y*y);
}


double dist(double a2,double b2,double a,double b)
{
  double x,y;
  x=dd(a,b,a2,b2);
  a += 180;
  if (a > 360.0) a -= 360;
  if (a < 0.0) a += 360;

  b = 180 - b;
  if (b > 180.0) b -= 180;
  if (b < 0.0) b += 180;
    
  y=dd(a,b,a2,b2);
  if (x < y) return x;
  return y;
}


void tclcmd_restrict(Tcl_Interp* interp) {
   Tcl_CreateCommand(interp,"resshift",tclResShift,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
   Tcl_CreateCommand(interp,"resdipole",tclResDipole,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
   Tcl_CreateCommand(interp,"resand",tclResAND,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
   Tcl_CreateCommand(interp,"resxor",tclResXOR,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
   Tcl_CreateCommand(interp,"resor",tclResOR,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
   Tcl_CreateCommand(interp,"restiltlines",tclResTiltLines,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
   Tcl_CreateCommand(interp,"resmakeps",tclResMakePs,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
   Tcl_CreateCommand(interp,"resmakeps4",tclResMakePs4,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
}












