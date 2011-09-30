/*
    Fast routines for CSA calculation
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
    

    This performs a fast CSA calculation based on the COMPUTE algoritm.
    They are stand alone Tcl routines and allows for a simulation like
    this:
    
    proc main {} {
        global par stop

        csainit 128 128 500 rep168
        csaspec -4214 0.25

        set lst {}
        foreach i [csareturnint -4133  31.88] {
           lappend lst [concat [lrange $i 1 2] 100 0]
        }
        set f [fcreate -np 2048 -sw 20000 -type spe]
        faddpeaks $f 0 $lst
        fsave $f $par(name).spe
        funload $f
    }


*/

#include <stdlib.h>
#include <string.h>
#ifdef __APPLE__
#  include <malloc/malloc.h>
#else
#  include <malloc.h>
#endif
#include <ctype.h>
#include <tcl.h>
#include "matrix_new.h"

#ifndef PI2
#define PI2 (M_PI*2.0)
#endif

#ifndef D2R
#define D2R (M_PI/180.0)
#endif


void csa_gammid(complx* A22,double* fidsum,double weight,int first)
{
  int ct,j,k,na,nsw;
  complx c;
  static complx* qt = NULL;
  static complx** fexp = NULL;
  double FAK;
  complx* f;
  int na1,na2;

  nsw=LEN(fidsum);
  na=LEN(A22);
     
  if (first) {
    if (qt && LEN(qt) != na) {
      free_complx_vector(qt);
      qt=NULL;
    }
    if (!qt)
      qt=complx_vector(na);

    if (fexp && ROWS(fexp) != na) {
      free_complx_matrix(fexp);
      fexp=NULL;
    }
    if (!fexp)
      fexp=complx_matrix(na,na);

    if ((na%2)==0) {
      for (k=1;k<=na;k++) {
        for (j=0;j<na;j++) {
          fexp[k][j+1]= Cexpi(-(double)j*(PI2*(double)(k-na/2.0)/(double)na));
        }
      }
    } else {
      for (k=1;k<=na;k++) {
        for (j=0;j<na;j++) {
          fexp[k][j+1]= Cexpi(-(double)j*(PI2*(double)(k-(na-1)/2-1)/(double)na));
        }
      }
    }
  }
  qt[0].re=1.0;
  qt[0].im=0.0;
  for (ct=1;ct<=na;ct++) {
      qt[ct] = Cmul(A22[ct],A22[ct]);
  }
  FAK=weight/(double)(na*na);
  /* this is 66% of the full computation time */
  if (nsw > na) {
    fprintf(stderr,"error: csainit: na must be larger than or equal to nsw\n");
    exit(1);
  } else if (nsw == na) {
    na1=1;
    na2=na;
  } else  {
    na1= (na - nsw)/2+1;
    na2= na - na1;
  }
  for (k=na1;k<=na2;k++) {
    c.re=0.0;
    c.im=0.0;
    f=fexp[k];
	  for (j=0;j<na;j++) {
      /* c = Cadd(c,Cmul(qt[j],f[j+1])); */
       complx *a=&qt[j];
       complx *b=&f[j+1];
       c.re += a->re * b->re - a->im * b->im;
       c.im += a->im * b->re + a->re * b->im;
    }
    fidsum[k-na1 + 1] +=  Cnorm(c)*FAK;
	}
}

void csa_cartrot(mv_double* Urot,mv_double* U,double alpha,double beta,double gamma)
{
  double ca,cb,cg,sa,sb,sg;
  int m,n,j,k;
  double R[3][3];
 
  alpha *= D2R;
  beta *= D2R;
  gamma *= D2R;
  ca=cos(alpha);
  cb=cos(beta);
  cg=cos(gamma);
  sa=sin(alpha);
  sb=sin(beta);
  sg=sin(gamma);
  R[0][0]= ca*cb*cg-sa*sg;
  R[0][1]= sa*cb*cg+ca*sg;
  R[0][2]= -sb*cg;
  R[1][0]= -ca*cb*sg-sa*cg;
  R[1][1]= -sa*cb*sg+ca*cg;
  R[1][2]= sb*sg;
  R[2][0]= ca*sb;
  R[2][1]= sa*sb;
  R[2][2]= cb;

  memset(Urot->data,0,9*sizeof(double));
  for (m=0;m<3;m++) {
    for (n=0;n<3;n++) {
      for (j=0;j<3;j++) {
        for (k=0;k<3;k++) {
          Urot->data[m+n*3] += U->data[j+k*3]*R[m][j]*R[n][k];
        }
      } 
    }
  } 
}

void csa_makesigma(mv_double* s,double iso,double aniso,double eta)
{
 
/* defined on the shift/deshielding axis
  (same as publications and spectrometer software)
  iso=(xx+yy+zz)/3
  |zz-iso| >= |xx-iso| >= |yy-iso|
  aniso = zz - iso
  eta = (yy-xx)/aniso
*/
   iso *= PI2;
   aniso *= PI2;
   s->data[0]=iso+aniso/2.0*(1.0+eta);
   s->data[0+1*3]=0.0;
   s->data[0+2*3]=0.0;
   s->data[1]=0.0;
   s->data[1+1*3]=iso+aniso/2.0*(1.0-eta);
   s->data[1+2*3]=0.0;
   s->data[2]=0.0;
   s->data[2+1*3]=0.0;
   s->data[2+2*3]=iso-aniso;
}

#define K  1
#define C1 2
#define C2 3
#define S1 4
#define S2 5

/*sin(2 brl)*/
#define S2MA 0.942809041585

/* sin(brl)*sin(brl) */
#define SMA2 0.666666666666666666666666666

typedef double CONSTS[6];

void csa_makeconst(double* cvec,mv_double* s)
{
  cvec[K ]=SMA2*(0.5*(s->data[0+0*3]+s->data[1+1*3])-s->data[2+2*3])+s->data[2+2*3];
  cvec[C1]=S2MA*s->data[0+2*3];
  cvec[C2]=SMA2*0.5*(s->data[0]-s->data[1+1*3]);
  cvec[S1]=S2MA*s->data[1+2*3];
  cvec[S2]=SMA2*s->data[0+1*3];
}


void csa_integcalc(double* f,double t,double wr)
{
  double ph;

  ph= -wr*t;
  f[1]= t;
  f[2]= sin(ph)/wr;
  f[3]= 0.5*(sin(2.0*ph))/wr;
  f[4]= (cos(ph)-1.0)/wr;
  f[5]= 0.5*(cos(2.0*ph)-1.0)/wr;

}

double csa_integmul(double* f,double* c)
{
  return c[K]*f[1] + c[C1]*f[2] + c[C2]*f[3] - c[S1]*f[4] - c[S2]*f[5];
}


int loaded=0;
int nsw,na,ncr;
double crystintensity;
double srate,wr,sw,dt;
char crystname[256];

double *fid = NULL;
mv_double *omega = NULL;
complx* A22 = NULL;
CONSTS* f=NULL,c;
/*double **sc = NULL,**sr = NULL; */
mv_double *sc=NULL, *sr=NULL;

mv_double * read_crystfile_byname(char* crystname);

/*
   csainit <gamma-angles> <n-spec-width> <spin-rate> <crystal-file>
   
   gamma-angles is the number of Euler gamma angles in the powder averaging
   gamma-angles must be larger than or equal to n-spec-width
   n-spec-width is the spectral width in multiples of the spin-rate
   spin-rate is the sample spinning speed given in Hz 
   crystal-file is a powder averagning file containing alpha and beta crystallites

*/
int tclCsaInit(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int i,n,a;
  double *dw;

  if (argc != 5) {
    interp->result="usage: csainit <gamma-angles> <n-spec-width> <spin-rate> <crystal-file>";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp,argv[1],&na) == TCL_ERROR) {
     interp->result="csainit: argument 1 must be integer <gamma-angles>";
     return TCL_ERROR;
  }
  if (Tcl_GetInt(interp,argv[2],&nsw) == TCL_ERROR) {
     interp->result="csainit: argument 2 must be integer <n-spec-width>";
     return TCL_ERROR;
  }
     if (Tcl_GetDouble(interp,argv[3],&srate) == TCL_ERROR) {
     interp->result="csanit: argument 3 must be double <spin-rate>";
     return TCL_ERROR;
  }
  if (nsw > na) {
    interp->result="csainit: <gamma-angles> must be larger than or equal to <n-spec-width>\n";
    return TCL_ERROR;
  }

  strcpy(crystname,argv[4]);

  if (omega) double_matrix_free(omega);
  
  omega = read_crystfile_byname(crystname);
  ncr = omega->row;
  if (!sc) sc=double_matrix_alloc(3,3);
  if (!sr) sr=double_matrix_alloc(3,3);

  crystintensity=0.0;
  /*n=LEN(omega);
  for (i=1;i<=n;i++) {
    crystintensity += omega[i].weight;
  }*/
  dw = omega->data + 3*ncr;
  for (i=0; i<ncr; i++) {
     crystintensity += *dw;
     dw++;
  }
  
  if (fid && LEN(fid) != nsw) {
    free_double_vector(fid);
    fid=NULL;
  }
  if (!fid) fid=double_vector(nsw);

  loaded=1;

  wr=srate*M_PI*2.0;
  sw=srate*na;
  dt=1.0/(srate*na);

  if (A22 && LEN(A22) != na) {
    free_complx_vector(A22);
    A22=NULL;
  }
  if (!A22) A22=complx_vector(na);

  if (f) free(f);
  f=malloc(sizeof(CONSTS)*(na+1));
  if (!f) {
     sprintf(interp->result,"error: unable to allocate %d times the CONSTS structure\n",na);
     return TCL_ERROR;
  }
  for (a=1;a<=na;a++) {
   double t=a*dt;
   csa_integcalc(f[a],t,wr);
 }

  return TCL_OK; 
}

/*
  csasetspinrate <spin-rate>
  
  re-sets the spin-rate after 'csainit' is called.
*/
int tclCsaSetSpinRate(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int a;

  if (argc != 2) {
    interp->result="usage: csasetspinrate <spin-rate> ";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp,argv[1],&srate) == TCL_ERROR) {
    interp->result="argument 1 must be double <spin-rate>";
    return TCL_ERROR;
  }
  wr=srate*M_PI*2.0;
  sw=srate*na;
  dt=1.0/(srate*na);
  for (a=1;a<=na;a++) {
    double t=a*dt;
    csa_integcalc(f[a],t,wr);
  }
  return TCL_OK; 
}


/*
   csaspec <aniso> <eta> ?<alpha> <beta> <gamma>?
   
   calculates the intensities in a CSA powder pattern given the parameters
   initiated with the 'csainit' procedure and the anisotropy (aniso)
   asymmetry (eta). The alpha, beta and gamma Euler angles gives the
   tensor orientation, and is only used if one wants to check the 
   powder averagning for rotation equivalence.
   The result is returned with a call to 'csareturnint'
*/
int tclCsaSpec(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int icr,a,i,n,j;
  double aniso,eta;

  if (!loaded) {
    interp->result="csainit must be called before csaspec";
    return TCL_ERROR;
  }
  if (argc != 3 && argc != 6) {
    interp->result="usage: csaspec <aniso> <eta> ?<alpha> <beta> <gamma>?";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp,argv[1],&aniso) == TCL_ERROR) {
     interp->result="argument 1 must be double <aniso>";
     return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp,argv[2],&eta) == TCL_ERROR) {
     interp->result="argument 2 must be double <eta>";
     return TCL_ERROR;
  }
  csa_makesigma(sc,0.0,aniso,eta);

  if (argc == 6) {
    double a,b,g;
    if (Tcl_GetDouble(interp,argv[3],&a) == TCL_ERROR) {
       interp->result="argument 3 must be double <alpha>";
       return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp,argv[4],&b) == TCL_ERROR) {
       interp->result="argument 4 must be double <beta>";
       return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp,argv[5],&g) == TCL_ERROR) {
       interp->result="argument 5 must be double <gamma>";
       return TCL_ERROR;
    }
    csa_cartrot(sr,sc,a,b,g);
    /*for (i=1;i<=3;i++)
      for (j=1;j<=3;j++)
         sc[i][j]=sr[i][j]; */
    memcpy(sc->data,sr->data,9*sizeof(double));
  }

  csa_makeconst(c,sc);
  n=LEN(fid);
  for (i=1;i<=n;i++)
     fid[i]=0.0;

  for (icr=1;icr<=ncr;icr++) {
      double alpha = double_get_elem(omega,icr,1);
      double beta =  double_get_elem(omega,icr,2);
      double weight = double_get_elem(omega,icr,4);
      csa_cartrot(sr,sc,alpha,beta,0.0);
      csa_makeconst((double*)c,sr);
      for (a=1;a<=na;a++) {
        double omega = csa_integmul(f[a],c);
        A22[a]=Cexpi(omega*0.5);
      }
      csa_gammid(A22,fid,weight,(icr == 1));
  }
  n=LEN(fid);
  for (i=1;i<=n;i++)
     fid[i] /= crystintensity;

  return TCL_OK;
}

typedef struct _peak {
  int ssb;
  double weight;
  double freq;
  double intensity;
  double lb;
  double lbweight;
} peak;

int npeaks=0;
peak peaks[200];

/*
    csareturnint <iso> ?<scaling=1>?
    
    Returns the calculated CSA spectrum spinning sideband intensities in a list
    of lists each containing a frequency and intensity. Iso is the isotropic chemical
    shielding ie. the position of the central peak. The scaling parameter
    is multiplied to the intensities.
*/
int tclCsaReturnInt(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  char buf[256];
  double iso,scaling;
  int i,n,n2,ssb;

  if (argc != 3 && argc != 2) {
    interp->result="usage: csareturnint <iso> ?<scaling=1>? (returns the calculated spectrum intensities)";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp,argv[1],&iso) == TCL_ERROR) {
    interp->result="argument 1 must be double <iso>";
    return TCL_ERROR;
  }
  if (argc == 3) {
    if (Tcl_GetDouble(interp,argv[2],&scaling) == TCL_ERROR) {
      interp->result="argument 2 must be double <scaling>";
      return TCL_ERROR;
    }
  } else {
    scaling=1.0;
  }

  Tcl_ResetResult(interp);
  n=LEN(fid);
  n2=(n-n%2)/2;

  for (i=1;i<=n;i++) {
    if (fid[i] == 0.0) continue;
    ssb=i-n2-n%2;
    sprintf(buf,"%d %g %g",ssb,srate*(double)ssb+iso,fid[i]*scaling);
    Tcl_AppendElement(interp,buf);
  }
  return TCL_OK;
}

void tclcmd_csa(Tcl_Interp* interp)
{
  Tcl_CreateCommand(interp,"csainit",tclCsaInit,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"csasetspinrate",tclCsaSetSpinRate,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"csaspec",tclCsaSpec,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"csareturnint",tclCsaReturnInt,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
}

