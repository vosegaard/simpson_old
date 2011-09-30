/*
    Reading and setup of the spinsystem
    Copyright (C) 1999 Mads Bak, Jimmy T. Rasmussen
    2009 ZT modification with new matrix types
    
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
    
    Reads the spin system from the 'spinsys' section in the input
    file and creates the Hamiltonian for the spin system.
    
     Called from sim.c the setup and performs the simulation
*/

#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <tcl.h>
#include "matrix_new.h"
#include "spinsys.h"
#include "ham.h"
#include "sim.h"
#include "defs.h"
#include "cm_new.h"
#include "defs_blas_lapack.h"

/* Makes the multiplication and destroys the memory that's not used.
   D for dynamic allocation.
*/

complx* Dtensor2(double deltazz,double eta)
{
  complx* v;
  
  v=(complx*)malloc(5*sizeof(complx));
  v[0]=v[4]= Complx(-eta*deltazz/2.0,0.0);
  v[1]=v[3]= Complx(0.0,0.0);
  v[2]=      Complx(deltazz*sqrt(3.0/2.0),0.0);
  return v;
}


complx * Drot(complx* R,double alpha,double beta, double gamma)
{
  static mv_complx *wsp=NULL;
  complx *v;

  v = (complx*)malloc(5*sizeof(complx));
  wsp=cmv_static(wsp,5,5);
  wigner2(wsp->data,alpha,beta,gamma);
  wig2rot(v,R,wsp->data);
  free((char*)R);
  return v;
}


mv_complx * Dmul(mv_complx * m1,mv_complx * m2)
{
   mv_complx * m;
   m=complx_matrix_alloc(m1->row,m1->col);
   cm_mul(m,m1,m2);
   complx_matrix_free(m1);
   complx_matrix_free(m2);
   return m;
}

mv_complx * Dadd(mv_complx * m1,mv_complx * m2)
{
   cmv_addto(m1,m2);
   complx_matrix_free(m2);
   return m1;
}

mv_complx * Dsub(mv_complx * m1,mv_complx * m2)
{
   cmv_subfrom(m1,m2);
   complx_matrix_free(m2);
   return m1;
}

mv_complx * Dmulc(mv_complx * m,double re,double im)
{
   cmv_mulc(m,Complx(re,im));
   return m;
}

mv_complx * Dmuld(mv_complx * m,double re)
{
   cmv_muld(m,re);
   return m;
}

mv_complx * IxIx(SpinSys* ss,int n1,int n2) { return Dmul(Ix(ss,n1),Ix(ss,n2)); }
mv_complx * IyIy(SpinSys* ss,int n1,int n2) { return Dmul(Iy(ss,n1),Iy(ss,n2)); }
mv_complx * IzIz(SpinSys* ss,int n1,int n2) { 
  return Dmul(Iz(ss,n1),Iz(ss,n2)); 
}

mv_complx * II(SpinSys* ss,int n1,int n2)
{ return Dadd(Dadd(IxIx(ss,n1,n2),IyIy(ss,n1,n2)),IzIz(ss,n1,n2)); }

mv_complx * T20(SpinSys* ss,int n1,int n2)
{
   return Dmuld(
      Dsub( Dmuld(IzIz(ss,n1,n2),2.0), Dadd( IxIx(ss,n1,n2), IyIy(ss,n1,n2) ) ),
      1.0/sqrt(6.0));
}

void readsys(Tcl_Interp* interp,Sim* s)
{
  SpinSys* ss;
  Hamilton* h;
  mv_complx *T,*Ta,*Tb;
  double *ReT, *ReTa, *ReTb;
  complx *R,*R2;
  double C,v0,wq,v[100];
  int N;
  char name[256];
  int ver,nelem,nval[100];
  char **p0,**elem,**val[100];
  int  nspins=0;
  char **spins=NULL;
  char **cchan=NULL;
  int i,j,l,ok,n1,n2,kk;
  char *buf;
  char *liq;
  int is_liq, n_mix=0, mix[100], qn, dn;
  /* ZT: tables for interaction parameters. Used in mixing terms and could
   *     possibly be in relaxation - neds to be linked with Relax structure...
   */
  double **shifts, **quadrupoles, ***dipoles;

  ss=s->ss;
  h=s->H;
  ss_initialize(ss);

  ver = (verbose & VERBOSE_SPINSYS);
  
  if (ver) printf("Reading spinsystem:\n");

  /* ZT: handling of 'liquid on' in par section - only isotropic interactions 
   *     should be read here
   */
  if ((liq=Tcl_GetVar2(interp,"par","liquid",TCL_GLOBAL_ONLY)) == NULL) {
     is_liq = 0;
  } else {
     if (!strncmp(liq,"on",2)) {
        is_liq = 1;
     } else {
        is_liq = 0;
     }
  }
  /* printf("liquid is %s\n",(is_liq ? "on" : "off")); */
   
  if (Tcl_Eval(interp,"array names spinsysres") != TCL_OK) {
    fprintf(stderr,"%s\n",interp->result);
    exit(1);
  }
  /* RA+ZT: this is always enough for buf! 
     Reuse of buf is subsets of this content 
  */
  buf = (char *) malloc(sizeof(char) * (strlen(interp->result)+1));
  strcpy(buf,interp->result);

  if (Tcl_SplitList(interp,buf,&nelem,&p0) != TCL_OK) {
    fprintf(stderr,"%s\n",interp->result);
    exit(1);
  }
  elem = p0 - 1;
  for (i=1;i<=nelem;i++) {     
    char* pval=Tcl_GetVar2(interp,"spinsysres",elem[i],0);
    if (ver) {
      strcpy(buf,elem[i]);
      l=strlen(buf)-1;
      while (buf[l] == 'X') l--;
      buf[l+1]=0;
      printf("  %s ",buf);
    }
    if (!pval) {
      fprintf(stderr,"%s\n",interp->result);
      exit(1);
    }
    if (Tcl_SplitList(interp,pval,&nval[i],&p0) != TCL_OK) {
      fprintf(stderr,"%s\n",interp->result);
      exit(1);
    }
    val[i] = p0 - 1;
    if (ver) {
      for (j=1;j<=nval[i];j++) {
        printf(" %s",val[i][j]);
      }
      printf("\n");
    }
  }

  ok=0;
  for (i=1;i<=nelem;i++) {
    if (!strcmp(elem[i],"nuclei")) {
      ok=1;
      nspins=nval[i];
      spins=val[i];
      break;
    }
  }
  if (!ok) {
    fprintf(stderr,"error: the spinsystem must contain a 'nuclei' field\n");
    exit(1);
  }

  ok=0;
  for (i=1;i<=nelem;i++) {
    if (!strcmp(elem[i],"channels")) {
      ok=1;
      ss->nchan=nval[i];
      cchan=val[i];
      break;
    }
  }
  if (!ok) {
    fprintf(stderr,"error: the spinsystem must contain a 'channels' field\n");
    exit(1);
  }

  if (ver) printf( "Number of spins : %d\n",nspins);

  
  for (i=1;i<=nspins;i++) {
    if (ver) printf("  nucleus %d  : '%s'\n",i,spins[i]);
    ss_addspin(ss,spins[i]);
  }

  for (i=1;i<=ss->nchan;i++) {
    strcpy(ss->channames[i],cchan[i]);
    ss->nchanelem[i]=0;
    for (j=1;j<=nspins;j++) {
       if (!strcmp(spins[j],ss->channames[i])) {
         ss->nchanelem[i]++;
         ss->chan[i][ss->nchanelem[i]]=j;
       }
    }
  }
  if (ver) {
    printf("Number of channels : %d\n",ss->nchan);
    for (i=1;i<=ss->nchan;i++) {
      printf( "  channel %d  : '%s' contain",i,ss->channames[i]);
      for (j=1;j<=ss->nchanelem[i];j++) {
         printf(" nucleus(%d,%s)",ss->chan[i][j],spins[ss->chan[i][j]]);
      }
      printf("\n");
    }
  }


  /* ZT: here we are about to read spin system parametrs
   *     -> initialize interaction maps
   */
   shifts = (double**)malloc((nspins+1)*sizeof(double*));
   quadrupoles = (double**)malloc((nspins+1)*sizeof(double*));
   dipoles = (double***)malloc((nspins+1)*sizeof(double**));
   dipoles[1] = (double**)malloc((nspins*nspins+1)*sizeof(double*));
      for (i=2; i<=nspins; i++) dipoles[i]=dipoles[i-1]+nspins;
   for (i=1; i<=nspins; i++) {
      shifts[i]=NULL;
      quadrupoles[i]=NULL;
      for (j=1; j<=nspins; j++) {
	 dipoles[i][j] = NULL;
      }
   }

  N=ss->matdim;
  ham_initialize(h,N);

  h->isdiag=1;    
  for (i=1;i<=nelem;i++) {
     if (!strcmp(elem[i],"nuclei")) continue;
     if (!strcmp(elem[i],"channels")) continue;
     l=strlen(elem[i])-1;
     while (elem[i][l] == 'X') 
       elem[i][l--]=0;

     for (j=1;j<=nval[i];j++) {
       v[j] = atof(val[i][j]);                
       //printf("v[%i] is %d \n", j, v[j]);
       errno = 0;
       if (errno) {
	 sprintf(buf,"failure in spinsys line %d word %d\n",i,j+1);
	 perror(buf);
	 exit(-1);
       }
     }
     if (!strcmp(elem[i],"shift")) {
       n1=(int)v[1];
/*
       if (v[4] < 0.0) v[4]=0.0; else
       if (v[4] > 1.0) v[4]=1.0;
*/
       sprintf(name,"shift_%d",n1);
//printf("   readsys doing %s\n",name);
       if (ham_exists(h,name)) {
          fprintf(stderr,"error: spinsys: interaction '%s' does already exists\n",name);
          exit(1);
       }
/* No need to bail out when shift is zero.
       if (v[2] == 0.0 && v[3] == 0.0) {
         fprintf(stderr,"error: spinsys: interaction '%s' has zero coupling strength\n",name);
         exit(1);
       }
*/

       /* ZT: modification for liquid on */
       if ( (ver) && (!is_liq) ) {
         /* solid state situation */
         printf( "Chemical shift on nucleus %d\n",n1);
         printf( "  isotropic shift        : %g Hz\n",v[2]);
         printf( "  anisotropic shift      : %g Hz\n",v[3]);
         printf( "  assymmetry parameter   : %g\n",v[4]);
         printf( "  euler angles of tensor :  (%g,%g,%g) degrees\n",v[5],v[6],v[7]); 
       }
       if ( (ver) && (is_liq) ) {
         /* liquid state situation */
         printf( "Chemical shift on nucleus %d\n",n1);
         printf( "  isotropic shift        : %g Hz\n",v[2]);
         printf( "  anisotropic part ignored for spectra calculations in liquid\n"); 
       }
       v[2] *= (2.0*M_PI);
       v[3] *= (2.0*M_PI);

       if (ss_gamma(ss,n1) > 0) {
         v[2] = -v[2];
         v[3] = -v[3];
       }
       
       T=Iz(ss,n1);
       kk = 0;
       if (v[2] != 0.0){
         /* ham_add_static(h,Complx(v[2],0.0),m_dup(T),name); */
         ReT = get_real_diag(T);
	 ham_add_static(h,v[2],ReT,name,1);
       }
       if ( (v[3] != 0.0) && (!is_liq) ) {       
         C= v[3]*sqrt(2.0/3.0);
         R=Dtensor2(C,v[4]);
         R = Drot(R,v[5],v[6],v[7]);       
         ReT = get_real_diag(T);
         ham_add(h,R,ReT,name,1);
       }
       complx_matrix_free(T);
              
       /* ZT: store original values in shifts */
       shifts[n1] = (double*)malloc(6*sizeof(double));
       for (kk=0; kk<6; kk++) shifts[n1][kk] = v[kk+2];
//printf("               ------> done\n");
       
     } else if (!strcmp(elem[i],"dipole")) {
       int bsign;
       n1=(int)v[1];
       n2=(int)v[2];
       
       /* ZT: modification for liquid on */
       if (is_liq) {
          if (ver) printf("Dipolar coupling between nucleus %d and %d ignored for specra calculation in liquid\n", n1, n2);
       } else {
         sprintf(name,"dipole_%d_%d",n1,n2);
//printf("   readsys doing %s\n",name);
         if (ham_exists(h,name)) {
            fprintf(stderr,"error: spinsys: interaction '%s' does already exists\n",name);
            exit(1);
         }
/*
       if (v[3] == 0.0) {
         fprintf(stderr,"error: spinsys: interaction '%s' has zero coupling strength\n",name);
         exit(1);
       }
*/
         if (ss_issame(ss,n1,n2)) {
           if (ver) printf( "Homonuclear ");
           T=T20(ss,n1,n2);
           h->isdiag = kk = 0;
	   ReT = get_real(T);   
         } 

	 else {
           if (ver) printf( "Heteronuclear ");
	   T=Dmuld( IzIz(ss,n1,n2), 2.0/sqrt(6.0) );
	   kk = 1;
	   ReT = get_real_diag(T);
         }
         if (ver) {
           printf( "dipolar coupling between nucleus %d and %d\n",n1,n2);
           printf( "  dipolar coupling       : %g Hz\n",v[3]);
           printf( "  euler angles of tensor :  (%g,%g,%g) degrees\n",v[4],v[5],v[6]); 
         }
         if (s->dipole_check) {
           bsign = ( -ss_gamma(ss,n1)*ss_gamma(ss,n2) > 0 ? 1 : -1);
           if (v[3] > 0.0 && bsign < 0) {
             fprintf(stderr,"error: the dipolar coupling between nucleus %d and %d must be negative to comply\n"
                            "       to the conventions used in this program. Set 'dipole_check' to 'false' to\n"
                            "       override this sign check.\n",n1,n2);
             exit(1);
           } else if (v[3] < 0.0 && bsign > 0) {
             fprintf(stderr,"error: the dipolar coupling between nucleus %d and %d must be positive to comply\n"
                            "       to the conventions used in this program. Set 'dipole_check' to 'false' to\n"
                            "       override this sign check.\n",n1,n2);
             exit(1);
           }
         }
         v[3] *= (2.0*M_PI);
         R=Dtensor2(2.0*v[3],0.0);
         R = Drot(R,v[4],v[5],v[6]);

         ham_add(h,R,ReT,name,kk);
	 complx_matrix_free(T);
       }

       /* ZT: store original values in dipoles */
       if (n1>n2) {
	    int dum;
	    dum = n2;
	    n2 = n1;
	    n1 = dum;
       }
       dipoles[n1][n2] = (double*)malloc(4*sizeof(double));
       for (kk=0; kk<4; kk++) dipoles[n1][n2][kk] = v[kk+3];
//printf("               ------> done\n");

     } else if (!strcmp(elem[i],"quadrupole")) {
       n1=(int)v[1];
       
       /* ZT: modification for liquid on */
       if (is_liq) {
          if (ver) printf("Quadrupolar coupling on nucleus %d ignored for specra calculation in liquid\n", n1);
       } else {
/*
       if (v[4] < 0.0) v[4]=0.0; else
       if (v[4] > 1.0) v[4]=1.0;
*/
         sprintf(name,"quadrupole_%d",n1);
         if (ham_exists(h,name)) {
            fprintf(stderr,"error: spinsys: interaction '%s' does already exists\n",name);
            exit(1);
         }
/*
       if (v[3] == 0.0) {
         fprintf(stderr,"error: spinsys: interaction '%s' has zero coupling strength\n",name);
         exit(1);
       }
*/
         if (v[2] == 2) {
           if (ver) {
             printf( "First and second order quadrupolar coupling on nucleus %d\n",n1);
             printf( "  quadrupolar constant         :  %g Hz\n",v[3]);
             printf( "  quadrupolar assymetry        :  %g\n",v[4]);
             printf( "  euler angles of tensor       :  (%g,%g,%g) degrees\n",v[5],v[6],v[7]); 
             printf( "  proton frequency (with sign) :  %g Hz\n",s->specfreq);
           }
           v[3] *= (2.0*M_PI);
           T=T20(ss,n1,n1);
           wq=v[3]/(4.0*ss_qn(ss,n1)*(2.0*ss_qn(ss,n1)-1));
           /* Iz ( 2 I^2 - 2 Iz^2 - 1 )*/
           Ta=Dmul(Dsub(Dsub(Dmuld(II(ss,n1,n1),2.0),Dmuld(IzIz(ss,n1,n1),2.0)),Ie(ss)),Iz(ss,n1));
           /* 1/2 Iz ( 4 I^2 - 8 Iz^2 - 1 )*/
           Tb=Dmuld(Dmul(Dsub(Dsub(Dmuld(II(ss,n1,n1),4.0),Dmuld(IzIz(ss,n1,n1),8.0)),Ie(ss)),Iz(ss,n1)),0.5);

           v0 = ss_gamma(ss,n1)*s->specfreq/ss_gamma1H()*2.0*M_PI;
           C= 2.0*wq*wq/v0;
           R=Dtensor2(1.0,v[4]);
           cmv_muld(T,2.0*wq);
           cmv_muld(Ta,C);
           cmv_muld(Tb,C);
           R = Drot(R,v[5],v[6],v[7]);
	   ReT = get_real_diag(T);
	   ReTa = get_real_diag(Ta);
	   ReTb = get_real_diag(Tb);
           ham_add_Q2(h,R,ReT,ReTa,ReTb,name);
	   complx_matrix_free(T);
           complx_matrix_free(Ta);
	   complx_matrix_free(Tb);

         } else if (v[2] == 1) {       
           if (ver) {
             printf( "First order quadrupolar coupling on nucleus %d\n",n1);
             printf( "  quadrupolar constant   :  %g Hz\n",v[3]);
             printf( "  quadrupolar assymetry  :  %g\n",v[4]);
             printf( "  euler angles of tensor :  (%g,%g,%g) degrees\n",v[5],v[6],v[7]); 
           }
           v[3] *= (2.0*M_PI);
           T=T20(ss,n1,n1);
           wq=v[3]/(4.0*ss_qn(ss,n1)*(2.0*ss_qn(ss,n1)-1));
           R=Dtensor2(2.0*wq,v[4]);
           R = Drot(R,v[5],v[6],v[7]);
	   ReT = get_real_diag(T);
           ham_add(h,R,ReT,name,1);
	   complx_matrix_free(T);

         } else {
           fprintf(stderr,"error: spinsys: argument 2 to 'quadrupole' "
                           "must be 1 or 2 (i.e. first and second order interaction)");
           exit(1);
         }
       }
       /* ZT: store original values in quadrupoles */
       quadrupoles[n1] = (double*)malloc(5*sizeof(double));
       for (kk=0; kk<5; kk++) quadrupoles[n1][kk] = v[kk+3];

     } else if (!strcmp(elem[i],"jcoupling")) {
       n1=(int)v[1];
       n2=(int)v[2];

       sprintf(name,"jcoupling_%d_%d",n1,n2);
//printf("   readsys doing %s\n",name);
       if (ham_exists(h,name)) {
         fprintf(stderr,"error: spinsys: interaction '%s' does already exists\n",name);
         exit(1);
       }
       if (v[3] == 0.0 && v[4] == 0.0) {
/*
         fprintf(stderr,"warning: spinsys: interaction '%s' has zero coupling strength\n",name);
         exit(1);
*/
       }

       if (ss_issame(ss,n1,n2)) {
         if (ver) printf( "Homonuclear ");
         Ta=II(ss,n1,n2);
         T=T20(ss,n1,n2);
         h->isdiag = kk = 0; 
	 ReT = get_real(T);
	 ReTa = get_real(Ta);     
       } else {
         if (ver) printf( "Heteronuclear ");
         Ta=IzIz(ss,n1,n2);
         T=Dmuld( IzIz(ss,n1,n2), 2.0/sqrt(6.0) );
	 kk = 1;
	 ReT = get_real_diag(T);
	 ReTa = get_real_diag(Ta);
       }
       complx_matrix_free(T);
       complx_matrix_free(Ta);
       if ( (ver) && (!is_liq) ) {
         /* ZT: solid state situation */
         printf( "J-coupling between nucleus %d and %d\n",n1,n2);
         printf( "  isotropic value        : %g Hz\n",v[3]);
         printf( "  anisotropic value      : %g Hz\n",v[4]);
         printf( "  assymmetry             : %g Hz\n",v[5]);
         printf( "  euler angles of tensor :  (%g,%g,%g) degrees\n",v[6],v[7],v[8]); 
       }
       if ( (ver) && (is_liq) ) {
         /* ZT: liquid state situation */
         printf( "J-coupling between nucleus %d and %d\n",n1,n2);
         printf( "  isotropic value        : %g Hz\n",v[3]);
         printf( "  anisotropic part ignored for spectra calculation in liquid\n");
       }

       if (v[3] != 0.0) {
         v[3] *= (2.0*M_PI);
         ham_add_static(h,v[3],ReTa,name,kk);
       } else {
         free((char*)ReTa);
       }
       if ( (v[4] != 0.0) && (!is_liq) ) {
         v[4] *= (2.0*M_PI);
         R=Dtensor2(2.0*v[4],v[5]);
         R = Drot(R,v[6],v[7],v[8]);
         ham_add(h,R,ReT,name,kk);
       } else {
         free((char*)ReT);  
       }
//printf("               ------> done\n");
     } else if (!strncmp(elem[i],"mixing",6)) {
       n_mix++;
       mix[n_mix] = i;
     } else {
       fprintf(stderr,"error: unknown identifier '%s' in spinsys, must be one of\n",elem[i]);
       fprintf(stderr,"       channels, nuclei, shift, dipole, jcoupling, quadrupole,\n");
       fprintf(stderr,"       or mixing_quad_shift, mixing_quad_dipole\n");
       exit(1);
     }
  }

  /* ZT: focus on mixing terms */
  for (kk=1; kk<=n_mix; kk++) {
     i = mix[kk];
     if (!strncmp(elem[i],"mixing_quad_shift",17)) {
       n1=atoi(val[i][1]);
       if (errno) {
          sprintf(buf,"failure in spinsys line %s\n",elem[i]);
	  perror(buf);
       }
       sprintf(name,"mixing_quad_shift_%d",n1);
       if (ham_exists(h,name)) {
         fprintf(stderr,"error: spinsys: mixingterm '%s' does already exists\n",name);
         exit(1);
       }
       if (ver) {
         printf( "Including second-order mixing term between quadrupole and shift for nucleus %d\n",
	     n1);
       }
       if (quadrupoles[n1] == NULL) {
          fprintf(stderr,"error: spinsys - mixing_quad_shift: no quadrupole definitions for nucleus %d\n",n1);
	  exit(1);
       } 
       if (shifts[n1] == NULL) {
          fprintf(stderr,"error: spinsys - mixing_quad_shift: no shift definitions for nucleus %d\n",n1);
	  exit(1);
       } 
       /* operator is */
       T = T20(ss,n1,n1);
       wq=quadrupoles[n1][0]/(4.0*ss_qn(ss,n1)*(2.0*ss_qn(ss,n1)-1));
       v0 = ss_gamma(ss,n1)*s->specfreq/ss_gamma1H()*2.0*M_PI;
       cmv_muld(T,3.0*wq/v0);
       /* spatial components are */
       R=Dtensor2(1.0,quadrupoles[n1][1]);
       R = Drot(R,quadrupoles[n1][2],quadrupoles[n1][3],quadrupoles[n1][4]);
       R2 = Dtensor2(sqrt(2.0/3.0)*shifts[n1][1],shifts[n1][2]);
       R = Drot(R2,shifts[n1][3],shifts[n1][4],shifts[n1][5]);
       /* done. Hopefully correct... */
       ReT = get_real(T);
       ham_add_QC(h,R,R2,ReT,name);
       complx_matrix_free(T);

     } else if (!strncmp(elem[i],"mixing_quad_dipole",18)) {
       qn=atoi(val[i][1]);
       if (errno) {
          sprintf(buf,"failure in spinsys line %s\n",elem[i]);
	  perror(buf);
       }
       dn=atoi(val[i][2]);
       if (errno) {
          sprintf(buf,"failure in spinsys line %s\n",elem[i]);
	  perror(buf);
       }
       sprintf(name,"mixing_quad_dipole_%d_%d",qn,dn);
       if (ham_exists(h,name)) {
         fprintf(stderr,"error: spinsys: mixingterm '%s' does already exist\n",name);
         exit(1);
       }
       if (ver) {
         printf( "Including second-order mixing term between quadrupole nucleus %d and dipole between nuclei %d and %d\n",
	     qn, qn, dn);
       }
       if (quadrupoles[qn] == NULL) {
          fprintf(stderr,"error: spinsys - mixing_quad_dipole: no quadrupole definitions for nucleus %d\n",n1);
	  exit(1);
       } 
       if (qn>dn) {
          n1 = dn;
	  n2 = qn;
       } else {
          n1 = qn;
	  n2 = dn;
       }
       double *dips = dipoles[n1][n2];
       if (dips == NULL) {
          fprintf(stderr,"error: spinsys - mixing_quad_dipole: no dipole definitions between nuclei %d and %d\n",n1, n2);
	  exit(1);
       } 
       wq=quadrupoles[qn][0]/(4.0*ss_qn(ss,qn)*(2.0*ss_qn(ss,qn)-1));
       v0 = ss_gamma(ss,qn)*s->specfreq/ss_gamma1H()*2.0*M_PI;
       /* operators are */
       T = Dmul(T20(ss,qn,qn), Iz(ss,dn));
       cmv_muld(T,sqrt(6.0)*wq/v0);
       ReT = get_real_diag(T);
       if (ss_issame(ss,qn,dn)) {
          /*homocuclear*/
	  Ta = Dmul(Dmul(Ip(ss,qn),Im(ss,dn)),Dadd(Dmuld(Iz(ss,qn),2.0),Ie(ss)));
	  Tb = Dmul(Dmul(Ip(ss,dn),Im(ss,qn)),Dsub(Dmuld(Iz(ss,qn),2.0),Ie(ss)));
	  cmv_muld(Ta,-0.5*wq/v0);
	  cmv_muld(Tb,-0.5*wq/v0);
	  if (!cmv_isreal(Ta)) {
	     fprintf(stderr,"readsys error: mixing term DD/Q matrix Ta is not real\n");
	     exit(1);
	  } 
	  if (!cmv_isreal(Tb)) {
	     fprintf(stderr,"readsys error: mixing term DD/Q matrix Tb is not real\n");
	     exit(1);
	  } 
	  ReTa = get_real(Ta);
	  ReTb = get_real(Tb);
	  complx_matrix_free(Ta);
	  complx_matrix_free(Tb);
       } else {
          /*heteronuclear*/
	  ReTa = NULL;
	  ReTb = NULL;
       }

       /* spatial components are */
       R=Dtensor2(1.0,quadrupoles[qn][1]);
       R = Drot(R,quadrupoles[qn][2],quadrupoles[qn][3],quadrupoles[qn][4]);
       R2 = Dtensor2(2.0*dips[0],0.0);
       R = Drot(R2,dips[1],dips[2],dips[3]);
       /* done. Hopefully correct... */
       ham_add_QD(h,R,R2,ReT,ReTa,ReTb,name);
       
     } else {
       fprintf(stderr,"error: unknown mixing identifier '%s' in spinsys, must be\n",elem[i]);
       fprintf(stderr,"       mixing_quad_shift or mixing_quad_dipole\n");
       exit(1);
     }

  }



//printf("   readsys done all interactions, now ham ini\n");

  /* this will allocate structure for the total Hamiltonian */
  ham_ini_hamilton(h);



//printf("   Done reading spinsystem file.\n");
//printf("   Spin system matrix dimension: %d\n",ss->matdim);

  if (ver) printf("Done reading spinsystem file.\n");
  if (ver) printf("Spin system matrix dimension: %d\n",ss->matdim);
  if (ver) printf("The hamiltonian is%s diagonal\n",(h->isdiag ? "" : " not"));
  if ((various & VARIOUS_NODIAGPROP) && h->isdiag) {
     h->isdiag=0;
     if (ver) printf("Not optimizing for diagonal propagators because of the 'various' parameter setting\n");
  }
  if (s->wr == 0) {
     if (ver) printf("Integration of delays is not done in the static case\n");
     /* ZT: removed the line below and changed code arround _delay */
     /* h->isdiag=0; */
  }
  for (i=1;i<=nelem;i++) {
    /* free(val[i]+1); */
    Tcl_Free((char *)(val[i]+1));
  }
  /* free(elem+1); */
  Tcl_Free((char *)(elem+1));
 
  /* ZT: free interaction maps at the moment. 
   *     REMOVE this if want to use them for relaxation 
   */
  for (i=1; i<=nspins; i++) {
      if (shifts[i]) {
         free((char *)(shifts[i]));
      }
      if (quadrupoles[i]) {
         free((char *)(quadrupoles[i]));
      }
      for (j=1; j<=nspins; j++) {
         if (dipoles[i][j]) {
            free((char *)(dipoles[i][j]));
         }
      }
   }
   free((char *)(dipoles[1]));
   free((char *)(dipoles));
   free((char *)(quadrupoles));
   free((char *)(shifts));
 
}

