/*
    Hamiltonian calculation routines
    Copyright (C) 1999 Mads Bak, Jimmy T. Rasmussen

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
    
    Setup and calculations for the Hamiltonian for the spinsystem.
    Creates the spin and space tensors, rotates them and calculates
    fourier components.
    
    Called from readsys.c where the spin-system is set up,
    sim.c where the crystallite averaging takes place,
    and pulse.c where the pulse propagation is performed.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "matrix_new.h"
#include "cm_new.h"
#include "ham.h"
#include "defs.h"

/*
  Allocation strategy:

  matrices and vector added with 
    ham_add* _ham_setdiagmatrix
  overgoes to ham_hamilton to destroy

  *rot matrices is allocated internal in Hamilton

  The return matrix of ham_hamilton is internal and must NOT be destroyed
  
  Every vector and matrix is owned by Hamilton and destroyed with ham_destroy
*/

void ham_initialize(Hamilton* h,int matdim)
{
  h->n=0;
  h->nstatic=0;
  h->nq=0;
  h->nqc=0;
  h->nqd=0;
  h->isdiag=0;
  h->matdim=matdim;

  /* the first static interaction is DC offset 
     it is allocated, but turned off by default */
  double *ptr;
  ptr = (double*)malloc(matdim*sizeof(double));
  ham_add_static(h,0.0,ptr,"offset",1);
  ham_set_offset(h,NULL,0); 
} 
/* this I want to call after readsys is done */
void ham_ini_hamilton(Hamilton *h)
{
 int matdim = h->matdim;
 
 if (h->isdiag) {
     h->H = double_vector_alloc(matdim);
  } else {
     h->H = double_matrix_alloc(matdim,matdim);
  }
}

void ham_destroy(Hamilton* h)
{
  int i;

  for (i=1;i<=h->nstatic;i++) {
    free((char*)h->Hstatic[i]);
  }
  for (i=1;i<=h->n;i++) {
    free((char*)h->R[i]);
    free((char*)h->Rrot[i]);
    free((char*)h->T[i]);
  }
  for (i=1;i<=h->nq;i++) {
    free((char*)h->RQ[i]);
    free((char*)h->RQrot[i]);
    free((char*)h->WQ[i]);
    free((char*)h->TQ[i]);
    free((char*)h->TQa[i]);
    free((char*)h->TQb[i]);
  }
  for (i=1;i<=h->nqc;i++) {
    free((char*)h->RQCq[i]);
    free((char*)h->RQCqrot[i]);
    free((char*)h->WQCq[i]);
    free((char*)h->RQCc[i]);
    free((char*)h->RQCcrot[i]);
    free((char*)h->WQCc[i]);
    free((char*)h->TQC[i]);
  }
  for (i=1;i<=h->nqd;i++) {
    free((char*)h->RQDq[i]);
    free((char*)h->RQDqrot[i]);
    free((char*)h->WQDq[i]);
    free((char*)h->RQDd[i]);
    free((char*)h->RQDdrot[i]);
    free((char*)h->WQDd[i]);
    free((char*)h->TQD[i]);
    free((char*)h->TQDa[i]);
    free((char*)h->TQDb[i]);
  }
  /* this works both for matrix and vector */
  double_matrix_free(h->H);
  
}

#define check_dimensions(T) if ((ROWS((T))) != h->matdim) {\
    fprintf(stderr,"error: ham.c: conflict in matrix dimensions");\
    exit(1);\
  }

void ham_add_static(Hamilton* h,double iso,double* T,char* name, int isdiag) 
{
/*  check_dimensions(T); */
  int len;
  
  if (isdiag) 
     len = h->matdim;
  else
     len = h->matdim*h->matdim;

  h->nstatic++;
  h->nstatic_used[h->nstatic]=1;
  /* m_mulmc(T,T,iso); */
  dscal_(&len,&iso,T,&INTONE);
  h->Hstatic[h->nstatic]=T;
  h->Hstatic_iso[h->nstatic]=iso;
  strcpy(h->nstatic_names[h->nstatic],name);
  h->nstatic_isdiag[h->nstatic]=isdiag;
}

void ham_add(Hamilton* h,complx* R,double *T,char* name, int isdiag) 
{
  /* check_dimensions(T); */

  h->n++;
  h->n_used[h->n]=1;
  h->R[h->n]=R;
  h->T[h->n]=T;
  h->Rrot[h->n]=(complx*)malloc(5*sizeof(complx));
  memset(h->Rrot[h->n],0,5*sizeof(complx));
  strcpy(h->n_names[h->n],name);
  h->n_isdiag[h->n]=isdiag;
/*int NN=(h->matdim*h->matdim-h->matdim)*(1-isdiag)+h->matdim;
int kk; for (kk=0; kk<NN; kk++) printf("H++ %s: (%i) = %g\n",name,kk,T[kk]);
for (kk=0; kk<5; kk++) printf("        R[%i] = %g\n",kk,R[kk]); */
}

/* Two new functions for dynamics */

complx* ham_get(Hamilton* h,char* name) 
{
  int i;

  for (i=1; i<=h->n; i++) {
    if (strcmp(h->n_names[h->n],name)) continue;
    return h->R[i];
  }
  return NULL;
}

void ham_modify(Hamilton* h,complx* R,char* name) 
{
  int i;

  for (i=1; i<=h->n; i++) {
    if (strcmp(h->n_names[h->n],name)) continue;
    h->n_used[i]=1;
    h->R[i]=R;
    return;
  }
}

/* End of new stuff */

void ham_add_Q2(Hamilton* h,complx* R,double* T,double* Ta,double* Tb,char* name) 
{
  /* check_dimensions(T); */
  
  /* this is always diagonal interaction */
  h->nq++;
  h->RQ[h->nq]=R;
  h->TQ[h->nq]=T;
  h->TQa[h->nq]=Ta;
  h->TQb[h->nq]=Tb;
  h->RQrot[h->nq]=(complx*)malloc(5*sizeof(complx));
  memset(h->RQrot[h->nq],0,5*sizeof(complx));
  h->WQ[h->nq]=(complx*)malloc(5*sizeof(complx));
  memset(h->WQ[h->nq],0,5*sizeof(complx));
  strcpy(h->nq_names[h->nq],name);
/*int kk; for (kk=0; kk<h->matdim; kk++) printf("H++ %s: (%i) = %g | %g | %g\n",name,kk,T[kk],Ta[kk],Tb[kk]);
for (kk=0; kk<5; kk++) printf("        R[%i] = %g, %g\n",kk,R[kk].re,R[kk].im); */
}

/* ZT: mixing quadrupole/CSA */
void ham_add_QC(Hamilton* h,complx* Rq,complx *Rc,double* T,char* name) 
{
  /* check_dimensions(T);  */

  /* this is always diagonal interaction */
  h->nqc++;
  h->RQCq[h->nqc]=Rq;
  h->RQCc[h->nqc]=Rc;
  h->TQC[h->nqc]=T;
  h->RQCqrot[h->nqc]=(complx*)malloc(5*sizeof(complx));
  memset(h->RQCqrot[h->nqc],0,5*sizeof(complx));
  h->WQCq[h->nqc]=(complx*)malloc(5*sizeof(complx));
  memset(h->WQCq[h->nqc],0,5*sizeof(complx));
  h->RQCcrot[h->nqc]=(complx*)malloc(5*sizeof(complx));
  memset(h->RQCcrot[h->nqc],0,5*sizeof(complx));
  h->WQCc[h->nqc]=(complx*)malloc(5*sizeof(complx));
  memset(h->WQCc[h->nqc],0,5*sizeof(complx));
  strcpy(h->nqc_names[h->nqc],name);
}

/* ZT: mixing quadrupole/dipole */
void ham_add_QD(Hamilton* h,complx* Rq,complx *Rd,double* T,double* Ta,double* Tb,char* name) 
{
  /* check_dimensions(T); */

  /* TQD is diagonal, TQDa and TQDb are not diagonal */
  h->nqd++;
  h->RQDq[h->nqd]=Rq;
  h->RQDd[h->nqd]=Rd;
  h->TQD[h->nqd]=T;
  h->TQDa[h->nqd]=Ta;
  h->TQDb[h->nqd]=Tb;
  h->RQDqrot[h->nqd]=(complx*)malloc(5*sizeof(complx));
  memset(h->RQDqrot[h->nqd],0,5*sizeof(complx));
  h->WQDq[h->nqd]=(complx*)malloc(5*sizeof(complx));
  memset(h->WQDq[h->nqd],0,5*sizeof(complx));
  h->RQDdrot[h->nqd]=(complx*)malloc(5*sizeof(complx));
  memset(h->RQDdrot[h->nqd],0,5*sizeof(complx));
  h->WQDd[h->nqd]=(complx*)malloc(5*sizeof(complx));
  memset(h->WQDd[h->nqd],0,5*sizeof(complx));
  strcpy(h->nqd_names[h->nqd],name);
}

void ham_set_offset(Hamilton* h,double* offset,int used) 
{
  h->nstatic_used[1]=used;
  if (used) {
    /* check_dimensions(offset); */ 
    memcpy(h->Hstatic[1],offset,h->matdim*sizeof(double));
  }
}

int ham_exists(Hamilton* h,char* name)
{
  int i;

  for (i=1;i<=h->n;i++) {
    if (!strcmp(name,h->n_names[i])) return 1;
  }
  for (i=1;i<=h->nstatic;i++) {
    if (!strcmp(name,h->nstatic_names[i])) return 1;
  }
  for (i=1;i<=h->nq;i++) {
    if (!strcmp(name,h->nq_names[i])) return 1;
  }
  for (i=1;i<=h->nqc;i++) {
    if (!strcmp(name,h->nqc_names[i])) return 1;
  }
  for (i=1;i<=h->nqd;i++) {
    if (!strcmp(name,h->nqd_names[i])) return 1;
  }
  return 0;
}

int ham_ischanged(Hamilton* h)
{
  int i;

  /* return true if offset is set*/
  /* ZT: static field inhomogeneity key: 
   *       used = 1 offset set by offset, then return true
   *       used = 2 offset set to nominal value, return false
   */
  if (h->nstatic_used[1] == 1) return 1;

  /* return true if any interactions are disabled */
  for (i=1;i<=h->n;i++) {
    if (!h->n_used[i]) return 1;
  }
  for (i=2;i<=h->nstatic;i++) {
    if (!h->nstatic_used[i]) return 1;
  }
  for (i=1;i<=h->nq;i++) {
    if (!h->nq_used[i]) return 1;
  }
  for (i=1;i<=h->nqc;i++) {
    if (!h->nqc_used[i]) return 1;
  }
  for (i=1;i<=h->nqd;i++) {
    if (!h->nqd_used[i]) return 1;
  }
  return 0;
}

void ham_turnoff(Hamilton* h,char* name)
{
  int i,was;

  /* the same interaction name can be present in all three kinds
     of interactions, i.e. csa gives a contribution to n and nstatic*/
  if (!strcmp(name,"all")) {
    for (i=1;i<=h->n;i++)
      h->n_used[i]=0;

    /* DC offset (in nstatic[1]) is not changed because
       it isnt an interaction specified in the spinsys section */
    for (i=2;i<=h->nstatic;i++)
      h->nstatic_used[i]=0;
    for (i=1;i<=h->nq;i++)
      h->nq_used[i]=0;
    for (i=1;i<=h->nqc;i++)
      h->nqc_used[i]=0;
    for (i=1;i<=h->nqd;i++)
      h->nqd_used[i]=0;
    return;
  }
  was=0;
  for (i=1;i<=h->n;i++) {
    if (!strcmp(name,h->n_names[i])) {
      h->n_used[i]=0;
      was=1;
      break;
    }
  }
  for (i=1;i<=h->nstatic;i++) {
    if (!strcmp(name,h->nstatic_names[i])) {
      h->nstatic_used[i]=0;
      was=1;
      break;
    }
  }
  for (i=1;i<=h->nq;i++) {
    if (!strcmp(name,h->nq_names[i])) {
      h->nq_used[i]=0;
      was=1;
      break;
    }
  }
  for (i=1;i<=h->nqc;i++) {
    if (!strcmp(name,h->nqc_names[i])) {
      h->nqc_used[i]=0;
      was=1;
      break;
    }
  }
  for (i=1;i<=h->nqd;i++) {
    if (!strcmp(name,h->nqd_names[i])) {
      h->nqd_used[i]=0;
      was=1;
      break;
    }
  }
  if (!was) {
    fprintf(stderr,"error: turnoff: unknown interaction name '%s'\n",name);
    exit(1);
  }
}

void ham_turnon(Hamilton* h,char* name)
{
  int i,was;

  /* the same interaction name can be present in all three kinds
     of interactions, i.e. csa gives a contribution to n and nstatic*/
  if (!strcmp(name,"all")) {
    for (i=1;i<=h->n;i++)
      h->n_used[i]=1;

    /* DC offset (in nstatic[1]) is not changed because
       it isnt an interaction specified in the spinsys section */
    for (i=2;i<=h->nstatic;i++)
      h->nstatic_used[i]=1;
    for (i=1;i<=h->nq;i++)
      h->nq_used[i]=1;
    for (i=1;i<=h->nqc;i++)
      h->nqc_used[i]=1;
    for (i=1;i<=h->nqd;i++)
      h->nqd_used[i]=1;

    return;
  }
  was=0;
  for (i=1;i<=h->n;i++) {
    if (!strcmp(name,h->n_names[i])) {
      h->n_used[i]=1;
      was=1;
      break;
    }
  }
  for (i=1;i<=h->nstatic;i++) {
    if (!strcmp(name,h->nstatic_names[i])) {
      h->nstatic_used[i]=1;
      was=1;
      break;
    }
  }
  for (i=1;i<=h->nq;i++) {
    if (!strcmp(name,h->nq_names[i])) {
      h->nq_used[i]=1;
      was=1;
      break;
    }
  }
  for (i=1;i<=h->nqc;i++) {
    if (!strcmp(name,h->nqc_names[i])) {
      h->nqc_used[i]=1;
      was=1;
      break;
    }
  }
  for (i=1;i<=h->nqd;i++) {
    if (!strcmp(name,h->nqd_names[i])) {
      h->nqd_used[i]=1;
      was=1;
      break;
    }
  }
  if (!was) {
    fprintf(stderr,"error: turnoff: unknown interaction name '%s'\n",name);
    exit(1);
  }
}

void ham_rotate(Hamilton* h,double alpha,double beta,double gamma)
{
  int i;
  static mv_complx *d2=NULL;
  
  d2=cmv_static(d2,5,5);
  
  wigner2(d2->data,alpha,beta,gamma);
  /* This procedure is called outside the pulse sequence and
     must therefore be performed even if the interactions are turned off. */
  for (i=1;i<=h->n;i++) {
     wig2rot(h->Rrot[i],h->R[i],d2->data);
  }
  for (i=1;i<=h->nq;i++) {
     wig2rot(h->RQrot[i],h->RQ[i],d2->data);
  }
  for (i=1;i<=h->nqc;i++) {
     wig2rot(h->RQCqrot[i],h->RQCq[i],d2->data);
     wig2rot(h->RQCcrot[i],h->RQCc[i],d2->data);
  }
  for (i=1;i<=h->nqd;i++) {
     wig2rot(h->RQDqrot[i],h->RQDq[i],d2->data);
     wig2rot(h->RQDdrot[i],h->RQDd[i],d2->data);
  }
  
}

void ham_rotate2(Hamilton* h,double alpha,double beta) 
{
  int i;
  static mv_complx *d20=NULL,*d2=NULL;

  d20= cmv_static(d20,5,1);
  d2 = cmv_static(d2,5,5);
  wigner20(d20->data,alpha,beta);
  for (i=1;i<=h->n;i++) {
    if (!h->n_used[i]) continue;
  	h->W[i]=wig20rot(h->Rrot[i],d20->data);
  }
  if (h->nq) {
    wigner2(d2->data,alpha,beta,0.0);
    for (i=1;i<=h->nq;i++) {
      if (!h->nq_used[i]) continue;
      wig2rot(h->WQ[i],h->RQrot[i],d2->data);
    }
  }
  if (h->nqc) {
    wigner2(d2->data,alpha,beta,0.0);
    for (i=1;i<=h->nqc;i++) {
      if (!h->nqc_used[i]) continue;
      wig2rot(h->WQCq[i],h->RQCqrot[i],d2->data);
      wig2rot(h->WQCc[i],h->RQCcrot[i],d2->data);
    }
  }
  if (h->nqd) {
    wigner2(d2->data,alpha,beta,0.0);
    for (i=1;i<=h->nqd;i++) {
      if (!h->nqd_used[i]) continue;
      wig2rot(h->WQDq[i],h->RQDqrot[i],d2->data);
      wig2rot(h->WQDd[i],h->RQDdrot[i],d2->data);
    }
  }
}

void ham_rotate2_dor(Hamilton* h,double alpha1,double beta1,
                     double alpha2, double beta2)
{
  int i;
  static mv_complx *d20=NULL,*d2_1=NULL,*d2_2=NULL,*d2=NULL;

  d20 = cmv_static(d20,5,1);
  d2_1 = cmv_static(d2_1,5,5);
  d2_2 = cmv_static(d2_2,5,5);
  d2 = cmv_static(d2,5,5);
  
  wigner2(d2_1->data,alpha1,beta1,0.0);
  wigner2(d2_2->data,alpha2,beta2,0.0);
  cm_mul(d2,d2_1,d2_2);
  memcpy(d20->data,d2->data + 10,5*sizeof(complx));
  /*for (i=1;i<=5;i++) {
    d20[i]=d2[i][3];
  }*/

  for (i=1;i<=h->n;i++) {
    if (!h->n_used[i]) continue;
  	h->W[i]=wig20rot(h->Rrot[i],d20->data);
  }
  if (h->nq) {
    for (i=1;i<=h->nq;i++) {
      if (!h->nq_used[i]) continue;
      wig2rot(h->WQ[i],h->RQrot[i],d2->data);
    }
  }
  if (h->nqc) {
    for (i=1;i<=h->nqc;i++) {
      if (!h->nqc_used[i]) continue;
      wig2rot(h->WQCq[i],h->RQCqrot[i],d2->data);
      wig2rot(h->WQCc[i],h->RQCcrot[i],d2->data);
    }
  }
  if (h->nqd) {
    for (i=1;i<=h->nqd;i++) {
      if (!h->nqd_used[i]) continue;
      wig2rot(h->WQDq[i],h->RQDqrot[i],d2->data);
      wig2rot(h->WQDd[i],h->RQDdrot[i],d2->data);
    }
  }
}

/*
complx mvv(complx* a,complx *b)
{
   int i;
   complx c,*pa,*pb;

   c.re=0;
   c.im=0;
   pa = &a[1];
   pb = &b[1];
   
   for (i=1;i<=5;i++,pa++,pb++) {
     c.re += pa->re*pb->re-pa->im*pb->im;
     c.im += pa->im*pb->re+pa->re*pb->im;
   }
   return c;
}
*/

/*
   int_t1^t2 exp(-i m wr t) = -i/(m wr) (exp(-i m wr t2) - exp(-i m wr t1))
*/
void integ(complx* R,double t1,double t2,double wr)
{
   double wrm;
/*
   wrm= wr * (-2.0); R[1] = Cmul(Complx(0,-1.0/wrm),Csub(Cexpi(-t2*wrm),Cexpi(-t1*wrm))); 
   wrm= wr * (-1.0); R[2] = Cmul(Complx(0,-1.0/wrm),Csub(Cexpi(-t2*wrm),Cexpi(-t1*wrm)));
                     R[3] = Complx(t2-t1,0);
   wrm= wr * ( 1.0); R[4] = Cmul(Complx(0,-1.0/wrm),Csub(Cexpi(-t2*wrm),Cexpi(-t1*wrm))); 
   wrm= wr * ( 2.0); R[5] = Cmul(Complx(0,-1.0/wrm),Csub(Cexpi(-t2*wrm),Cexpi(-t1*wrm))); 
*/
   wrm= wr * (-2.0); R[0] = Cmul(Complx(0,1.0/wrm),Csub(Cexpi(-t2*wrm),Cexpi(-t1*wrm)));
   wrm= wr * (-1.0); R[1] = Cmul(Complx(0,1.0/wrm),Csub(Cexpi(-t2*wrm),Cexpi(-t1*wrm)));
                     R[2] = Complx(t2-t1,0);
 /*  wrm= wr * ( 1.0); R[3] = Cmul(Complx(0,1.0/wrm),Csub(Cexpi(-t2*wrm),Cexpi(-t1*wrm)));
   wrm= wr * ( 2.0); R[4] = Cmul(Complx(0,1.0/wrm),Csub(Cexpi(-t2*wrm),Cexpi(-t1*wrm)));
  */
  
  R[3] = Conj(R[1]);
  R[4] = Conj(R[0]);
}

void ham_rotate2_integrate(Hamilton* h, double t1, double dt,double wr,double beta) 
{
  int i,j,k,l,m;
  static mv_complx *d20=NULL, *d2=NULL;
  double t2;
  complx R[5],K[5],I[9],c1,c2,c4,c5,cc1,cc2;
  double wrm;
  int kk;

  d20 = cmv_static(d20,5,1);
  d2 = cmv_static(d2,5,5);

  t2=t1+dt;

  wigner20(d20->data,0.0,beta);
  integ(R,t1,t2,wr);    

  for (i=1;i<=h->n;i++) {

    if (!h->n_used[i]) continue;

    for (j=0;j<5;j++) {
      K[j]=Cmul(h->Rrot[i][j],d20->data[j]); 
    }
   /* h->W[i]=mvv(R,K); */
    h->W[i] = wig20rot(R,K);
    h->W[i] = h->W[i]/dt; /* Divide by dt to later multiply
                                      by the same in m_realexpdiag */
  }

  if (h->nq || h->nqc || h->nqd) {
     wigner2(d2->data,0,beta,0.0);
     for (l=0;l<4;l++) {
        kk= l-4;
        wrm  = wr * kk;
        I[l]=Cmul(Complx(0,1.0/wrm/dt),Csub(Cexpi(-t2*wrm),Cexpi(-t1*wrm))); 
     }
     I[4] = Complx(1.0,0.0);
     I[5] = Conj(I[3]);
     I[6] = Conj(I[2]);
     I[7] = Conj(I[1]);
     I[8] = Conj(I[0]);
  }

  if (h->nq) {
    for (i=1;i<=h->nq;i++) {

      if (!h->nq_used[i]) continue;
      memset(h->WQ[i],0,5*sizeof(complx));
      cc1=Complx(0,0);
      cc2=Complx(0,0);
      for (l=0;l<5;l++) {
        c2= Cmul(h->RQrot[i][l],d2->data[l+1*5]);
        c1= Cmul(h->RQrot[i][l],d2->data[l+0*5]);
        for (m=0;m<5;m++) {
          c4=Cmul(h->RQrot[i][m],d2->data[m+3*5]);
          c5=Cmul(h->RQrot[i][m],d2->data[m+4*5]);
          cc2=Cadd(cc2,Cmul(Cmul(c2,c4),I[l+m]));
          cc1=Cadd(cc1,Cmul(Cmul(c1,c5),I[l+m]));
        }
      }
      if (fabs(cc1.im)+fabs(cc2.im)>1e-8) {
        fprintf(stderr,"ham_rotate2_integrate error: cc1 and cc2 are not real\n");
        exit(1);
      }
      /* ZT: this would make sense to me */
         /* h->WQ[i][2].re = cc2.re/dt;        
            h->WQ[i][1].re = cc1.re/dt;     */     
         /* but original simpson had this   */
         h->WQ[i][2].re = cc2.re;        
         h->WQ[i][1].re = cc1.re; 

        for (j=0;j<5;j++) {
          K[j]=Cmul(h->RQrot[i][j],d20->data[j]); 
        }    
        /* ZT: this would make sense to me    */
        /* h->WQ[i][0].re = wig20rot(R,K)/dt; */
        /* but original simpson had this      */
        /* h->WQ[i][0].re = -wig20rot(R,K)/dt; */
      h->WQ[i][0].re = wig20rot(R, K)/dt;
    }
  }

  if (h->nqc) {

    for (i=1;i<=h->nqc;i++) {

      if (!h->nqc_used[i]) continue;
      
         memset(h->WQCc[i],0,5*sizeof(complx));
         memset(h->WQCq[i],0,5*sizeof(complx));
	 cc1 = Complx(0.0,0.0);
         for (l=0;l<5;l++) {
           K[0] = Cmul(h->RQCqrot[i][l],d2->data[l+1*5]); /* d2(l,-1) */
           K[1] = Cmul(h->RQCqrot[i][l],d2->data[l+3*5]); /* d2(l,+1) */
           for (m=0;m<5;m++) {
              K[2] = Cmul(h->RQCcrot[i][m],d2->data[m+1*5]); /* d2(m,-1) */
              K[3] = Cmul(h->RQCcrot[i][m],d2->data[m+3*5]); /* d2(m,+1) */
              cc2 = Cadd(Cmul(K[1],K[2]),Cmul(K[0],K[3]));
              cc1 = Cadd(cc1,Cmul(cc2,I[l+m]));
           }
         }
	 if (fabs(cc1.im)>1e-8) {
	    fprintf(stderr,"ham_rotate2_integrate error: cc1 is not real (QC mixing)\n");
	    exit(1);
	 }
         h->WQCq[i][0].re = cc1.re/dt;
    }
  }


  if (h->nqd) {
    for (i=1;i<=h->nqd;i++) {

      if (!h->nqc_used[i]) continue;
         complx KK[4], cc3, cc4;
	 
         memset(h->WQDd[i],0,5*sizeof(complx));
         memset(h->WQDq[i],0,5*sizeof(complx));
         for (l=0;l<5;l++) {
           K[0] = Cmul(h->RQDqrot[i][l],d2->data[l+1*5]); /* d2(l,-1) */
           K[1] = Cmul(h->RQDqrot[i][l],d2->data[l+3*5]); /* d2(l,+1) */
           K[2] = Cmul(h->RQDqrot[i][l],d2->data[l+0*5]); /* d2(l,-2) */
           K[3] = Cmul(h->RQDqrot[i][l],d2->data[l+4*5]); /* d2(l,+2) */
           for (m=0;m<5;m++) {
              KK[0] = Cmul(h->RQDdrot[i][m],d2->data[m+1*5]); /* d2(m,-1) */
              KK[1] = Cmul(h->RQDdrot[i][m],d2->data[m+3*5]); /* d2(m,+1) */
              KK[2] = Cmul(h->RQDdrot[i][m],d2->data[m+0*5]); /* d2(m,-2) */
              KK[3] = Cmul(h->RQDdrot[i][m],d2->data[m+4*5]); /* d2(m,+2) */
              cc1 = Cmul(K[0],KK[1]);
	      cc2 = Cmul(K[1],KK[0]);
              cc3 = Cmul(K[2],KK[3]);
	      cc4 = Cmul(K[3],KK[2]);
	      h->WQDq[i][0] = Cadd(h->WQDq[i][0],Cmul(Cadd(cc1,cc2),I[l+m]));
	      h->WQDq[i][1] = Cadd(h->WQDq[i][1],Cmul(Cadd(cc1,cc3),I[l+m]));
	      h->WQDq[i][2] = Cadd(h->WQDq[i][2],Cmul(Cadd(cc2,cc4),I[l+m]));
           }
         }
	 if (fabs(h->WQDq[i][0].im)+fabs(h->WQDq[i][1].im)+fabs(h->WQDq[i][2].im)>1e-8) {
	    fprintf(stderr,"ham_rotate2_integrate error: coefs are not real (QD mixing)\n");
	    exit(1);
	 }
         h->WQDq[i][0].re /= dt;
         h->WQDq[i][1].re /= dt;
         h->WQDq[i][2].re /= dt;
    }
  }
}


mv_double * ham_hamilton(Hamilton* h) 
{
  int i, N, NN;
  static mv_double * hdg=NULL;
  const double done=1.0;
  double dw;
  complx cdw;
  double *dm, *dv, *stop;
  
  N = h->matdim;
  NN = N*N;
  if (h->isdiag) {
     hdg = h->H;
  } else {
     hdg = dmv_static(hdg,N,1);
     dmv_zero(hdg);
  }
  dmv_zero(h->H);
  
  for (i=1;i<=h->nstatic;i++) {
    if (!h->nstatic_used[i]) continue;
    if (h->nstatic_isdiag[i]) {
       daxpy_(&N,&done,h->Hstatic[i],&INTONE,hdg->data,&INTONE);
    } else {
       daxpy_(&NN,&done,h->Hstatic[i],&INTONE,h->H->data,&INTONE);
    }
  }
  for (i=1;i<=h->n;i++) {
    if (!h->n_used[i]) continue;
    if (h->n_isdiag[i]) {
       daxpy_(&N,&(h->W[i]),h->T[i],&INTONE,hdg->data,&INTONE);
    } else {
       daxpy_(&NN,&(h->W[i]),h->T[i],&INTONE,h->H->data,&INTONE);
    }
  }
  for (i=1;i<=h->nq;i++) {
    if (!h->nq_used[i]) continue;
    if (h->WQ[i][2].im > 1e-8) {
       fprintf(stderr,"ham_hamilton error: WQ[%i][2] is not real\n",i);
       exit(1);
    }
    dw = h->WQ[i][2].re;
    daxpy_(&N,&dw,h->TQ[i],&INTONE,hdg->data,&INTONE);
    cdw = Cmul(h->WQ[i][0],h->WQ[i][4]);
    if (cdw.im > 1e-8) {
       fprintf(stderr,"ham_hamilton error: TQa coef. is not real\n",i);
       exit(1);
    }
    daxpy_(&N,&(cdw.re),h->TQa[i],&INTONE,hdg->data,&INTONE);
    cdw = CRmul(Cmul(h->WQ[i][1],h->WQ[i][3]),2.0);
    if (cdw.im > 1e-8) {
       fprintf(stderr,"ham_hamilton error: TQb coef. is not real\n",i);
       exit(1);
    }
    daxpy_(&N,&(cdw.re),h->TQb[i],&INTONE,hdg->data,&INTONE);
  }
  for (i=1;i<=h->nqc;i++) {
    if (!h->nqc_used[i]) continue;
    cdw = Cmul(h->WQCq[i][1],h->WQCc[i][3]);
    cdw = Cadd(cdw,Cmul(h->WQCq[i][3],h->WQCc[i][1]));
    if (cdw.im > 1e-8) {
       fprintf(stderr,"ham_hamilton error: TQC coef. is not real\n",i);
       exit(1);
    }
    daxpy_(&N,&(cdw.re),h->TQC[i],&INTONE,hdg->data,&INTONE);
    printf("ham_hamilton: Q/CSA\n");
  }
  for (i=1;i<=h->nqd;i++) {
    if (!h->nqd_used[i]) continue;
    cdw = Cmul(h->WQDq[i][1],h->WQDd[i][3]);
    cdw = Cadd(cdw,Cmul(h->WQDq[i][3],h->WQDd[i][1]));
    if (cdw.im > 1e-8) {
       fprintf(stderr,"ham_hamilton error: TQD coef. is not real\n",i);
       exit(1);
    }
    daxpy_(&N,&(cdw.re),h->TQD[i],&INTONE,hdg->data,&INTONE);
    printf("ham_hamilton: Q/DD common-nuclear\n");
    if (h->TQDa[i]) {
       cdw = Cmul(h->WQDq[i][1],h->WQDd[i][3]);
       cdw = Cadd(cdw,Cmul(h->WQDq[i][0],h->WQDd[i][4]));
       if (cdw.im > 1e-8) {
          fprintf(stderr,"ham_hamilton error: TQDa coef. is not real\n",i);
          exit(1);
       }
       daxpy_(&NN,&(cdw.re),h->TQDa[i],&INTONE,h->H->data,&INTONE);
       cdw = Cmul(h->WQDq[i][3],h->WQDd[i][1]);
       cdw = Cadd(cdw,Cmul(h->WQDq[i][4],h->WQDd[i][0])); 
       if (cdw.im > 1e-8) {
          fprintf(stderr,"ham_hamilton error: TQDb coef. is not real\n",i);
          exit(1);
       }
       daxpy_(&NN,&(cdw.re),h->TQDb[i],&INTONE,h->H->data,&INTONE);
       printf("ham_hamilton: Q/DD homonuclear\n");
    }
  }
  
  if (!h->isdiag) {
     /* sum diagonal and full contributions */
     dm = h->H->data;
     dv = hdg->data;
     stop = dv + N;
     do {
        *dm += *dv;
        dv++;
        dm += N+1;
     } while ( dv != stop);
  }
  return h->H;
}


mv_double * ham_hamilton_integrate(Hamilton* h) 
{
  int i, N, NN;
  static mv_double * hdg;
  double *dm, *dv, *stop, dw;
  const double done=1.0;

  N = h->matdim;
  NN = N*N;
  if (h->isdiag) {
     hdg = h->H;
  } else {
     hdg = dmv_static(hdg,N,1);
     dmv_zero(hdg);
  }
  dmv_zero(h->H);

  for (i=1;i<=h->nstatic;i++) {
    if (!h->nstatic_used[i]) continue;
    if (h->nstatic_isdiag[i]) {
       daxpy_(&N,&done,h->Hstatic[i],&INTONE,hdg->data,&INTONE);
    } else {
       daxpy_(&NN,&done,h->Hstatic[i],&INTONE,h->H->data,&INTONE);
    }
  }
  for (i=1;i<=h->n;i++) {
    if (!h->n_used[i]) continue;
    if (h->n_isdiag[i]) {
       daxpy_(&N,&(h->W[i]),h->T[i],&INTONE,hdg->data,&INTONE);
    } else {
       daxpy_(&NN,&(h->W[i]),h->T[i],&INTONE,h->H->data,&INTONE);
    }
  }
  for (i=1;i<=h->nq;i++) {
    if (!h->nq_used[i]) continue;

    dw = h->WQ[i][0].re;
    daxpy_(&N,&dw,h->TQ[i],&INTONE,hdg->data,&INTONE);
    dw = h->WQ[i][1].re;
    daxpy_(&N,&dw,h->TQa[i],&INTONE,hdg->data,&INTONE);
    dw = h->WQ[i][2].re*2.0;
    daxpy_(&N,&dw,h->TQb[i],&INTONE,hdg->data,&INTONE);
  }
  for (i=1;i<=h->nqc;i++) {
    if (!h->nqc_used[i]) continue;
    dw = h->WQCq[i][0].re;
    daxpy_(&N,&dw,h->TQC[i],&INTONE,hdg->data,&INTONE);
  }
  for (i=1;i<=h->nqd;i++) {
    dw = h->WQDq[i][0].re;
    daxpy_(&N,&dw,h->TQD[i],&INTONE,hdg->data,&INTONE);
    if (h->TQDa[i]) {
       dw = h->WQDq[i][1].re;
       daxpy_(&NN,&dw,h->TQDa[i],&INTONE,h->H->data,&INTONE);
       dw = h->WQDq[i][2].re;
       daxpy_(&NN,&dw,h->TQDb[i],&INTONE,h->H->data,&INTONE);
    }
  }

  if (!(h->isdiag)) {
     /* sum diagonal and full contributions */
     dm = h->H->data;
     dv = hdg->data;
     stop = dv + N;
     do {
        *dm += *dv;
        dv++;
        dm += N+1;
     } while ( dv != stop);
  }

  return h->H;
}






