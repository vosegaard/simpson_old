/*
    Hamiltonian calculation
    Copyright (C) 1999 Mads Bak, Jimmy T. Rasmussen
                  2002 Thomas Vosegaard

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

#ifndef __HAM_H
#define __HAM_H

#ifndef __COMPLX_H
#include "complx.h"
#endif

#ifndef __SPINSYS_H
#include "spinsys.h"
#endif

#define MAXSPIN 1024
#define MAXQSPIN 20

typedef struct _Hamilton {
  
  /* number of interactions, excluding second order quadrupolar */
  int n;
  char n_names[MAXSPIN][16];
  int n_used[MAXSPIN];
  int n_isdiag[MAXSPIN];
  complx *R[MAXSPIN],*Rrot[MAXSPIN];
  double W[MAXSPIN];
  double *T[MAXSPIN];

  /* number of second order quadrupolar interactions */
  int nq;
  int nq_used[MAXSPIN];
  char nq_names[MAXQSPIN][16];
  complx *RQ[MAXQSPIN],*RQrot[MAXQSPIN],*WQ[MAXQSPIN];
  double *TQ[MAXQSPIN],*TQa[MAXQSPIN],*TQb[MAXQSPIN];

  /* number of second order Q/CSA mixing-term interactions */
  int nqc;
  int nqc_used[MAXSPIN];
  char nqc_names[MAXQSPIN][32];
  complx *RQCq[MAXQSPIN],*RQCqrot[MAXQSPIN],*WQCq[MAXQSPIN];
  complx *RQCc[MAXQSPIN],*RQCcrot[MAXQSPIN],*WQCc[MAXQSPIN];
  double *TQC[MAXQSPIN];

  /* number of second order Q/DD mixing-term interactions */
  int nqd;
  int nqd_used[MAXSPIN];
  char nqd_names[MAXQSPIN][32];
  complx *RQDq[MAXQSPIN],*RQDqrot[MAXQSPIN],*WQDq[MAXQSPIN];
  complx *RQDd[MAXQSPIN],*RQDdrot[MAXQSPIN],*WQDd[MAXQSPIN];
  double *TQD[MAXQSPIN],*TQDa[MAXQSPIN],*TQDb[MAXQSPIN];

  /* number of static contributions, from jiso, csa, offset */
  /* offset is the first static interaction */
  int nstatic;
  int nstatic_used[MAXSPIN];
  char nstatic_names[MAXSPIN][16];
  int nstatic_isdiag[MAXSPIN]; 
  double *Hstatic[MAXSPIN];
  double Hstatic_iso[MAXSPIN]; 
  
   
  /* matrix dimension. whether the hamiltonian is diagonal*/
  int matdim,isdiag;

  /* final hamiltonian */
  mv_double *H;
  
  /* TV: Current euler angles */
  double alpha,beta,gamma;
} Hamilton;


/*
  Allocation strategy:

  matrices and vector added with 
    Hadd* Hsetdiagmatrix
  overgoes to Hamilton to destroy

  The return matrix of Hamilton is internal and must NOT be destroyed
  
  Every vector and matrix is owned by Hamilton and destroyed with Hdestroy
  ham_add* takes over the allocated arguments
  ham_set_offset makes a copy (ham_offset is called during the pulse sequence)
*/

void ham_print(Hamilton* h);
void ham_initialize(Hamilton* h,int matdim);
void ham_destroy(Hamilton* h);
void ham_ini_hamilton(Hamilton *h);

void ham_add_static(Hamilton* h,double iso,double* T,char* name, int isdiag);
void ham_add(Hamilton* h,complx* R,double* T,char* name,int isdiag); 
complx* ham_get(Hamilton* h,char* name); 
void ham_modify(Hamilton* h,complx* R,char* name); 
void ham_add_Q2(Hamilton* h,complx* R,double* T,double* Ta,double* Tb,char* name);
/* ZT: mixing quadrupole/CSA */
void ham_add_QC(Hamilton* h,complx* Rq,complx *Rc,double* T,char* name);
/* ZT: mixing quadrupole/dipole */
void ham_add_QD(Hamilton* h,complx* Rq,complx *Rd,double* T,double* Ta,double* Tb,char* name);


int ham_ischanged(Hamilton* h);
void ham_set_offset(Hamilton* h, double* offset,int used);

void ham_turnoff(Hamilton* h,char* name);
void ham_turnon(Hamilton* h,char* name);
int ham_exists(Hamilton* h,char* name);

void ham_rotate(Hamilton* h,double alpha,double beta,double gamma);
void ham_rotate2(Hamilton* h,double alpha,double beta);
void ham_rotate2_integrate(Hamilton* h, double t1, double dt, double wr,double beta);

mv_double* ham_hamilton(Hamilton* h);
mv_double* ham_hamilton_integrate(Hamilton* h);

void wigner2(complx* d2,double alpha,double beta,double gamma);
void wigner20(complx* d20,double alpha,double beta);
void wig2rot(complx* res, complx* vec, complx *d2);
double wig20rot(complx *vec, complx *d20);

/* DOR is not really finished implementation */
void ham_rotate2_dor(Hamilton* h,double alpha1,double beta1,double alpha2,double beta2);
void ham_rotate2_integrate_dor(Hamilton* h, double t1, double t2,
 double wr1,double beta1,double wr2,double beta2,double dtmax);

 
#endif
