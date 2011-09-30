/*
    FID calculation procedures
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
    
    
    Here are located to ways to make the free induction decay evolution.
    - gcompute is an gamma averaged COMPUTE algorithm.
    - gammarep reuses propagators calculated for each gamma angle.
    
    They are called from sim.c
*/

#include <stdlib.h>
#include <stdio.h>
#include "complx.h"
#include "cm_new.h"
#include "fidcalc.h"

/* RHO symmetry is if 
     fstart = 0.5 (fdetect + ajoint(fdetect))
*/
int is_rhosymmetry(mv_complx* fstart,mv_complx* fdetect)
{
  int N,i,j;
  double sum;
  mv_complx *tmp;

  tmp = complx_matrix_alloc(fdetect->row,fdetect->col);
  cm_adjoint(tmp,fdetect);
  cmv_addto(tmp,fdetect);
  cmv_muld(tmp,-0.5);
  cmv_addto(tmp,fstart);
  sum = cmv_sumnorm1(tmp);
  complx_matrix_free(tmp);   
  return (sum < 0.0001 ? 1 : 0);
}

#define NAMAX 1024

void gcompute(mv_complx** ustore,
	      int na,mv_complx* rho,mv_complx* obs,double tt,
	      complx* fid,int rhosymmetry,int realspec)
{
  int nt,N,m,n,i,j,k,l,p,ct,yy;
  int ja,jb,nareal=na;
  double tr,fac;
  complx ff,kk,Itr,*e;

  static mv_complx *tmp=NULL, *tmp2=NULL;
  static mv_complx *qt[NAMAX],*favg[NAMAX],*rhot[NAMAX];
  static int first=1;
  static mv_complx *uc=NULL,*evec=NULL,*evec_inv=NULL,*omega=NULL;
  static mv_int *fok=NULL;
  static mv_complx *evalv=NULL;
  complx dumz;
  
  tr=(double)na*tt;
  nt=LEN(fid);
  N = rho->row;
  
  m_zerov(fid);

  if (realspec) {
    if ((na%2)==0) {
      nareal=na/2;
    } else {
      nareal=((na+1)/2);
    }
  }
  
  uc = cmv_static(uc,N,N);
  tmp = cmv_static(tmp,N,N);
  tmp2 = cmv_static(tmp2,N,N);
  fok = imv_static(fok,N,N);
  evec = cmv_static(evec,N,N);
  evec_inv = cmv_static(evec_inv,N,N);
  omega = cmv_static(omega,N,N);
  evalv = cmv_static(evalv,N,1);
  
  if (na >= NAMAX) {  
    fprintf(stderr,"static pointer allocation overflow (NAMAX)");
    exit(1);
  }
  if (first) {
    for (i=0;i<NAMAX;i++) {
      qt[i]=NULL;
      rhot[i]=NULL;
      favg[i]=NULL;
    }
    first=0;  
  }

  for (i=0;i<=na;i++) {
    qt[i] = cmv_static(qt[i],N,N);
    cmv_zero(qt[i]);
    rhot[i] = cmv_static(rhot[i],N,N);
    cmv_zero(rhot[i]);
    favg[i] = cmv_static(favg[i],N,N);
    cmv_zero(favg[i]);
  }

  /* DIAGONALIZE TOTAL PROPAGATOR uc */
  
  cmv_copy(uc,ustore[na]);
  cm_diag(uc,evalv,evec);
  /* ZT: for propagator made of Hamiltonian, the inverse transformation is the adjoint */
  cm_adjoint(evec_inv,evec);

  /* CREATE TRANSITION-FREQUENCY MATRIX omega[][] */
  
  Itr=Complx(0.0,-tr);

  for (i=0;i<N;i++) {
    complx l1 = Clog(evalv->data[i]);
    for (j=0;j<N;j++) {
      complx l2 = Clog(evalv->data[j]);
      omega->data[i+j*N] = Cdiv( Csub(l1,l2), Itr);
    }
  }
  
  /* CREATE TRANSFORMED Q-MATRICES
     Q^T0=X^-1 Q X
     Q^Tj = X^-1 Aj^+ Q Aj X
  */
  
  /* Q^T0=X^-1 Q X */
  cmv_copy(qt[0],obs);
  simtrans_adj(qt[0],evec);
  /* Q^Tj = X^-1 Aj^+ Q Aj X */
  for (ct=1; ct<na; ct++) {
     cmv_copy(qt[ct],obs);
     simtrans_adj(qt[ct],ustore[ct]);
     simtrans_adj(qt[ct],evec);
  }
  

  /* CREATE TRANSFORMED rho FOR ALL p */
  
  for (p=0;p<na;p++) {
  
    if (rhosymmetry) { /* rhoT=0.5*(QT + QT+) */
      /*      adjoint(tmp,qt[p]); */
      /*      rhot[p]=qt[p]+tmp; */
       cm_adjoint(rhot[p],qt[p]);
       cmv_addto(rhot[p],qt[p]);
    } else {

      /* rhoT0 = X^-1 rho X
         rhoTp = X^-1 Ap+ rho Ap X */

      if (p==0) {
         cmv_copy(rhot[p],rho);
	 simtrans_adj(rhot[p],evec);
      } else {
         cmv_copy(rhot[p],rho);
	 simtrans_adj(rhot[p],ustore[p]);
	 simtrans_adj(rhot[p],evec);
      }         
    } /* end if symmetry */
  } /* gamma loop 1 */
  
  
  /* CREATE GAMMAAVERAGED f-MATRICES "the hard way" */

  if (!realspec) {
      for (p=0;p<na;p++) {
        complx *rhot_p = rhot[p]->data;
        for (j=0;j<na;j++) {
          complx cc,*favg_j=favg[j]->data;
          complx *qt_ja;
	  ja= ((j+p) % na);
          if ((j+p) >= na) jb=j-na; else jb=j;
          cc.re= 0.0;
          cc.im= -(double)jb*tt;
          qt_ja = qt[ja]->data;  
          for (m=0;m<N;m++) {
            for (n=0;n<N;n++) {
              double f;
              complx ll,pp;
              complx *omega_mn = &(omega->data[m+n*N]);
              complx *qt_jamn = &(qt_ja[m+n*N]);
              complx *favg_jmn = &(favg_j[m+n*N]);
              complx *rhot_pnm = &(rhot_p[n+m*N]);

              ll.re = omega_mn->re * cc.re - omega_mn->im * cc.im;
              ll.im = omega_mn->im * cc.re + omega_mn->re * cc.im;
              f=exp(ll.re);
              pp.re = f*cos(ll.im);
              pp.im = f*sin(ll.im);

              ll.re = rhot_pnm->re * pp.re - rhot_pnm->im * pp.im;
              ll.im = rhot_pnm->im * pp.re + rhot_pnm->re * pp.im;

              favg_jmn->re += qt_jamn->re * ll.re - qt_jamn->im * ll.im;
              favg_jmn->im += qt_jamn->im * ll.re + qt_jamn->re * ll.im;
            }
        }
      }
    }
    
  } else { /* realspec assumed ; only nareal f-matrices needed; */

    for (m=0;m<N;m++) {
      for (n=0;n<N;n++) {
        complx *favg_0mn = &(favg[0]->data[m+n*N]);
	int *fok_mn = &(fok->data[m+n*N]);
        for (p=0;p<na;p++) {
          *favg_0mn = Cadd(*favg_0mn, Cmul(qt[p]->data[m+n*N],rhot[p]->data[n+m*N]));
        }
        /*	favg[na]=favg[0]; */
	*fok_mn = !is_small(*favg_0mn);
        if (*fok_mn) {
          for (j=1;j<=nareal;j++) {
            for (p=0;p<na;p++) {
              double f;
              complx ll,pp,cc;
              complx *omega_mn= &(omega->data[m+n*N]);
              complx *qt_jamn;
              complx *favg_jmn= &(favg[j]->data[m+n*N]);
              complx *rhot_pnm = &(rhot[p]->data[n+m*N]);

              ja= ((j+p) % na);
              if ((j+p) >= na) jb=j-na; else jb=j;

              cc.re=0.0;
              cc.im= -(double)jb*tt;
              qt_jamn= &(qt[ja]->data[m+n*N]);

              ll.re = omega_mn->re * cc.re - omega_mn->im * cc.im;
              ll.im = omega_mn->im * cc.re + omega_mn->re * cc.im;
              f=exp(ll.re);
              pp.re = f*cos(ll.im);
              pp.im = f*sin(ll.im);

              ll.re = rhot_pnm->re * pp.re - rhot_pnm->im * pp.im;
              ll.im = rhot_pnm->im * pp.re + rhot_pnm->re * pp.im;

              favg_jmn->re += qt_jamn->re * ll.re - qt_jamn->im * ll.im;
              favg_jmn->im += qt_jamn->im * ll.re + qt_jamn->re * ll.im;
            }
          }
        }
      }
    }


    if ((na%2)==0) {
      for (m=0;m<N;m++)
        for (n=0;n<N;n++)
          for (j=1;j<nareal;j++)
            favg[na-j]->data[m+n*N] = Conj(favg[j]->data[m+n*N]);
    } else {
      for (m=0;m<N;m++)
        for (n=0;n<N;n++)
          for (j=1;j<(nareal-1);j++)
            favg[na-j]->data[m+n*N] = Conj(favg[j]->data[m+n*N]);
    }

  }

  /* GENERATE THE gamma-averaged FID ...... */
    
  /*  slow but transparent... 
  
      for (i=0;i<nt;i++){
      j = (i % na);
      for (m=1;m<=N;m++)
      for (n=1;n<=N;n++) {
      fid[i+1] += favg[j][m][n]*exp(omega[m][n]*complx(0.0,double(i)*tt));
      }
      }
  */  
  
  /* FAST !
     tmp is used as temp variable.
    omega as the sum variable.
  */
  for (k=0;k<N;k++)
    for (l=0;l<N;l++)
      omega->data[k+l*N] = tmp->data[k+l*N] = Cexp(Cmul(omega->data[k+l*N],Complx(0.0,tt)));
  
  for (k=0;k<N;k++)
    for (l=0;l<N;l++)
      fid[1] = Cadd(fid[1],favg[0]->data[k+l*N]);

  /* this is *the* optimization */
  for (i=1;i<nt;i++){
    j = i % na;
    for (k=0;k<N;k++)
      for (l=0;l<N;l++) {
        double okl_re;
        complx *okl = &(omega->data[k+l*N]);
        complx *fi = &fid[i+1];
        complx *fjkl = &(favg[j]->data[k+l*N]);
        complx *tkl = &(tmp->data[k+l*N]);
        fi->re += fjkl->re * okl->re - fjkl->im * okl->im;
        fi->im += fjkl->im * okl->re + fjkl->re * okl->im;
        okl_re=okl->re;
        okl->re = okl_re*tkl->re - okl->im*tkl->im;
        okl->im = okl->im*tkl->re + okl_re*tkl->im;
      }
  }


  fac=1.0/(double)na;
  if (rhosymmetry) fac *= 0.5;
  k=LEN(fid);
  for (i=1;i<=k;i++) {
    fid[i].re *= fac;
    fid[i].im *= fac;
  }
}



int mul_number(int N,int na,int rep,int nt)
{
  int i,j,k,p,ii,jj;
  int M,Nq,N2,N3,Ntrace,Nmul,Nsim_trans,Nsim_trans_adj_Q;
  int narep = na * rep;

  N2=N*N;
  N3=N2*N;
  Nq=m_getnqlist();

  Ntrace = N2;
  Nmul = N3;
  Nsim_trans = 2*N3;
  Nsim_trans_adj_Q = Nq*N*(N+1);

  M=0;

  k=1;
  for (i=1;i<=na;i++) {      
    M += Nsim_trans_adj_Q;
    k++;
  }

  for (j=2;j<=rep;j++) {
    for (i=1;i<=na;i++) {
      M += Nsim_trans_adj_Q + Nmul;
      k++;
    }
    M += Nmul; 
  }
  M += N2;
  for (i=1;i<nt;i++) {
    ii=((i-1) % narep)+1;
    M += Ntrace;  
    if (ii == narep) {
      M += Nsim_trans;
    }
  }

  for (p=1;p<na;p++) {
    i=2;
    for (k=p;k<na;k++){
      M += Nsim_trans + Nmul + Ntrace;
      i++;
    }
    ii=1;
    for (i=na-p+2;i<=nt;i++,ii++){
      jj=(ii-1)%narep+1;
      M += Ntrace;
      if (jj==narep){
        M += Nsim_trans;
      }
    }
  } /* p (gamma) loop */
  return M;
}


int rep_minimize(int N,int na,int nt)
{
  int m0,m,r,r0=0,n;
  m0= 99999999L;
  n = nt/na-1;
  for (r=1;r<=n;r++) {
    m=mul_number(N,na,r,nt);
    if (m < m0) {
      r0=r;
      m0=m;      
    } else {
      break;
    }
  }
  return (int)r0;
}

int rep_minimize_estimate(int nsig,int nsampr,int N,int nqlist)
{
  int i;
  /* N : matrix dimension
     minimize the equation :
       2 N^3 nsig/(x nsampr) + x nsampr nqlist(N^2+N)
     with respect to x 
  */
  i= (int) (sqrt(
		     (double)(nsig*2*N*N*N)
		     /
		     (double)(nsampr*nsampr*nqlist*(N*N+N))
		     ));
  return (i > 0 ? i : 1);

  /* original: nsig/(x nsampr) + x nsampr
     int i= (int) (sqrt((double)nsig)/(double)nsampr);
     return (i > 0 ? i : 1); */
}


#define NAREPMAX 4096

void gammarep(mv_complx** vn,mv_complx** un,
  int na,mv_complx* rho,mv_complx* obs,double dummy,
  complx* fid,int rr)
{
  static int rep = 0;
  int i,ii,j,jj,k,p;
  int nt,N,narep;
  double f;
  static mv_complx *tmp=NULL,*us=NULL,*rhoa=NULL,*kn[NAREPMAX],*utmp[NAREPMAX],*vjj=NULL;
  static int first=1;

  nt=LEN(fid);
  N = rho->row;

/*  rr=-1; */
  if (!rep) {
    if (rr <= 0) {
      rep=rep_minimize(N,na,nt);
    } else {
      i= nt/na-1;
      rep=(rr > i ? i : rr);
    }
  }

  narep = na * rep;

  if (narep >= NAREPMAX) {
    fprintf(stderr,"errro: overflow in gamma_rep_fast (increase NAREPMAX)");
    exit(1);
  }
  
  if (first) {
    for (i=0;i<NAREPMAX;i++) {    
      kn[i]=NULL;
      utmp[i]=NULL;
    }
    first=0;
  }

  tmp = cmv_static(tmp,N,N);
  us = cmv_static(us,N,N);
  rhoa = cmv_static(rhoa,N,N);
  vjj = cmv_static(vjj,N,N);

  for (i=1;i<=narep;i++) {
    kn[i] = cmv_static(kn[i],N,N);
    utmp[i] = cmv_static(utmp[i],N,N);
  }

  k=1;
  for (i=1;i<=na;i++) {      
    cmv_copy(utmp[i],vn[i]);
    m_simtransadjQ(kn[k],vn[i]);
    k++;
  }
  
  cmv_copy(us,vn[na]);
  for (j=2;j<=rep;j++) {
    for (i=1;i<=na;i++) {
      cmv_copy(tmp,utmp[i]);
      cm_mul(utmp[i],tmp,vn[na]);
      m_simtransadjQ(kn[k],utmp[i]);
      k++;
    }
    cmv_copy(tmp,us);
    cm_mul(us,tmp,vn[na]);
  }
  cmv_copy(rhoa,rho);
  fid[1]=Cmul(cm_trace(obs,rhoa),Complx(na,0)); /* *na=ngamma  */

  for (i=1;i<nt;i++) {
    ii=((i-1) % narep)+1;
    fid[i+1]=cm_trace(kn[ii],rhoa);
    if (ii == narep) {
      simtrans(rhoa,us);
    }
  }

  /* INTELLIGENT GAMMA-AVERAGING STARTS HERE */

  for (p=1;p<na;p++) {

    cm_unit(vjj);

    i=2;
    for (k=p;k<na;k++){
      cmv_copy(rhoa,rho);
      cmv_copy(tmp,vjj);
      cm_mul(vjj,un[k+1],tmp);        /* vjj=U[k+1]*vjj */

      simtrans(rhoa,vjj);

      fid[i]=Cadd(fid[i],cm_trace(obs,rhoa));
      i++;
    }
    ii=1;

    for (i=na-p+2;i<=nt;i++,ii++){
      jj=(ii-1)%narep+1;
      fid[i]=Cadd(fid[i],cm_trace(kn[jj],rhoa));
      if (jj==narep){
	simtrans(rhoa,us);
      }
    }
  } /* p (gamma) loop */
  f= 1.0/(double)na;
  for (i=1;i<=nt;i++) {
     fid[i].re *= f;
     fid[i].im *= f;
  }
}

  
  




















