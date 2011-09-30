/*
    Complex matrix blockdiagonalization
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
    
*/

#include <stdio.h>
#include <stdlib.h>
#include "cm_new.h"
#include "defs.h"


typedef struct _Qlist {
  int i,j;
  complx e;
} Qlist;

Qlist* qlist=NULL;
int nqlist;

int m_getnqlist() { return nqlist; }

void m_destroyQlist() {
  free(qlist);
  qlist=NULL;
}

void m_makeQlist(mv_complx* Q)
{
  int i,j,N;
  
  N = Q->row;
  if (qlist != NULL) {
    /* fprintf(stderr,"warning: Qlist was already initalized\n"); */
    m_destroyQlist();
  }
  nqlist=0;
  for (i=0;i<N;i++) {
    for (j=0;j<N;j++) {
      if (!is_small(Q->data[i+j*N])) {
         nqlist++;
      }
    }
  }
  qlist=(Qlist*)malloc(sizeof(Qlist)*(nqlist+1));
  if (!qlist) {
    fprintf(stderr,"error: allocation failure in m_makeQList\n");
    exit(1);
  }  
  nqlist=0;

  for (i=0;i<N;i++) {
    for (j=0;j<N;j++) {
      if (!is_small(Q->data[i+j*N])) {
        nqlist++;
        qlist[nqlist].i=i;
        qlist[nqlist].j=j;
        qlist[nqlist].e=Q->data[i+j*N];
      }
    }
  }  
}

void m_simtransadjQ(mv_complx * L,mv_complx* U)
{
  int f,i,j,k,l,N;  
  complx e,x;
  
  /* L = U+  * Q * U
     Call m_makeQlist with the argument Q, first. */

  N = U->row;

  cmv_zero(L);
  for (f=1;f<=nqlist;f++) {
    i=qlist[f].i;
    j=qlist[f].j;
    e=qlist[f].e;
    for (k=1;k<=N;k++) {
      x=Cmul(Conj(U->data[i+k*N]),e);
      for (l=1;l<=N;l++) {
	 L->data[k+l*N] = Cadd( L->data[k+l*N], Cmul(x,U->data[j+l*N]));
      }
    }
  }
}

/* BLOCK DIAGONALIZED OPTIMIZATIONS */
/* ZT: removed...    */

/*
void (*PDmul)(complx** ab,complx** a, complx** b) = m_mul;
void (*PDrealexp)(complx** m,complx** c,double dt)= m_realexp;
void (*Pmul)(complx** ab, complx** a, complx** b) = m_mul;
complx (*Ptrace)(complx** a, complx** b)            = m_trace;
void (*Psimtrans)(complx** S, complx** U)         = m_simtrans;

void m_setmatrixtype_block(double* blkdiag,int N)
{
  m_set_blk(blkdiag,N);
  PDmul=m_mul_blk;
  PDrealexp=m_realexp_blk;
  Pmul=m_mul_blk;
  Ptrace=m_trace_blk;
  Psimtrans=m_simtrans_blk;
}

void m_setmatrixtype_normal()
{
  PDmul=m_mul;
  PDrealexp=m_realexp;
  Pmul=m_mul;
  Ptrace=m_trace;
  Psimtrans=m_simtrans;
}
*/
