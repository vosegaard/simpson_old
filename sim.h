/*
    Simulation setup and calculation
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
    
*/

#ifndef __SIM_H
#define __SIM_H

#include <tcl.h>
#include "cryst.h"
#include "ham.h"
#include "spinsys.h"
#include "pulse.h"

typedef struct _Sim {
  Pulse* P;
  Hamilton* H;
  SpinSys* ss;

  char method[32],crystfile[256],rfproffile[256];
  int ni,np,imethod,gammethod,nstepr,nsampr_rep;
  int ng,ns,ntot,ngamma,matdim,conjugate_fid,obs_nuc;
  int rhosymmetry,realspec,dipole_check,dor;
/* AB::beg */
/* Flag to determine, whether a 2 or 3 angle set is used for powder averaging. */
/* I guess in the ideal case this should be detected automatically from the    */
/* angle file. However this is not implemented yet.                            */
  int use_3_angle_set;
/* AB::end */
  double tstepr,wr,brl,sw,specfreq,sw1,gamma_zero,wr1,wr2,brl1,brl2,gamma1,gamma2;
  double* blkdiag;
  mv_complx *fstart,*fdetect;
  
  /* Temporary variables, which are too big to be allocated
     each time the 'calcfid' procedure is called */
  mv_complx **un,**vn,*vi,*vj,*ku;
  complx* fid;
  mv_complx * sigma;
  /* ZT: tricks for new direct method */
  double dw, dw1, taur;

} Sim;

void sim_initialize(Tcl_Interp* interp,Sim* sim);
void sim_destroy(Sim* sim);
/*int sim_calcfid(Sim* sim,double alpha,double beta,double* rfscalefact,complx* fidsum); */
int sim_calcfid(Sim* sim,Omega omega,double* rfscalefact, complx* fidsum);


#define    M_GAMMAREP        2000
#define    M_GCOMPUTE        2001
#define    M_DIRECT          2002
#define    M_DIRECT_NEW      2009
#define    M_GCOMPUTE_NEW    2010
#define    M_GCOMPUTE2_NEW   2011


#endif /* __SIM_H */
