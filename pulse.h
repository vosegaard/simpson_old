/*
    Pulse propagation
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

#ifndef __PULSE_H
#define __PULSE_H

#include <tcl.h>
#include "ham.h"
#include "spinsys.h"

#define MAXMATRIX 1023
#define MAXSTO  1023
#define MAXSPINS 100

typedef struct _Pulse {

  int    dor,N,nspins,nchan;
  double phv[30];
  double rfv[30];
  double dtmax,t_usec,wr,sw,dwellt,brl,brlorig,dt,brl1,brl2,gamma1,gamma2;
  double tpropstart_usec,tpulsestart_usec,gamma_add,wr1,wr2;
  char pulsename[64];
  double* rfscalefactors; 
    
  int curr_nsig;
  complx *fid;
  mv_complx *sigma, *fdetect, *fstart;
  int acq_adjoint, is_relax;

  Hamilton* H;
  SpinSys *ss;
  
  double *chan_Ix[30], *chan_Iz[30], *Ix[100], *Iz[100], *tmpIz;
  complx *chan_Iy[30];
  double *sumUph, *sumHrf;
  mv_double * Htot;
  
  mv_complx *U, *dU, *tmpU;
  
  mv_complx *STO[MAXSTO+1];
  mv_complx *matrix[MAXMATRIX+1];
  double   STO_tproplength_usec[MAXSTO+1];
  double   STO_tpropstart_usec[MAXSTO+1];
  double   STO_brl[MAXSTO+1];
  int Uisunit,waspulse,STO_waspulse[MAXSTO+1];
  int isselectivepulse,cannotbestored,hamchanged,spinused[MAXSPINS+1];
  Tcl_Interp* interp;
  
  double zcoor;
  double *inhom_offset;
  int propmethod;
  /* ZT: tricks for new pulseq scheme */
  int Npsq, *acqnp, cacqs_pos, new_direct_mxpos, new_gcompute;
  double *acqph, check_dwelltime_t0;
  mv_int *pacqs;
  Tcl_Obj **psq, **cacqs;
  
} Pulse;


void pulse_initialize(Pulse* P,Tcl_Interp* interp,SpinSys* ss,double sw,double wr,double brl,int dor,double wr1,double brl1,double wr2,double brl2);
void pulse_setpulsename(Pulse* P,char* pulsename);
void pulse_setprop(Pulse* P,mv_complx * fstart,mv_complx * fdetect);
void pulse_propagate(Pulse* P,Hamilton* H,mv_complx * U,double dt,double t,complx* fid);

/*
void pulse_propagatemultiacq(Pulse* P,Hamilton* H,complx*** U,double dt,double _t);
*/
void pulse_destroy(Pulse* P);

/********* ZT: few extra definitions for OC to work ************/
#define check_pulse() if (puls == NULL) return TclError(interp,"error: pulse sequence not initiated\n");

void _evolve_with_prop(void);
void _reset_prop(void);
void _ph(int channel,double phase);
void _rf(int channel,double rffield);
void _pulse(double duration);

#endif /* __Pulse_H */
