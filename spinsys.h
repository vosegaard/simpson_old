/*
    Spinsystem declaration and routines
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

#ifndef __SPINSYS_H
#define __SPINSYS_H

#ifndef __COMPLX_H
#include "complx.h"
#endif
#ifndef __MATRIX_H
#include "matrix_new.h"
#endif

#define MAXCHAN  30
#define MAXSPINS 100

typedef struct _ISOTOPE {
  int  number;
  char name[8];
  double spin;
  double gamma;
} ISOTOPE;     

extern ISOTOPE isotopes[];

typedef struct _SpinSys { 
  ISOTOPE* iso[MAXSPINS+1];
  int nspins;
  int matdim;
  
  int chan[MAXCHAN][32];
  int nchanelem[MAXCHAN];
  char channames[MAXCHAN][8];
  int nchan;
} SpinSys;

void ss_showspins(SpinSys* S);

ISOTOPE* ss_findisotope(char* name);
double ss_qn(SpinSys* S,int spin);
double ss_gamma1H();
double ss_gamma(SpinSys* S,int spin);
ISOTOPE* ss_isotope(SpinSys* S,int spin);
int ss_matdim(SpinSys* S);
void ss_addspin(SpinSys* S,char* name);
void ss_initialize(SpinSys* S);
mv_complx * ss_qdiag(SpinSys* S);
int ss_issame(SpinSys* S,int spin1,int spin2);

mv_complx * Iq(SpinSys* S);
mv_complx * Iqdelta(SpinSys* S);

mv_complx * Ie(SpinSys* S);
mv_complx * Ic(SpinSys* S,int spin);
mv_complx * Ip(SpinSys* S,int spin);
mv_complx * Im(SpinSys* S,int spin);
mv_complx * Ix(SpinSys* S,int spin);
mv_complx * Iy(SpinSys* S,int spin);
mv_complx * Iz(SpinSys* S,int spin);
mv_complx * Icoherence(SpinSys* S, double* coh);

mv_complx * ss_oper(SpinSys* S,char* name);

mv_complx * ss_readoper(SpinSys* S,char* sop);
mv_complx * ss_eq_sigma(SpinSys* S);

void ss_hashinit();
void ss_hashinsert(char* ptr,complx c);
int  ss_hashlookup(char* ptr,complx* c);
void ss_hashdestroy();

#endif


