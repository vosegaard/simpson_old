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
    
*/

#ifndef __FIDCALC_H
#define __FIDCALC_H

void gcompute(mv_complx** ustore,int na, mv_complx* rho, mv_complx* obs,double tt,complx* fid,int rhosymmetry,int realspec);
void gammarep(mv_complx** vn,mv_complx** un,int na, mv_complx* rho, mv_complx* obs,double,complx* fid,int rep);

int rep_minimize_estimate(int nsig,int nsampr,int N,int nqlist);
int rep_minimize(int N,int na,int nt);
int is_rhosymmetry(mv_complx* fstart,mv_complx* fdetect);


#endif
