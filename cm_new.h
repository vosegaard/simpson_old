/*
    Complex matrix operations using BLAS and LAPACK
    Copyright (C) 2009 Zdenek Tosner

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

#ifndef __CM_H
#define __CM_H

#ifndef __MATRIX_H
#include "matrix_new.h"
#endif

#define SMALL 1.0e-10

#define is_small(c) (!((c).re > SMALL || (c).re < -SMALL  || (c).im > SMALL || (c).im < -SMALL))



/* useful constants defined as global variables */
extern const complx CPLX1, CPLXm1, CPLX0;
extern const int INTONE;

mv_complx * cmv_dup(mv_complx * obj);

complx * get_complx_wsp(int len, complx * wsp, int * wsp_size, int * wsp_id);
double * get_double_wsp(int len, double * wsp, int * wsp_size, int * wsp_id);
int * get_int_wsp(int len, int * wsp, int * wsp_size, int * wsp_id);
void free_wsp();
mv_complx * cmv_static(mv_complx * obj, int row, int col);
mv_double * dmv_static(mv_double * obj, int row, int col);
mv_int * imv_static(mv_int * obj, int row, int col);
void free_mv_static();

void cmv_real(mv_complx *C, mv_double *D);
void cmv_imag(mv_complx *C, mv_double *D);
void dmv_complx(mv_double *D, mv_complx *C);
int cmv_isreal(mv_complx *C);
double * get_real_diag(mv_complx *C);
double * get_real(mv_complx *C);
int cm_ishermit(mv_complx *a);
void cv_conj(mv_complx *conja, mv_complx * a);
void cv_conj_in(mv_complx * a);
void cm_conj(mv_complx *conja, mv_complx * a);
double cmv_sumnorm1(mv_complx *tmp);

complx cm_trace(mv_complx * a,mv_complx * b);
complx cm_trace_adjoint(mv_complx * a,mv_complx * b);
void cm_and(mv_complx * a_and_b, mv_complx * a,mv_complx * b);
void cm_or(mv_complx * a_and_b, mv_complx * a,mv_complx * b);

void cmv_zero(mv_complx * obj); 
void dmv_zero(mv_double * obj); 
mv_complx * cm_direct(mv_complx * a,mv_complx * b);
mv_complx * cm_directadd(mv_complx * a,mv_complx * b);
void cm_copydiagmv(mv_complx * diagm,mv_complx * v);
void cmv_copy(mv_complx * dst, mv_complx * src);
void cm_unit(mv_complx * obj);
void cv_d1(mv_complx *obj);

void cmv_mulc(mv_complx * obj, complx z);
void cmv_muld(mv_complx * obj, double d);
void cmv_multoc(mv_complx * y, mv_complx * x, complx z); /* y = y + z*x */
void cmv_multod(mv_complx * y, mv_complx * x, double d);
/* this simply multiplies two objects element wise */ 
void cmv_mul_elem(mv_complx * res, mv_complx * a, mv_complx * b);
void cm_mul_fd(mv_complx * res, mv_complx * a, mv_complx * b); /* a is full, b is diagonal */
void cm_mul_df(mv_complx * res, mv_complx * a, mv_complx * b); /* a is diagonal, b is full */

void cm_mul(mv_complx * res, mv_complx * a, mv_complx * b);
void cm_commutator(mv_complx * res, mv_complx * a, mv_complx * b);

complx cv_dotmul(mv_complx * a, mv_complx * b);
void cm_mulmv(mv_complx * y, mv_complx * A, mv_complx * x);
void cm_mulvm(mv_complx * y, mv_complx * x, mv_complx * A);
void cm_adjoint(mv_complx *aplus, mv_complx * a);

void cm_print(mv_complx * m,char* name);
void dm_print(mv_double * m,char* name);
void cv_print(mv_complx * v,char* name);
void dv_print(mv_double * v,char* name);

  /* S = U  * S * U+  where U is diagonal, two versions */
void simtrans_diagprop1(mv_complx * A, mv_complx * U);
void simtrans_diagprop2(mv_complx * A, mv_complx * U);
void simtrans_diagprop3(mv_complx * A, mv_complx * U);
void simtrans_adj_diagprop3(mv_complx * A, mv_complx * U);
void simtrans_zrot(mv_complx * U, double * P);
void simtrans_zrot2(mv_complx * U, double * P);
void ham_zrot_real(mv_complx *Ham, double * P);
void simtrans(mv_complx * S, mv_complx * U);
void simtransh(mv_complx * H, mv_complx * U);
void simtrans_adj(mv_complx * S, mv_complx * U);
void simtransh_adj(mv_complx * S, mv_complx * U);

void prop_realdiag(mv_complx * res, mv_double * dv, double dt);
/* propagator using diagonalization of real Hamiltonian */
void prop_real(mv_complx * prop, mv_double * ham, double dt); 
/* propagator using Chebyshev approx. and complex Hamiltonian */
void prop_cheb(mv_complx * prop, mv_complx * ham, double dt); 
/* propagator using Pade approx. + Sc.&Sq., real Hamiltonian */
void prop_pade_real(mv_complx * prop, mv_double * ham, double dt); 
/* propagator using Pade approx. + Sc.&Sq., complex Hamiltonian */
void prop_pade_complx(mv_complx * prop, mv_complx * ham, double dt);

void cm_multo(mv_complx * ab, mv_complx * b);
void cm_multo_rev(mv_complx * ba, mv_complx * b);
void cmv_addto(mv_complx *ab, mv_complx * b);
void cmv_subfrom(mv_complx *ab, mv_complx * b);

mv_complx * cm_ln(mv_complx * m);
void dm_symdiag(mv_double *A, mv_double *eigs, mv_double *T);
void cm_diag(mv_complx *A, mv_complx *eigs, mv_complx *T);


#include "cmblock.h"

#endif /* __CM_H */
