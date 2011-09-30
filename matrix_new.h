/*
    Matrix and vector definition
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
  
    Matrix definition and implementation. 
*/

#ifndef __MATRIX_H
#define __MATRIX_H

/* AB::beg */
/* Define struct for the Euler angles used for powder averaging.  */
/* This allows the use of 2 or 3 angle sets for powder averaging. */
/* In case of 2 angle sets gamma is set to zero.                  */ 
/*typedef struct _Omega {
	double alpha, beta, gamma, weight;
} Omega; */
/* AB::end */

#ifndef __COMPLX_H
#include "complx.h"
#endif


#ifndef M_PI
#define M_PI            3.14159265358979323846
#endif
 

typedef struct _mv_complx {
       int row, col;
       complx * data;
} mv_complx;

typedef struct _mv_double {
       int row, col;
       double * data;
} mv_double;

typedef struct _mv_float {
       int row, col;
       float * data;
} mv_float;

typedef struct _mv_int {
       int row, col;
       int * data;
} mv_int;

/***
 * complex 
 ***/
mv_complx * complx_matrix_alloc(int rows, int cols);
void complx_matrix_free(mv_complx * obj);
mv_complx * complx_vector_alloc(int len);
void complx_vector_free(mv_complx * obj);
/* is it good idea to have this 1-based??? */
void complx_put_elem(mv_complx * obj, int row, int col, complx z);
/* is it good idea to have this 1-based??? */
complx complx_get_elem(mv_complx * obj, int row, int col);
/***
 * double 
 ***/
mv_double * double_matrix_alloc(int rows, int cols);
void double_matrix_free(mv_double * obj);
mv_double * double_vector_alloc(int len);
void double_vector_free(mv_double * obj);
/* is it good idea to have this 1-based??? */
void double_put_elem(mv_double * obj, int row, int col, double z);
/* is it good idea to have this 1-based??? */
double double_get_elem(mv_double * obj, int row, int col);
double * dm_row(mv_double *obj, int row);
double * dm_col(mv_double *obj, int col);


/***
 * float is missing...
 ***/

/***
 * int 
 ***/
mv_int * int_matrix_alloc(int rows, int cols);
void int_matrix_free(mv_int * obj);
mv_int * int_vector_alloc(int len);
void int_vector_free(mv_int * obj);
/* is it good idea to have this 1-based??? */
void int_put_elem(mv_int * obj, int row, int col, int z);
/* is it good idea to have this 1-based??? */
int int_get_elem(mv_int * obj, int row, int col);

/* this I keep for ease */
#define LEN(v)   (*(int*)v)
#define ROWS(m)  (*((int*)*m))
int * int_vector(int len);
double * double_vector(int len);
complx * complx_vector(int len);
complx** complx_matrix(int nrow, int ncol);
void free_int_vector(int *);
void free_double_vector(double *);
void free_complx_vector(complx *);
void free_complx_matrix(complx **);
void m_zerov(complx * v);

/* AB::beg */
/* Omega is the struct containing the Euler angles for powder averaging
PROTO_VECTOR_MATRIX_TYPE(Omega)
   AB::end */

#endif /* __MATRIX_H */

