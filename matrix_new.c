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
  
    Matrix definition and implementation. Data structure according to BLAS, 
    i.e. Column major. Vector is just a single column. 
*/

#include <stdlib.h>
#include <stdio.h>
#include "complx.h"
#include "matrix_new.h"

/* this I keep for ease...
         1-based vectors with the length saved in the 0-element */
int * int_vector(int len)
{
   int * v;
   
   v = (int*)malloc((len+1)*sizeof(int));
   if (!v) {
      fprintf(stderr,"Can not allocate int vector (%i)\n",len);
      exit(1);
   }
   *(int*)v = len;
   return v;
}

double * double_vector(int len)
{
   double * v;
   
   v = (double*)malloc((len+1)*sizeof(double));
   if (!v) {
      fprintf(stderr,"Can not allocate double vector (%i)\n",len);
      exit(1);
   }
   *(int*)v = len;
   return v;
}

complx * complx_vector(int len)
{
   complx * v;
   
   v = (complx*)malloc((len+1)*sizeof(complx));
   if (!v) {
      fprintf(stderr,"Can not allocate complx vector (%i)\n",len);
      exit(1);
   }
   *(int*)v = len;
   return v;
}

/* I would like to get rid of this but don't have time... */
complx ** complx_matrix(int nrow, int ncol)
{
   int i;
   complx **m;
   
   m = (complx**)malloc(((nrow+1)*sizeof(complx*)));
   if (!m) {
      fprintf(stderr,"old complx_matrix error: unable to allocate pointers\n");
      exit(1);
   }
   m[0] = (complx*)malloc(2*sizeof(int));
   ((int*)m[0])[0] = nrow;
   ((int*)m[0])[1] = ncol;
   m[1] = (complx*)malloc(((nrow*ncol+1)*sizeof(complx)));
   if (!m[1]) {
      fprintf(stderr,"old complx_matrix error: unable to allocate matrix\n");
      exit(1);
   }
   for (i=2; i<=nrow; i++) m[i] = m[i-1]+ncol;
   return m;
}

void free_complx_matrix(complx** m)
{
   free((char*)m[0]);
   free((char*)m[1]);
   free((char*)m);
}

void free_int_vector(int *v)
{
   free((char*)v);
}

void free_double_vector(double *v)
{
   free((char*)v);
}

void free_complx_vector(complx *v)
{
   free((char*)v);
}

void m_zerov(complx * v)
{
  memset(&v[1],0,LEN(v)*sizeof(complx));
}



/***
 * complex 
 ***/
mv_complx * complx_matrix_alloc(int rows, int cols)
{
   mv_complx * obj;
   
   obj = (mv_complx*)(malloc(sizeof(mv_complx)));
   obj->row = rows;
   obj->col = cols;
   obj->data = (complx*)(malloc(rows*cols*sizeof(complx)));
   if (!(obj->data)) {
      fprintf(stderr,"Can not allocate complex matrix (%i x %i)\n",rows, cols);
      exit(1);
   }
   return obj;
}

void complx_matrix_free(mv_complx * obj)
{
   free((char*)(obj->data));
   obj->row = obj->col = 0;
   free((char*)obj);
}

mv_complx * complx_vector_alloc(int len)
{
   mv_complx * obj;
   
   obj = (mv_complx*)(malloc(sizeof(mv_complx)));
   obj->col = 1;
   obj->row = len;
   obj->data = (complx*)(malloc(len*sizeof(complx)));
   if (!(obj->data)) {
      fprintf(stderr,"Can not allocate complex vector (%i)\n",len);
      exit(1);
   }
   return obj;
}
   
void complx_vector_free(mv_complx * obj)
{
   free((char*)(obj->data));
   obj->row = obj->col = 0;
   free((char*)obj);
}

/* is it good idea to have this 1-based??? */
void complx_put_elem(mv_complx * obj, int row, int col, complx z)
{
   int k;
   
   if ( (obj->row < row) || (obj->col < col) ) {
      fprintf(stderr,"complx_put_elem error: wrong element indeces\n");
      exit(1);
   }
   k = row-1 + (col-1)*(obj->row);
   obj->data[k].re = z.re;
   obj->data[k].im = z.im;
}

/* is it good idea to have this 1-based??? */
complx complx_get_elem(mv_complx * obj, int row, int col)
{
   complx z;
   
   if ( (obj->row < row) || (obj->col < col) ) {
      fprintf(stderr,"complx_get_elem error: wrong element indeces\n");
      exit(1);
   }
   z = obj->data[row-1 + (col-1)*(obj->row)];
   return z;
}


/***
 * double
 */
mv_double * double_matrix_alloc(int rows, int cols)
{
   mv_double * obj;
   
   obj = (mv_double*)(malloc(sizeof(mv_double)));
   obj->row = rows;
   obj->col = cols;
   obj->data = (double*)(malloc(rows*cols*sizeof(double)));
   if (!(obj->data)) {
      fprintf(stderr,"Can not allocate double matrix (%i x %i)\n",rows, cols);
      exit(1);
   }
   return obj;
}

void double_matrix_free(mv_double * obj)
{
   free((char*)(obj->data));
   obj->row = obj->col = 0;
   free((char*)obj);
}

mv_double * double_vector_alloc(int len)
{
  mv_double * obj;
  
   obj = (mv_double*)(malloc(sizeof(mv_double)));
   obj->col = 1;
   obj->row = len;
   obj->data = (double*)(malloc(len*sizeof(double)));
   if (!(obj->data)) {
      fprintf(stderr,"Can not allocate double vector (%i)\n",len);
      exit(1);
   }
   return obj;
}
   
void double_vector_free(mv_double * obj)
{
   free((char*)(obj->data));
   obj->row = obj->col = 0;
   free((char*)obj);
}

/* is it good idea to have this 1-based??? */
void double_put_elem(mv_double * obj, int row, int col, double z)
{
   int k;
   
   if ( (obj->row < row) || (obj->col < col) ) {
      fprintf(stderr,"double_put_elem error: wrong element indeces\n");
      exit(1);
   }
   k = row-1 + (col-1)*(obj->row);
   obj->data[k] = z;
}

/* is it good idea to have this 1-based??? */
double double_get_elem(mv_double * obj, int row, int col)
{
   double z;
   
   if ( (obj->row < row) || (obj->col < col) ) {
      fprintf(stderr,"double_get_elem error: wrong element indeces\n");
      exit(1);
   }
   z = obj->data[row-1 + (col-1)*(obj->row)];
   return z;
}

/* is it good idea to have this 1-based??? */
double * dm_row(mv_double *obj, int row)
{
   double *res, *sl1, *sl2;
   int i, Nr, Nc;
   
   Nr = obj->row;
   if ( row > Nr ) {
      fprintf(stderr,"dm_row error: row out of obj. size\n");
      exit(1);
   } 
   Nc = obj->col;
   res = (double*)malloc(Nc*sizeof(double));
   sl1 = res;
   sl2 = obj->data + row-1;
   for (i=0; i<Nc; i++) {
      *sl1 = *sl2;
      sl1++;
      sl2 += Nr;
   }
   return res;
}

/* is it good idea to have this 1-based??? */
double * dm_col(mv_double *obj, int col)
{
   double *res;
   int N;
   
   if (col > obj->col) {
      fprintf(stderr,"dm_col error: col out of obj. range\n");
      exit(1);
   }
   N = obj->row;
   res = (double*)malloc(N*sizeof(double));
   memcpy(res,obj->data+N*(col-1),N*sizeof(double));
   return res;
}




/***
 * float is missing here...
 ***/


/***
 * int
 */
mv_int * int_matrix_alloc(int rows, int cols)
{
   mv_int * obj;
   
   obj = (mv_int*)(malloc(sizeof(mv_int)));
   obj->row = rows;
   obj->col = cols;
   obj->data = (int*)(malloc(rows*cols*sizeof(int)));
   if (!(obj->data)) {
      fprintf(stderr,"Can not allocate int matrix (%i x %i)\n",rows, cols);
      exit(1);
   }
   return obj;
}

void int_matrix_free(mv_int * obj)
{
   free((char*)(obj->data));
   obj->row = obj->col = 0;
   free((char*)obj);
}

mv_int * int_vector_alloc(int len)
{
   mv_int * obj;
   
   obj = (mv_int*)(malloc(sizeof(mv_int)));
   obj->col = 1;
   obj->row = len;
   obj->data = (int*)(malloc(len*sizeof(int)));
   if (!(obj->data)) {
      fprintf(stderr,"Can not allocate int vector (%i)\n",len);
      exit(1);
   }
   return obj;
}
   
void int_vector_free(mv_int * obj)
{
   free((char*)(obj->data));
   obj->row = obj->col = 0;
   free((char*)obj);
}

/* is it good idea to have this 1-based??? */
void int_put_elem(mv_int * obj, int row, int col, int z)
{
   int k;
   
   if ( (obj->row < row) || (obj->col < col) ) {
      fprintf(stderr,"int_put_elem error: wrong element indeces\n");
      exit(1);
   }
   k = row-1 + (col-1)*(obj->row);
   obj->data[k] = z;
}

/* is it good idea to have this 1-based??? */
int int_get_elem(mv_int * obj, int row, int col)
{
   int z;
   
   if ( (obj->row < row) || (obj->col < col) ) {
      fprintf(stderr,"int_get_elem error: wrong element indeces\n");
      exit(1);
   }
   z = obj->data[row-1 + (col-1)*(obj->row)];
   return z;
}
