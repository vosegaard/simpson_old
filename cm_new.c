/*
    Complex matrix allocation and calculations
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
    
    
    Sort of naming convention:
      cm_<name>     : complex matrix procedure
      cv_<name>     : complex vector procedure
      dm_<name>     : double matrix procedure
      dv_<name>     : double vector procedure
      
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "cm_new.h"

#ifdef MKL_INTEL
#include "mkl.h"
#elif defined(__APPLE__)
#include <Accelerate/Accelerate.h>
#else
#include "defs_blas_lapack.h"
#include <cblas.h>
#endif

#define DEBUGPRINT noprintf

void noprintf(const char* format, ...)
{
}

/* useful constants */
const complx CPLX1={1.0,0.0}, CPLXm1={-1.0,0.0}, CPLX0={0.0,0.0};
/*   CPLX1.re = 1.0; CPLX1.im = 0.0;
   CPLXm1.re = -1.0; CPLXm1 = 0.0;
   CPLX0.re = 0.0; CPLX0 = 0.0; */
const int INTONE=1;
   

#define MAX_STATIC_WSP 1024

int n_cwsp=0;
complx* ptr_cwsp[MAX_STATIC_WSP];

complx * get_complx_wsp(int len, complx* wsp, int * wsp_size, int * wsp_id)
{

  if (wsp != NULL) {
    if (wsp != ptr_cwsp[*wsp_id]) {
       fprintf(stderr,"get_complx_wsp error: mismatch in workspace identifiers!\n");
       exit(1);
    }
    if (*wsp_size == len) {
       DEBUGPRINT("get_complx_wsp does nothing, sizes are OK\n");
       return wsp;
    }
    DEBUGPRINT("get_complx_wsp re-alloc slot %i: block (%d)->(%d)\n",*wsp_id,*wsp_size,len);
    free((char*)wsp);
    wsp=(complx*)(malloc(len*sizeof(complx)));
    *wsp_size = len;
    ptr_cwsp[*wsp_id] = wsp;
  } else {
    n_cwsp++;
    if (n_cwsp >= MAX_STATIC_WSP) {
      fprintf(stderr,"overflow error in cm.c: increase MAX_STATIC_WSP\n");
      exit(1);
    }
    DEBUGPRINT("complx_wsp alloc (%d)-block\n",len);
    ptr_cwsp[n_cwsp]=(complx*)(malloc(len*sizeof(complx)));
    wsp = ptr_cwsp[n_cwsp];
    *wsp_size = len;
    *wsp_id = n_cwsp;
  }
  return wsp;
}

/* static mv_complx objects */
int n_cmv=0;
mv_complx *ptr_cmv[MAX_STATIC_WSP];

mv_complx * cmv_static(mv_complx * obj, int row, int col)
{
   int len;
   
   len = row*col;
   if (obj != NULL) {
      if ( (obj->row)*(obj->col) == len) {
         obj->row = row;
	 obj->col = col;
	 DEBUGPRINT("cmv_static: size match\n");
	 return obj;
      }
      free((char*)obj->data);
      obj->data = (complx*)malloc(len*sizeof(complx));
      obj->row = row;
      obj->col = col;
      DEBUGPRINT("cmv_static: re-allocating\n");
      return obj;
   }
   n_cmv++;
   if (n_cmv > MAX_STATIC_WSP) {
      fprintf(stderr,"overflow errorin cm.c: increase MAX_STATIC_WSP\n");
      exit(1);
   }
   ptr_cmv[n_cmv] = obj = complx_matrix_alloc(row,col);
   return obj;
}

/* static mv_double objects */
int n_dmv=0;
mv_double *ptr_dmv[MAX_STATIC_WSP];

mv_double * dmv_static(mv_double * obj, int row, int col)
{
   int len;
   
   len = row*col;
   if (obj != NULL) {
      if ( (obj->row)*(obj->col) == len) {
         obj->row = row;
	 obj->col = col;
	 return obj;
      }
      free((char*)obj->data);
      obj->data = (double*)malloc(len*sizeof(double));
      obj->row = row;
      obj->col = col;
      return obj;
   }
   n_dmv++;
   if (n_dmv > MAX_STATIC_WSP) {
      fprintf(stderr,"overflow errorin cm.c: increase MAX_STATIC_WSP\n");
      exit(1);
   }
   ptr_dmv[n_dmv] = obj = double_matrix_alloc(row,col);
   return obj;
}

/* static mv_int objects */
int n_imv=0;
mv_int *ptr_imv[MAX_STATIC_WSP];

mv_int * imv_static(mv_int * obj, int row, int col)
{
   int len;
   
   len = row*col;
   if (obj != NULL) {
      if ( (obj->row)*(obj->col) == len) {
         obj->row = row;
	 obj->col = col;
	 return obj;
      }
      free((char*)obj->data);
      obj->data = (int*)malloc(len*sizeof(int));
      obj->row = row;
      obj->col = col;
      return obj;
   }
   n_imv++;
   if (n_imv > MAX_STATIC_WSP) {
      fprintf(stderr,"overflow errorin cm.c: increase MAX_STATIC_WSP\n");
      exit(1);
   }
   ptr_imv[n_imv] = obj = int_matrix_alloc(row,col);
   return obj;
}


int n_dwsp=0;
double* ptr_dwsp[MAX_STATIC_WSP+1];

double * get_double_wsp(int len,double* wsp, int *wsp_size, int * wsp_id)
{
  if (wsp != NULL) {
    if (wsp != ptr_dwsp[*wsp_id]) {
       fprintf(stderr,"get_complx_wsp error: mismatch in workspace identifiers!\n");
       exit(1);
    }
    if (*wsp_size == len) {
       DEBUGPRINT("get_double_wsp does nothing, sizes are OK\n");
       return wsp;
    }
    DEBUGPRINT("get_double_wsp re-alloc slot %i: block (%d)->(%d)\n",*wsp_id,*wsp_size,len);
    free((char*)wsp);
    wsp=(double*)(malloc(len*sizeof(double)));
    *wsp_size = len;
    ptr_dwsp[*wsp_id] = wsp;
  } else {
    n_dwsp++;
    if (n_dwsp >= MAX_STATIC_WSP) {
      fprintf(stderr,"overflow error in cm.c: increase MAX_STATIC_WSP\n");
      exit(1);
    }
    DEBUGPRINT("get_double_wsp alloc(%d)-block\n",len);
    ptr_dwsp[n_dwsp]=(double*)(malloc(len*sizeof(double)));
    wsp = ptr_dwsp[n_dwsp];
    *wsp_size = len;
    *wsp_id = n_dwsp;
  }
  return wsp;
}


int n_iwsp=0;
int* ptr_iwsp[MAX_STATIC_WSP+1];

int * get_int_wsp(int len, int * wsp, int * wsp_size, int * wsp_id)
{
  if (wsp != NULL) {
    if (wsp != ptr_iwsp[*wsp_id]) {
       fprintf(stderr,"get_int_wsp error: mismatch in workspace identifiers!\n");
       exit(1);
    }
    if (*wsp_size == len) {
       DEBUGPRINT("get_int_wsp does nothing, sizes are OK\n");
       return wsp;
    }
    DEBUGPRINT("get_int_wsp re-alloc slot %i: block (%d)->(%d)\n",*wsp_id,*wsp_size,len);
    free((char*)wsp);
    wsp=(int*)(malloc(len*sizeof(int)));
    *wsp_size = len;
    ptr_iwsp[*wsp_id] = wsp;
  } else {
    n_iwsp++;
    if (n_iwsp >= MAX_STATIC_WSP) {
      fprintf(stderr,"overflow error in cm.c: increase MAX_STATIC_WSP\n");
      exit(1);
    }
    DEBUGPRINT("get_int_wsp alloc(%d)-block\n",len);
    ptr_iwsp[n_iwsp]=(int*)(malloc(len*sizeof(int)));
    wsp = ptr_iwsp[n_iwsp];
    *wsp_size = len;
    *wsp_id = n_iwsp;
  }
  return wsp;
}


void free_wsp() 
{
  int i;

  if (n_cwsp) {
     DEBUGPRINT("free_wsp dealloc complx\n");
     for (i=1;i<=n_cwsp;i++) {
       DEBUGPRINT(" - deallocating %d\n",i);
       free((char*)(ptr_cwsp[i]));
       ptr_cwsp[i]=NULL;
     }
     n_cwsp=0;
  }
  if (n_dwsp) {
     DEBUGPRINT("free_wsp dealloc double\n");
     for (i=1;i<=n_dwsp;i++) {
       DEBUGPRINT(" - deallocating %d\n",i);
       free((char*)(ptr_dwsp[i]));
       ptr_dwsp[i]=NULL;
     }
     n_dwsp=0;
  }
  if (n_iwsp) {
     DEBUGPRINT("free_wsp dealloc int\n");
     for (i=1;i<=n_iwsp;i++) {
       DEBUGPRINT(" - deallocating %d\n",i);
       free((char*)(ptr_iwsp[i]));
       ptr_iwsp[i]=NULL;
     }
     n_iwsp=0;
  }
}

void free_mv_static()
{
   int i;
   
   if (n_cmv) {
      for (i=1; i<=n_cmv; i++) {
         DEBUGPRINT(" -deallocating cmv %i\n",i);
	 complx_matrix_free(ptr_cmv[i]);
	 ptr_cmv[i]=NULL;
      }
      n_cmv = 0;
   }
   
   if (n_dmv) {
      for (i=1; i<=n_dmv; i++) {
         DEBUGPRINT(" -deallocating dmv %i\n",i);
	 double_matrix_free(ptr_dmv[i]);
	 ptr_dmv[i]=NULL;
      }
      n_dmv = 0;
   }

   if (n_imv) {
      for (i=1; i<=n_imv; i++) {
         DEBUGPRINT(" -deallocating imv %i\n",i);
	 int_matrix_free(ptr_imv[i]);
	 ptr_imv[i]=NULL;
      }
      n_imv = 0;
   }
}

void cmv_real(mv_complx *C, mv_double *D)
{
   int len;
   const int itwo=2;
   const double k=1.0;
   
   len = (C->row)*(C->col);
   memset(D->data,0,len*sizeof(double));
   daxpy_(&len,&k,(double*)C->data,&itwo,D->data,&INTONE);
}

void cmv_imag(mv_complx *C, mv_double *D)
{
   int len;
   const int itwo=2;
   const double k=1.0;
   double * ptr;
   
   len = (C->row)*(C->col);
   ptr = (double*)(C->data) +1;
   memset(D->data,0,len*sizeof(double));
   daxpy_(&len,&k,ptr,&itwo,D->data,&INTONE);
}

void dmv_complx(mv_double *D, mv_complx *C)
{
   int len;
   const int itwo=2;
   const double k=1.0;
   
   len = (D->row)*(D->col);
   memset(C->data,0,len*sizeof(complx));
   daxpy_(&len,&k,D->data,&INTONE,(double*)(C->data),&itwo);
}

int cmv_isreal(mv_complx *C)
{
   int len;
   const int itwo=2;
   double sum;
   double *ptr;
   
   len = (C->row)*(C->col);
   ptr = (double*)(C->data) +1;
   sum = cblas_dasum(len,ptr,itwo);
   if (sum < 1.0e-8) {
      return 1;
   } else {
      return 0;
   }
}

double * get_real_diag(mv_complx *C)
{
  double *res, *sl1;
  int N,i;
  complx *sl2;
  
  N = C->row;
  if (N != C->col) {
     fprintf(stderr,"get_real_diag error: matrix must be rectangular\n");
     exit(1);
  }
  res = sl1 = (double*)malloc(N*sizeof(double));
  sl2 = C->data;
  for (i=0; i<N; i++) {
     *sl1 = sl2->re;
     sl1++;
     sl2 += N+1;
  }
  return res;
}

double * get_real(mv_complx *C)
{
   const int itwo=2;
   int len;
   double * res;
   
   len = (C->row)*(C->col);
   res = (double*)(malloc(len*sizeof(double)));
   dcopy_(&len,(double*)C->data,&itwo,res,&INTONE);
   return res;
}


complx cm_trace(mv_complx * a,mv_complx * b) 
{

   /* equation:       ___
                      \
         trace(A*B)=  /__  A   B
                      i,j   ij   ji
      equality:
          trace(A*B) = trace(B*A)
     
      trace(A,B) is a N^2 process in multiplication relative to
      trace(A*B) which is a N^3 process;
   */
   /* for this case we need to transpose one matrix, let say a */
   int i,j,k,len;
   complx tr;
   double *q, *qq, *qs;
   static complx * wsp=NULL;
   static int wsp_id=0;
   static int wsp_size=0;
   
   len = (a->row)*(a->col);
   wsp = get_complx_wsp(len, wsp, &wsp_size, &wsp_id);

   q = qs = (double*)(a->data);
   qq = (double*)wsp;
   j = (a->row) << 1;
   for (i=0; i<(a->row); i++) {
      for (k=0; k<(a->col); k++) {
         *qq = *q;
	 qq[1] = q[1];
	 qq += 2;
	 q += j;
      }
      q = (qs += 2);
   }
   /* calculate the trace now... */
   cblas_zdotu_sub(len,wsp,INTONE,b->data,INTONE,&tr);
   return tr;
}

complx cm_trace_adjoint(mv_complx * a, mv_complx * b) 
{

   /* equation:       ___
                      \     *   
         trace(A*B)=  /__  A   B
                      i,j   ij   ij
      equality:
          trace(A'*B) = trace(B*A')
     
      trace(A,B) is a N^2 process in multiplication relative to
      trace(A*B) which is a N^3 process;
   */
   /* no need to transpose (done due to column-major storage)
      nor conjugate (done inside zdotcsub_ )                 */
   int len;
   complx tr;
   
   len = (a->row)*(a->col);
   cblas_zdotc_sub(len,a->data,INTONE,b->data,INTONE,&tr);
   return tr;
}

/* logical AND of the real part of two matrices.
   All matrices can be the same pointer. */
   
void cm_and(mv_complx * a_and_b, mv_complx * a,mv_complx * b)
{
  int NN,i, ttrue;

  NN=(a->row)*(a->col);
  if ( (NN != (b->row)*(b->col)) || (NN != (a_and_b->row)*(a_and_b->col)) ) {
    fprintf(stderr,"error: cm_and: matrices must have the same dimensions\n");
    exit(-1);
  }

  for (i=0;i<NN;i++) {
       ttrue = (a->data[i].re != 0.0 && b->data[i].re != 0.0);
       a_and_b->data[i]= Complx( (ttrue != 0.0 ? 1.0 : 0.0) , 0.0);
  }
}


/* logical OR of the real part of two matrices.
   All matrices can be the same pointer. */
   
void cm_or(mv_complx * a_or_b, mv_complx * a,mv_complx * b)
{
  int NN,i, ttrue;

  NN=(a->row)*(a->col);
  if ( NN != (b->row)*(b->col) || NN != (a_or_b->row)*(a_or_b->col) ) {
    fprintf(stderr,"error: cm_or: matrices must have the same dimensions\n");
    exit(-1);
  }

  for (i=0;i<NN;i++) {
       ttrue = (a->data[i].re != 0.0 || b->data[i].re != 0.0);
       a_or_b->data[i]= Complx( (ttrue != 0.0 ? 1.0 : 0.0) , 0.0);
  }
}

void cmv_zero(mv_complx * obj) 
{
   memset(obj->data, 0, (obj->row)*(obj->col)*sizeof(complx));
}

void dmv_zero(mv_double * obj) 
{
   memset(obj->data, 0, (obj->row)*(obj->col)*sizeof(double));
}

/****
 * direct product of two complex matrices
 *     A (*) B = ( a_11*B a_12*B ...  )
 *               ( a_21*B a_22*B ...  )
 *               (  ...               )
 ****/ 
mv_complx * cm_direct(mv_complx * a,mv_complx * b)
{
  int Nar, Nac, Nbr, Nbc, i, j, k;
  mv_complx * res;
  complx * tmp, *da, *db, *dr;
  int s1, s2;

  Nar = a->row; Nac = a->col;
  Nbr = b->row; Nbc = b->col;
  s1 = Nar*Nbr;
  s2 = s1*Nbc;
  
  res = complx_matrix_alloc(s1, Nac*Nbc);
  cmv_zero(res);
  da = a->data;
  db = b->data;
  dr = res->data;
  for (i=0; i<Nbc; i++) {
     tmp = &(db[i*Nbr]); /* this is i-column of b */
     for (j=0; j<Nar; j++) {
        for (k=0; k<Nac; k++) {
	   zaxpy_(&Nbr, &(da[j+k*Nar]), tmp, &INTONE, &(dr[j*Nbr+i*s1+k*s2]), &INTONE);
	}
     }
  } 

  return res;
}

/****
 * sort of direct sum of two complex matrices
 *     A (+) B = ( a_11+B a_12+B ...  )
 *               ( a_21+B a_22+B ...  )
 *               (  ...               )
 *  note: this is different from direct sum and Kronecker sum in Mathematica...
 ****/ 
mv_complx * cm_directadd(mv_complx * a,mv_complx * b)
{
  int Na,Nb,NN,i,j,r,c;
  mv_complx * tmp;
  int s1, ptmp;
  complx aval;
  complx *da, *db, *dr;

  Na = a->row;
  Nb = b->row;
  s1 = Na*Nb;
  NN = s1*s1;
  
  tmp = complx_matrix_alloc(s1, s1);
  cmv_zero(tmp);
  da = a->data;
  db = b->data;
  dr = tmp->data;

  for (r=0; r<Na; r++) {
     for (c=0; c<Na; c++) {
        aval = da[r+c*Na];
        for (i=0; i<Nb; i++) {
	   for (j=0; j<Nb; j++) {
	      ptmp = r*Nb+c*s1*Nb+i+j*s1;
	      dr[ptmp] = Cadd(aval,db[i+j*Nb]);
	   }
	}
     }
  }
  return tmp;
}
    
mv_complx * cmv_dup(mv_complx * obj)
{
  int len;
  mv_complx * tmp;
  
  len = (obj->row)*(obj->col);
  tmp = complx_matrix_alloc(obj->row,obj->col);
  memcpy(tmp->data,obj->data,len*sizeof(complx));
  
  return tmp;
}

void cm_copydiagmv(mv_complx * diagm,mv_complx * v)
{
  /* copy from vector to diagonal of a  matrix */
  int len;
  complx *dm, *dv, *stop;
  
  len=(v->row)*(v->col);
  if ( (len != diagm->row) || (len != diagm->col) ) {
     fprintf(stderr,"cm_copydiagmv error: dimensions mismatch\n"); 
     exit(1);
  }
  cmv_zero(diagm);
  /* for (i=0;i<len;i++) diagm->data[i*(len+1)]=v->data[i]; */

  dm = diagm->data;
  dv = v->data;
  stop = dv + len;
  do {
     *dm = *dv;
     dv++;
     dm += len+1;
  } while ( dv != stop);

}

void cmv_copy(mv_complx * dst, mv_complx * src)
{
   int len,r,c;
   r = src->row;
   c = src->col;
   len = r*c;
   if ( (r != dst->row) || (c != dst->col) ) {
      fprintf(stderr,"cmv_copy error: dimension mismatch\n");
      exit(1);
   }
   memcpy(dst->data, src->data, len*sizeof(complx));
}

void cm_unit(mv_complx * obj)
{
   int N,i;
   
   N = obj->row;
   if (N != obj->col) {
      fprintf(stderr,"cm_unit error: non-square matrix\n");
      exit(1);
   }
   /*cmv_zero(obj);
     for (i=0; i<N; i++) obj->data[i*(N+1)].re = 1.0; */
   /*  CAN BE IMPROVED BY USING LAPACK AUXILIARY FUNCTION */
   zlaset_("A",&N,&N,&CPLX0,&CPLX1,obj->data,&N);
}

/* first column set to real 1.0 */
void cv_d1(mv_complx *obj)
{
   int i;
   complx *sl;
   
   sl = obj->data;
   for (i=0; i<obj->row; i++) {
      *sl = CPLX1;
      sl++;
   }
}
      

/* multiply by complex and overwrite */
void cmv_mulc(mv_complx * obj, complx z)
{
   int len;
   
   len = (obj->row)*(obj->col);
   zscal_(&len, &z, obj->data, &INTONE);
}

/* multiply by double and overwrite */
void cmv_muld(mv_complx * obj, double d)
{
   int len;

   len = (obj->row)*(obj->col);
   zdscal_(&len, &d, obj->data, &INTONE);
}

/* y = y + z*x */
void cmv_multoc(mv_complx * y, mv_complx * x, complx z)
{
   int len;
   
   len = (y->row)*(y->col);
   zaxpy_(&len, &z, x->data, &INTONE, y->data, &INTONE);
}

/* y = y + d*x */
void cmv_multod(mv_complx * y, mv_complx * x, double d)
{
   int len;
   complx tmp;
   
   len = (y->row)*(y->col);
   tmp = Complx(d, 0.0);
   zaxpy_(&len, &tmp, x->data, &INTONE, y->data, &INTONE);
}

/* this simply multiplies two objects element wise */ 
void cmv_mul_elem(mv_complx * res, mv_complx * a, mv_complx * b) 
{
   int len;
   complx *dr, *da, *db, *stop;
   
   len = (res->row)*(res->col);
   if ( (len != (a->row)*(a->col)) || (len != (b->row)*(b->col)) ) {
      fprintf(stderr,"cmv_mul_elem error: mismatch in dimensions\n");
      exit(1);
   }
   dr = res->data;
   da = a->data;
   db = b->data;
   if (dr == da || dr == db) {
      fprintf(stderr,"cmv_mul_elem error: you must use different matrices as arguments\n");
      exit(1);
   }
#ifdef MKL_INTEL
   vzMul(len,da,db,dr);
#else
   stop = dr+len;
   do {
      dr->re = da->re*db->re-da->im*db->im;
      dr->im = da->im*db->re+da->re*db->im;
      dr++;
      da++;
      db++;
   } while(dr != stop);
#endif
}

void cm_mul_fd(mv_complx * res, mv_complx * a, mv_complx * b)
{  /* a is full, b is diagonal */
   int i, r, c;
   complx *zval, *cola, *colres;
   
   r = a->row;
   c = a->col;
   if ( (r != res->row) || (c != res->col) || (r != c) || (r != b->row) ) {
      fprintf(stderr,"cm_mul_fd error: dimension mismatch\n");
      exit(1);
   }
   cmv_zero(res);
   zval = b->data;
   cola = a->data;
   colres = res->data;
   for (i=0; i<c; i++) {
      zaxpy_(&r,zval,cola,&INTONE,colres,&INTONE);
      cola += r;
      colres += r;
      zval++;
   }
}

void cm_mul_df(mv_complx * res, mv_complx * a, mv_complx * b)
{ /* a is diagonal, b is full */
   int i, r, c;
   complx *zval, *rowb, *rowres;
   
   r = b->row;
   c = b->col;
   if ( (r != res->row) || (c != res->col) || (r != c) || (r != a->row) ) {
      fprintf(stderr,"cm_mul_fd error: dimension mismatch\n");
      exit(1);
   }
   cmv_zero(res);
   zval = a->data;
   rowb = b->data;
   rowres = res->data;
   for (i=0; i<r; i++) {
      zaxpy_(&c,zval,rowb,&r,rowres,&r);
      rowb++;
      rowres++;
      zval++;
   }
}


#include <errno.h>
void cm_mul(mv_complx * res, mv_complx * a, mv_complx * b)
{
   complx *dr, *da, *db;
   int M, N, K;
   
   dr = res->data;
   da = a->data;
   db = b->data;
   if (dr == da || dr == db) {
      fprintf(stderr,"cm_mul error: you must use different matrices as arguments\n");
      exit(1);
   }
   M = a->row;
   N = b->col;
   K = a->col;
   if ( (M != res->row) || (N != res->col) || (K != b->row) ) {
      fprintf(stderr,"cm_mul error: dimensions mismatch\n");
      exit(1);
   }
   zgemm_("N","N",&M,&N,&K,&CPLX1,da,&M,db,&K,&CPLX0,dr,&M);
   errno = 0;
}

/* matrices must be square and all of equal dim, not checked here */
void cm_commutator(mv_complx * res, mv_complx * a, mv_complx * b)
{
   complx *dr, *da, *db;
   int N, len;
   static complx * wsp=NULL;
   static int wsp_id=0;
   static int wsp_size=0;
   
   dr = res->data;
   da = a->data;
   db = b->data;
   N = res->row;
   len = N*N;
   wsp = get_complx_wsp(len, wsp, &wsp_size, &wsp_id);
   zgemm_("N","N",&N,&N,&N,&CPLX1,da,&N,db,&N,&CPLX0,wsp,&N);
   zgemm_("N","N",&N,&N,&N,&CPLXm1,db,&N,da,&N,&CPLX1,wsp,&N);
   memcpy(dr,wsp,len*sizeof(complx));
}

complx cv_dotmul(mv_complx * a, mv_complx * b)
{
   complx res;
   int len;
   
   len = a->row;
   if ( len != b->row ) {
      fprintf(stderr,"cv_dotmul error: vector dimension mismatch\n");
      exit(1);
   }
   cblas_zdotu_sub(len,a->data,INTONE,b->data,INTONE,&res);
   return res;
}

void cm_mulmv(mv_complx * y, mv_complx * A, mv_complx * x)
{
  /*
    Only for quadratic matrices
    col_vector y = matrix A * col_vector x 
  */
   complx *dy, *dx;
   int r, c;
   
   dy = y->data;
   dx = x->data;
   if ( dx == dy ) {
      fprintf(stderr,"cm_mulmv error: you must use different vectors as arguments\n");
      exit(1);
   }
   r = A->row;
   c = A->col;
   if ( (c != x->row) || (c != y->row) ) {
      fprintf(stderr,"cm_mulmv error: dimension mismatch\n");
      exit(1);
   }
   zgemv_("N",&r,&c,&CPLX1,A->data,&r,dx,&INTONE,&CPLX0,dy,&INTONE);
} 

void cm_mulvm(mv_complx * y, mv_complx * x, mv_complx * A)
{
  /*
    Only for quadratic matrices
    row_vector y = row_vector x * matrix A 
  */
   complx *dy, *dx, *dA;
   int r, c;
   
   dy = y->data;
   dx = x->data;
   dA = A->data;
   if ( dx == dy ) {
      fprintf(stderr,"cm_mulvm error: you must use different vectors as arguments\n");
      exit(1);
   }
   r = A->row;
   c = A->col;
   if ( (r != x->row) || (c != y->row) ) {
      fprintf(stderr,"cm_mulmv error: dimension mismatch\n");
      exit(1);
   }
   zgemv_("T",&r,&c,&CPLX1,dA,&r,dx,&INTONE,&CPLX0,dy,&INTONE);
} 

/* I define this but try not to use it! Do adjoint in multiplication right way... */
void cm_adjoint(mv_complx *aplus, mv_complx * a)
{
   int i,j,k,r,c;
   double *q, *qq, *qs;
   complx *da, *dap;
   
   da = a->data;
   dap = aplus->data;
   if ( da == dap ) {
      fprintf(stderr,"cm_adjoint error: you must use diggerent matrices\n");
      exit(1);
   }
   r = a->row;
   c = a->col;
   if ( (r != aplus->col) || (c != aplus->row) ) {
      fprintf(stderr,"cm_adjoint error: dimension mismatch\n");
      exit(1);
   }

   q = qs = (double*)da;
   qq = (double*)dap;
   j = r << 1;
   for (i=0; i<r; i++) {
      for (k=0; k<c; k++) {
         *qq = *q;
	 qq[1] = -q[1];
	 qq += 2;
	 q += j;
      }
      q = (qs += 2);
   }
}

int cm_ishermit(mv_complx *a)
{
   int i, j, N;
   complx *ptr, v1, v2;
   
   N = a->row;
   if (N != a->col) return 0;  
   ptr = a->data;
   for (i=0; i<N; i++) {
      if ( fabs(ptr->im) > 1e-8 ) return 0;
      ptr += (N+1);
   }
   for (i=0; i<N; i++) {
      for (j=i+1; j<N; j++) {
         v1 = a->data[i+N*j];
	 v2 = a->data[j+N*i];
         if (fabs(v1.re - v2.re) > 1e-8) return 0;
         if (fabs(v1.im + v2.im) > 1e-8) return 0;
      }
   }
   return 1;
}

void cv_conj(mv_complx *conja, mv_complx * a)
{
  int N;
  const double dmone=-1.0;
  double *ptr;
  const int itwo=2;
  
  N = a->row;
  if (N != conja->row) {
     fprintf(stderr,"cv_conj error: dimension mismatch\n");
     exit(1);
  }
  memcpy(conja->data,a->data,N*sizeof(complx));
  ptr = (double*)(conja->data) + 1;
  dscal_(&N,&dmone,ptr,&itwo);
}

void cv_conj_in(mv_complx * a)
{
  int N;
  const double dmone=-1.0;
  double *ptr;
  const int itwo=2;
  
  N = a->row;
  ptr = (double*)(a->data) + 1;
  dscal_(&N,&dmone,ptr,&itwo);
}

void cm_conj(mv_complx *conja, mv_complx * a)
{
  int r,c,len;
  const double dmone=-1.0;
  double *ptr;
  const int itwo=2;
  
  r = a->row;
  c = a->col;
  if ( (r != conja->row) || (c != conja->col) ) {
     fprintf(stderr,"cm_conj error: dimension mismatch\n");
     exit(1);
  }
  len = r*c;
  memcpy(conja->data,a->data,len*sizeof(complx));
  ptr = (double*)(conja->data) + 1;
  dscal_(&len,&dmone,ptr,&itwo);
}

double cmv_sumnorm1(mv_complx *tmp)
{
   double res;
   int len;
   
   len = (tmp->row)*(tmp->col);
   res = cblas_dzasum(len,tmp->data,INTONE);
   return res;
}


#define FMT " %9.6g"

void cm_print(mv_complx * m,char* name)
{
  int rows,cols,i,j;
  complx * dm;
  
  rows=m->row;
  cols=m->col;
  dm = m->data;

  printf("%s : matrix(%i,%i)\nreal part:\n",name,rows,cols);
  for (i=0;i<rows;i++) {
    printf(" ");
    for (j=0;j<cols;j++) {
      printf(FMT,dm[i+j*rows].re);
    }
    printf("\n");
  }
  printf("imag part:\n");
  for (i=0;i<rows;i++) {
    printf(" ");
    for (j=0;j<cols;j++) {
      printf(FMT,dm[i+j*rows].im);
    }
    printf("\n");
  }
}

void dm_print(mv_double * m,char* name)
{
  int rows,cols,i,j;
  double * dm;
  
  rows=m->row;
  cols=m->col;
  dm = m->data;

  printf("%s : matrix(%i,%i)\n",name,rows,cols);
  for (i=0;i<rows;i++) {
    printf(" ");
    for (j=0;j<cols;j++) {
      printf(FMT,dm[i+j*rows]);
    }
    printf("\n");
  }
}

void cv_print(mv_complx * v,char* name)
{
  int len,i;
  complx *dv;
  
  len=v->row;
  dv = v->data;
  printf("%s : vector(%i)\nreal part:\n ",name,len);
  for (i=0;i<len;i++) printf(FMT,dv[i].re);
  printf("\nimag part:\n ");
  for (i=0;i<len;i++) printf(FMT,dv[i].im);
  printf("\n");
}

void dv_print(mv_double * v,char* name)
{
  int len,i;
  double *dv;
  
  len=v->row;
  dv = v->data;
  printf("%s : vector(%i)\n ",name,len);
  for (i=0;i<len;i++) printf(FMT,dv[i]);
  printf("\n");
}

void simtrans_diagprop1(mv_complx * A, mv_complx * U)
{
  /* this uses two for loops, direct approach without libraries */
  /*  U is stored as a vector holding the diagonal */
  int i,j,N;
  complx *dA, *dU, *dU2;
  
  N = U->row;
  dA = A->data;
  dU = U->data;
  for (i=0; i<N; i++) { /* columns */
     double uiire = dU->re;
     double uiiim = dU->im;
     dU2 = U->data;
     for (j=0; j<N; j++) { /* rows */
        double re = dA->re*uiire + dA->im*uiiim;
	double im = dA->im*uiire - dA->re*uiiim;
        dA->re = re*dU2->re - im*dU2->im;
	dA->im = im*dU2->re + re*dU2->im;
	dU2++;
	dA++;
     }
     dU++;
  }
}

void simtrans_diagprop2(mv_complx * A, mv_complx * U)
{
  /* this creates a mask matrix and element-wise multiplication */
  /*  U is stored as a vector holding the diagonal */

  /* version from Nov 2009, goes througn whole matrix */
   static complx * wsp=NULL;
   static int wsp_id=0;
   static int wsp_size=0;
   complx * dU, *dA, *sl, *stop;
   int N, len;
   
   dU = U->data;
   dA = A->data;
   N = U->row;
   len = N*N;
   wsp = get_complx_wsp(len, wsp, &wsp_size, &wsp_id);
   memset(wsp,0,len*sizeof(complx));
   zgerc_(&N,&N,&CPLX1,dU,&INTONE,dU,&INTONE,wsp,&N);
   sl = wsp;
   stop = sl+len;
   do {
      double re = dA->re;
      double im = dA->im;
      dA->re = re*sl->re - im*sl->im;
      dA->im = im*sl->re + re*sl->im;
      dA++;
      sl++;
   } while(sl != stop); 
}

void simtrans_diagprop3(mv_complx * A, mv_complx * U)
{
  
  /* version Dec 2009, goes through lower triangle */
  int i,j,N;
  complx *dA, *dA2, *dU, *dU2;
  
  N = U->row;
  dU = U->data;
  for (i=0; i<N; i++) { /* columns */
     double uiire = dU->re;
     double uiiim = dU->im;
     dU2 = dU;
     dA = dA2 = A->data + (N+1)*i;
     for (j=i+1; j<N; j++) { /* rows, lower triangle */
        dA++;
	dA2 += N;
	dU2++;
	/* if ( (fabs(dA->re) < 1e-8) && (fabs(dA->im) < 1e-8) )
	   continue;
	*/
        double u2iire = dU2->re;
        double u2iiim = dU2->im;
        double re = dA->re*uiire + dA->im*uiiim;
	double im = dA->im*uiire - dA->re*uiiim;
        dA->re = re*u2iire - im*u2iiim;
	dA->im = im*u2iire + re*u2iiim;

        re = dA2->re*uiire - dA2->im*uiiim;
	im = dA2->im*uiire + dA2->re*uiiim;
        dA2->re = re*u2iire + im*u2iiim;
	dA2->im = im*u2iire - re*u2iiim;
     }
     dU++;
  }
  
}

void simtrans_adj_diagprop3(mv_complx * A, mv_complx * U)
{
  /* version Dec 2009, goes through lower triangle */
  int i,j,N;
  complx *dA, *dA2, *dU, *dU2;
  
  N = U->row;
  dU = U->data;
  for (i=0; i<N; i++) { /* columns */
     double uiire = dU->re;
     double uiiim = dU->im;
     dU2 = dU;
     dA = dA2 = A->data + (N+1)*i;
     for (j=i+1; j<N; j++) { /* rows, lower triangle */
        dA++;
	dA2 += N;
	dU2++;
	/* if ( (fabs(dA->re) < 1e-8) && (fabs(dA->im) < 1e-8) )
	   continue;
	*/
        double u2iire = dU2->re;
        double u2iiim = dU2->im;
        double re = dA->re*uiire - dA->im*uiiim;
	double im = dA->im*uiire + dA->re*uiiim;
        dA->re = re*u2iire + im*u2iiim;
	dA->im = im*u2iire - re*u2iiim;

        re = dA2->re*uiire + dA2->im*uiiim;
	im = dA2->im*uiire - dA2->re*uiiim;
        dA2->re = re*u2iire - im*u2iiim;
	dA2->im = im*u2iire + re*u2iiim;
     }
     dU++;
  }
  
}


void simtrans_zrot(mv_complx * U, double * P)
{
   /*    exp{-iP} U exp{iP}, where P is diagonal                     */
   /*    P is passed as a real vector and defines the z-rotation     */
   /*    ( do exp(-iP), mask via zgerc, elem-wise mult mask*U        */
  
  int i,j,N;
  complx *dU;
  double * dP, *dP2, ph, c, s, rr, ii;
  
  N = U->row;
  dP = P;
  dU = U->data;
  for (i=0; i<N; i++) { /* columns */
     dP2 = P;
     for (j=0; j<N; j++) { /* rows */
        ph = *dP - *dP2;
	c = cos(ph);
	s = sin(ph);
        rr = dU->re;
	ii = dU->im;
	dU->re = rr*c - ii*s;
	dU->im = ii*c + rr*s;
	dP2++;
	dU++;
     }
     dP++;
  }
}

void simtrans_zrot2(mv_complx * U, double * P)
{
   /*    exp{-iP} U exp{iP}, where P is diagonal                     */
   /*    P is passed as a real vector and defines the z-rotation     */
  int i,j,N;
  complx *dU, *dU2;
  double *dP, *dP2, ph, c, s, rr, ii;
  
  N = U->row;
  dP = P;
  for (i=0; i<N; i++) { /* columns */
     dU = dU2 = U->data + (N+1)*i;
     dP2 = dP;
     for (j=i+1; j<N; j++) { /* rows, lower triangle */
        dU++;
	dU2 += N;
	dP2++;
        /* if ( (fabs(dU->re) < 1e-8) && (fabs(dU->im) < 1e-8) )
	   continue;
	*/
       	ph = *dP - *dP2;
	c = cos(ph);
	s = sin(ph);
        rr = dU->re;
	ii = dU->im;
	dU->re = rr*c - ii*s;
	dU->im = ii*c + rr*s;
        rr = dU2->re;
	ii = dU2->im;
	dU2->re = rr*c + ii*s;
	dU2->im = ii*c - rr*s;
     }
     dP++;
  }
}

void ham_zrot_real(mv_complx *Ham, double * P)
{
   /*    Ham = exp{-iP} real(Ham) exp{iP}, where P is diagonal       */
   /*    P is passed as a real vector and defines the z-rotation     */
  int i,j,N;
  complx *dH, *dH2;
  double *dP, *dP2, ph, c, s, rr, ii;
  
  N = Ham->row;
  dP = P;
  for (i=0; i<N; i++) { /* columns */
     dH = dH2 = Ham->data + (N+1)*i;
     dP2 = dP;
     for (j=i+1; j<N; j++) { /* rows, lower triangle */
        dH++;
	dH2 += N;
	dP2++;
        /* if ( (fabs(dU->re) < 1e-8) && (fabs(dU->im) < 1e-8) )
	   continue;
	*/
       	ph = *dP - *dP2;
	c = cos(ph);
	s = sin(ph);
        rr = dH->re;
	dH->re = rr*c;
	dH->im = rr*s;
        rr = dH2->re;
	dH2->re = rr*c;
	dH2->im = - rr*s;
     }
     dP++;
  }
}

void simtrans(mv_complx * S, mv_complx * U)
{
   /* S = U S U+  , general complex matrices */
   static complx * wsp=NULL;
   static int wsp_id=0;
   static int wsp_size=0;
   int N, len;
   complx * wsp2;
   
   N = U->row;
   len = N*N;
   wsp = get_complx_wsp(2*len, wsp, &wsp_size, &wsp_id);
   wsp2 = wsp+len;
   zgemm_("N","C",&N,&N,&N,&CPLX1,S->data,&N,U->data,&N,&CPLX0,wsp,&N);
   zgemm_("N","N",&N,&N,&N,&CPLX1,U->data,&N,wsp,&N,&CPLX0,wsp2,&N);
   memcpy(S->data,wsp2,len*sizeof(complx));
}

void simtransh(mv_complx * S, mv_complx * U)
{
   /* S = U S U+  , S is hermitian */
   static complx * wsp=NULL;
   static int wsp_id=0;
   static int wsp_size=0;
   int N, len;
   complx * wsp2;
   
   N = U->row;
   len = N*N;
   wsp = get_complx_wsp(2*len, wsp, &wsp_size, &wsp_id);
   wsp2 = wsp+len;
   zhemm_("R","U",&N,&N,&CPLX1,S->data,&N,U->data,&N,&CPLX0,wsp,&N);
   zgemm_("N","C",&N,&N,&N,&CPLX1,wsp,&N,U->data,&N,&CPLX0,wsp2,&N);
   memcpy(S->data,wsp2,len*sizeof(complx));
}

void simtrans_adj(mv_complx * S, mv_complx * U)
{
   /* S = U+ S U  , general complex matrices */
   static complx * wsp=NULL;
   static int wsp_id=0;
   static int wsp_size=0;
   int N, len;
   complx * wsp2;
   
   N = U->row;
   len = N*N;
   wsp = get_complx_wsp(2*len, wsp, &wsp_size, &wsp_id);
   wsp2 = wsp+len;
   zgemm_("N","N",&N,&N,&N,&CPLX1,S->data,&N,U->data,&N,&CPLX0,wsp,&N);
   zgemm_("C","N",&N,&N,&N,&CPLX1,U->data,&N,wsp,&N,&CPLX0,wsp2,&N);
   memcpy(S->data,wsp2,len*sizeof(complx));
}

void simtransh_adj(mv_complx * S, mv_complx * U)
{
   /* S = U+ S U  , S is hermitian */
   static complx * wsp=NULL;
   static int wsp_id=0;
   static int wsp_size=0;
   int N, len;
   complx * wsp2;
   
   N = U->row;
   len = N*N;
   wsp = get_complx_wsp(2*len, wsp, &wsp_size, &wsp_id);
   wsp2 = wsp+len;
   zhemm_("L","U",&N,&N,&CPLX1,S->data,&N,U->data,&N,&CPLX0,wsp,&N);
   zgemm_("C","N",&N,&N,&N,&CPLX1,wsp,&N,U->data,&N,&CPLX0,wsp2,&N);
   memcpy(S->data,wsp2,len*sizeof(complx));
}


void prop_realdiag(mv_complx * res, mv_double * dv, double dt)
{
   /*   both arguments are vectors!!!                  */
   
   int N;
   complx *dres;
   double *ddv, *stop, arg;
   
   N = dv->row;
   if (N != res->row) {
      fprintf(stderr,"prop_realdiag error: dimension mismatch\n");
      exit(1);
   }
   dres = res->data;
   ddv = dv->data;
   stop = ddv + N;
   do {
      arg = (*ddv)*dt;
      dres->re = cos(arg);
      dres->im = -sin(arg);
      ddv++;
      dres++;
   } while(ddv != stop);
}


void prop_real(mv_complx * prop, mv_double * ham, double dt)
{
   /* real matrix exponential using diagonalization of symmetric matrix

        exp( -i * ham * dt )
    */
   static double * wsp1=NULL;
   static int wsp1_id=0;
   static int wsp1_size=0;
   static double * wsp2=NULL;
   static int wsp2_id=0;
   static int wsp2_size=0;
   int N, len, lwsp2, info, i, j;
   double *eigs, *dT, *extrawsp, dum, cs, ss;

   static complx * wsp3=NULL;
   static int wsp3_id=0;
   static int wsp3_size=0;
   complx *dV;
   
   N = ham->row;
   len = N*N;
   wsp1 = get_double_wsp(3*len+N, wsp1, &wsp1_size, &wsp1_id);
   memcpy(wsp1,ham->data,len*sizeof(double));
   dscal_(&len,&dt,wsp1,&INTONE);
   eigs = wsp1+len;
   extrawsp = eigs+N;
   lwsp2 = -1;
   dsyev_("V","U",&N,wsp1,&N, eigs, &dum, &lwsp2, &info);
   lwsp2 = (int)dum;
   wsp2 = get_double_wsp(lwsp2, wsp2, &wsp2_size, &wsp2_id);
   dsyev_("V","U",&N,wsp1,&N, eigs, wsp2, &lwsp2, &info);
   if ( info != 0) {
      fprintf(stderr,"prop_real error: diagonalization failed\n");
      exit(1);
   }
   wsp3 = get_complx_wsp(len, wsp3, &wsp3_size, &wsp3_id);
   dT = wsp1;
   for (j=0; j<N; j++) {
      dV = wsp3 + j;
      cs = cos(*eigs);
      ss = -sin(*eigs);
      for (i=0;  i<N; i++) {
         dum = *dT;
         dV->re = dum*cs;
	 dV->im = dum*ss;
	 dT++;
	 dV += N;
      }
      eigs++;
   }
   /* result of this(^) should be wsp3=complx(exp(-i*eigs)*Tansf.Mat.+) */
   zlarcm_(&N,&N,wsp1,&N,wsp3,&N,prop->data,&N,extrawsp);
}

/* auxiliary for Pade */
void add_double_diag(double d, double *A, int N)
{
   double *s1;
   int i;
   s1 = A;
   
   for (i=0; i<N; i++) {
      *s1 += d;
      s1 += (N+1);
   }
}

/* auxiliary for Pade */
void cplx_add_double_diag(double d, complx *A, int N)
{
   complx *s1;
   int i;
   s1 = A;
   
   for (i=0; i<N; i++) {
      s1->re += d;
      s1 += (N+1);
   }
}

void prop_pade_real(mv_complx * prop, mv_double * ham, double dt)
{
   /* propagator by using Pade approximation with scaling & squaring 
        exp( -i * ham * dt )
      keeping real matrices as much as possible
    */
   double c[7], infnorm, scale, scale2;
   const double dtwo=2.0;
   const double done=1.0;
   const double dzero=0.0;
   static double * wsp1=NULL;
   static int wsp1_id=0;
   static int wsp1_size=0;
   int N, len, ns, info, i;
   double *pom,*yy,*hdt,*mb,*ma, *dumwsp;
   complx *cmx1, *cmx2, *sl1, *sl2, *pd;
   
   c[0] = 1.0;
   c[1] = 0.5;
   c[2] = 5.0/44.0;
   c[3] = 1.0/66.0;
   c[4] = 1.0/792.0;
   c[5] = 1.0/15840.0;
   c[6] = 1.0/665280.0;

   N = ham->row;
   len = N*N;
   pd = prop->data;
   wsp1 = get_double_wsp(7*len, wsp1, &wsp1_size, &wsp1_id);
   pom = wsp1;
   yy = pom+len;
   hdt = yy+len;
   mb = hdt+len;
   ma = mb+len;
   dumwsp = ma+len;
   memcpy(hdt,ham->data,len*sizeof(double));
   dscal_(&len,&dt,hdt,&INTONE);
   
   /* scaling */
   infnorm = -10.0;
   do {
      scale2 = cblas_dasum(N,hdt,INTONE);
      if (scale2>infnorm) infnorm = scale2;
      hdt += N;
   } while (hdt != mb);
   hdt = yy+len;
   if (infnorm <= 0) {
     fprintf(stderr,"prop_pade_real error: Hamiltonian has non-positive norm\n");
     exit(1);
   }
   ns = (int)(floor(log(infnorm)/log(2.0)))+2;
   if (ns<0) ns = 2;
   scale = -pow(2.0,-ns);
   scale2 = -scale*scale;
   dsymm_("L","U",&N,&N,&scale2,hdt,&N,hdt,&N,&dzero,yy,&N);  /* yy ~= ham*ham */
   
   memset(pom,0,len*sizeof(double));
   add_double_diag(c[6],pom,N);
   dgemm_("N","N",&N,&N,&N,&done,pom,&N,yy,&N,&dzero,mb,&N);
   add_double_diag(c[4],mb,N);
   dgemm_("N","N",&N,&N,&N,&done,mb,&N,yy,&N,&dzero,pom,&N);
   add_double_diag(c[2],pom,N);
   dgemm_("N","N",&N,&N,&N,&done,pom,&N,yy,&N,&dzero,mb,&N);
   add_double_diag(c[0],mb,N);
   
   memset(ma,0,len*sizeof(double));
   add_double_diag(c[5],ma,N);
   dgemm_("N","N",&N,&N,&N,&done,ma,&N,yy,&N,&dzero,pom,&N);
   add_double_diag(c[3],pom,N);
   dgemm_("N","N",&N,&N,&N,&done,pom,&N,yy,&N,&dzero,ma,&N);
   add_double_diag(c[1],ma,N);
   
   cmx1 = sl1 = (complx*)pom;
   cmx2 = sl2 = (complx*)dumwsp;
   pom = hdt;
   memset(cmx1,0,len*sizeof(complx));
   memset(cmx2,0,len*sizeof(complx));
   do {
      sl1->im = scale*(*hdt);
      sl2->re = *mb;
      sl1++;
      sl2++;
      hdt++;
      mb++;
   } while ( mb != ma);
   zlarcm_(&N,&N,ma,&N,cmx1,&N,pd,&N,pom);
   zaxpy_(&len,&CPLXm1,pd,&INTONE,cmx2,&INTONE);

   long int pvec[N+1];   
   zgesv_(&N,&N,cmx2,&N,pvec,pd,&N,&info);
   if (info != 0) {
      fprintf(stderr,"prop_pade_real error: zgesv failed with info=%d\n",info);
      exit(1);
   }
   zdscal_(&len,&dtwo,pd,&INTONE);
   cplx_add_double_diag(1.0,pd,N); /* hotovo */
   
   /* squaring */
   for (i=0; i<ns; i++) {
      if (i%2 == 0) {
         zgemm_("N","N",&N,&N,&N,&CPLX1,pd,&N,pd,&N,&CPLX0,cmx1,&N);
      } else {
         zgemm_("N","N",&N,&N,&N,&CPLX1,cmx1,&N,cmx1,&N,&CPLX0,pd,&N);
      }
   }
   /* printf("prop_pade_real info: # of squaring = %d\n",ns); */
   if (ns%2 != 0) memcpy(pd,cmx1,len*sizeof(complx));
}

void prop_pade_complx(mv_complx * prop, mv_complx * ham, double dt)
{
   /* propagator by using Pade approximation with scaling & squaring 
        exp( -i * ham * dt )
      everything is in complex matrices
    */
   double c[7], infnorm, rowsum, zabs;
   const double dtwo=2.0;
   int N, len, ns, info, i, j;
   static complx * wsp1=NULL;
   static int wsp1_id=0;
   static int wsp1_size=0;
   complx *mx1, *mx2, *mx3, *mx4, *pd;
   complx scale, scale2;
   
   c[0] = 1.0;
   c[1] = 0.5;
   c[2] = 5.0/44.0;
   c[3] = 1.0/66.0;
   c[4] = 1.0/792.0;
   c[5] = 1.0/15840.0;
   c[6] = 1.0/665280.0;

   N = ham->row;
   len = N*N;
   wsp1 = get_complx_wsp(4*len, wsp1, &wsp1_size, &wsp1_id);
   mx1 = wsp1;
   mx2 = mx1+len;
   mx3 = mx2+len;
   mx4 = mx3+len;
   pd = prop->data;
   memcpy(mx4,ham->data,len*sizeof(complx));
   zdscal_(&len,&dt,mx4,&INTONE);

   infnorm = -10.0;
   double *dsl, *dslst;
   dsl = dslst = (double*)mx1;
   memset(dsl,0,N*sizeof(double));
   for (i=0; i<N; i++) {
      mx4 = mx3+len+i*(N+1); dsl = dslst+i;
      rowsum = fabs(mx4->re) + (*dsl);
      for (j=i+1; j<N; j++) {
         mx4++; dsl++;
         zabs = cblas_dznrm2(INTONE,mx4,INTONE);
	 rowsum += zabs;
	 *dsl += zabs;
      }
      if (rowsum > infnorm) infnorm=rowsum;
   }
   mx4 = mx3+len;
   if (infnorm <= 0) {
     fprintf(stderr,"prop_pade_complx error: Hamiltonian has non-positive norm\n");
     exit(1);
   }
   ns = (int)(floor(log(infnorm)/log(2.0)))+2;
   if (ns<0) ns = 2;
   scale.im = -pow(2.0,-ns);
   scale2.re = -scale.im*scale.im;
   scale.re = scale2.im = 0.0;
   zhemm_("L","L",&N,&N,&scale2,mx4,&N,mx4,&N,&CPLX0,pd,&N);  /* pd ~= ham*ham */

   memset(mx2,0,len*sizeof(complx));
   cplx_add_double_diag(c[6],mx2,N);
   zgemm_("N","N",&N,&N,&N,&CPLX1,mx2,&N,pd,&N,&CPLX0,mx1,&N);
   cplx_add_double_diag(c[4],mx1,N);
   zgemm_("N","N",&N,&N,&N,&CPLX1,mx1,&N,pd,&N,&CPLX0,mx2,&N);
   cplx_add_double_diag(c[2],mx2,N);
   zgemm_("N","N",&N,&N,&N,&CPLX1,mx2,&N,pd,&N,&CPLX0,mx1,&N);
   cplx_add_double_diag(c[0],mx1,N);

   memset(mx2,0,len*sizeof(complx));
   cplx_add_double_diag(c[5],mx2,N);
   zgemm_("N","N",&N,&N,&N,&CPLX1,mx2,&N,pd,&N,&CPLX0,mx3,&N);
   cplx_add_double_diag(c[3],mx3,N);
   zgemm_("N","N",&N,&N,&N,&CPLX1,mx3,&N,pd,&N,&CPLX0,mx2,&N);
   cplx_add_double_diag(c[1],mx2,N);
   zgemm_("N","N",&N,&N,&N,&scale,mx2,&N,mx4,&N,&CPLX0,pd,&N);
   
   zaxpy_(&len,&CPLXm1,pd,&INTONE,mx1,&INTONE);
   long int pvec[N+1];   
   zgesv_(&N,&N,mx1,&N,pvec,pd,&N,&info);
   if (info != 0) {
      fprintf(stderr,"prop_pade_complx error: zgesv failed with info=%d\n",info);
      exit(1);
   }
   zdscal_(&len,&dtwo,pd,&INTONE);
   cplx_add_double_diag(1.0,pd,N);
   
   /* squaring */
   for (i=0; i<ns; i++) {
      if (i%2 == 0) {
         zgemm_("N","N",&N,&N,&N,&CPLX1,pd,&N,pd,&N,&CPLX0,mx3,&N);
      } else {
         zgemm_("N","N",&N,&N,&N,&CPLX1,mx3,&N,mx3,&N,&CPLX0,pd,&N);
      }
   }
   if (ns%2 != 0) memcpy(pd,mx3,len*sizeof(complx));

}


/* auxiliary for Chebyshev */
void cplx_add_complx_diag(complx d, complx *A, int N)
{
   complx *s1;
   int i;
   s1 = A;
   
   for (i=0; i<N; i++) {
      s1->re += d.re;
      s1->im += d.im;
      s1 += (N+1);
   }
}

void prop_cheb(mv_complx * prop, mv_complx * ham, double dt)
{
   /*    exp( -i * ham * dt )                                    */
   /* propagator using Chebyshev approx. and complex Hamiltonian */
   /*  there is probably no advantage in using real hamiltonian..*/
   complx alpha[14],theta[14], alpha0, dum;
   int i, N, len, info;
   static complx * wsp1=NULL;
   static int wsp1_id=0;
   static int wsp1_size=0;
   complx *mx1, *mx2, *mx3, *pd, *aa;

   alpha0.re  =  0.183216998528140087e-11; alpha0.im = 0.0;
   alpha[0].re = alpha[7].re = 0.557503973136501826e+02*0.5;
   alpha[1].re = alpha[8].re =  -0.938666838877006739e+02*0.5;
   alpha[2].re = alpha[9].re = 0.469965415550370835e+02*0.5;
   alpha[3].re = alpha[10].re =  -0.961424200626061065e+01*0.5;
   alpha[4].re = alpha[11].re = 0.752722063978321642e+00*0.5;
   alpha[5].re = alpha[12].re =  -0.188781253158648576e-01*0.5;
   alpha[6].re = alpha[13].re = 0.143086431411801849e-03*0.5;
   alpha[0].im = -0.204295038779771857e+03*0.5;
   alpha[1].im = 0.912874896775456363e+02*0.5;
   alpha[2].im = -0.116167609985818103e+02*0.5;
   alpha[3].im = -0.264195613880262669e+01*0.5;
   alpha[4].im = 0.670367365566377770e+00*0.5;
   alpha[5].im = -0.343696176445802414e-01*0.5;
   alpha[6].im = 0.287221133228814096e-03*0.5;
   alpha[7].im = -alpha[0].im;
   alpha[8].im = -alpha[1].im;
   alpha[9].im = -alpha[2].im;
   alpha[10].im = -alpha[3].im;
   alpha[11].im = -alpha[4].im;
   alpha[12].im = -alpha[5].im;
   alpha[13].im = -alpha[6].im;

   /* Warning!!! All theta's have oposite sign... */
   theta[0].re = theta[7].re = 0.562314417475317895e+01;
   theta[0].im =  -0.119406921611247440e+01;
   theta[1].re = theta[8].re = 0.508934679728216110e+01;
   theta[1].im =  -0.358882439228376881e+01;
   theta[2].re = theta[9].re = 0.399337136365302569e+01;
   theta[2].im =  -0.600483209099604664e+01;
   theta[3].re = theta[10].re = 0.226978543095856366e+01;
   theta[3].im =  -0.846173881758693369e+01;
   theta[4].re =  theta[11].re = -0.208756929753827868e+00;
   theta[4].im =  -0.109912615662209418e+02;
   theta[5].re =  theta[12].re = -0.370327340957595652e+01;
   theta[5].im =  -0.136563731924991884e+02;
   theta[6].re =  theta[13].re = -0.889777151877331107e+01;
   theta[6].im =  -0.166309842834712071e+02; 
   theta[7].im =  0.119406921611247440e+01;
   theta[8].im =  0.358882439228376881e+01;
   theta[9].im =  0.600483209099604664e+01;
   theta[10].im =  0.846173881758693369e+01;
   theta[11].im =  0.109912615662209418e+02;
   theta[12].im =  0.136563731924991884e+02;
   theta[13].im =  0.166309842834712071e+02; 

   N = ham->row;
   len = N*N;
   wsp1 = get_complx_wsp(3*len, wsp1, &wsp1_size, &wsp1_id);
   mx1 = wsp1;
   mx2 = mx1+len;
   mx3 = mx2+len;
   pd = prop->data;
   dum.re = 0.0; dum.im = dt;
   /* this can be improved using intel MKL out-place-scaling */
   memcpy(mx3,ham->data,len*sizeof(complx));
   zscal_(&len,&dum,mx3,&INTONE);
   
   zlaset_("A",&N,&N,&CPLX0,&alpha0,pd,&N);
   
   aa = alpha;
   long int pvec[N+1];
   for (i=0; i<14; i++) {
      memcpy(mx1,mx3,len*sizeof(complx));
      cplx_add_complx_diag(theta[i],mx1,N);
      zlaset_("A",&N,&N,&CPLX0,aa,mx2,&N);
      zgesv_(&N,&N,mx1,&N,pvec,mx2,&N,&info);
      if (info != 0) { 
         fprintf(stderr,"prop_cheb error in step %d : zgesv failed with code %d\n",i,info);
	 exit(1);
      }
      zaxpy_(&len,&CPLX1,mx2,&INTONE,pd,&INTONE);
      aa++;
   }
}



/* this does A = A*B */
void cm_multo(mv_complx * ab, mv_complx * b)
{
   static complx * wsp1=NULL;
   static int wsp1_id=0;
   static int wsp1_size=0;
   int N,len;
   
   N = b->row;
   len = N*N;
   wsp1 = get_complx_wsp(len, wsp1, &wsp1_size, &wsp1_id);
   
   zgemm_("N","N",&N,&N,&N,&CPLX1,ab->data,&N,b->data,&N,&CPLX0,wsp1,&N);
   memcpy(ab->data,wsp1,len*sizeof(complx));
}

/* this does A = B*A */
void cm_multo_rev(mv_complx * ba, mv_complx * b)
{
   static complx * wsp1=NULL;
   static int wsp1_id=0;
   static int wsp1_size=0;
   int N,len;
   
   N = b->row;
   len = N*N;
   wsp1 = get_complx_wsp(len, wsp1, &wsp1_size, &wsp1_id);
   
   zgemm_("N","N",&N,&N,&N,&CPLX1,b->data,&N,ba->data,&N,&CPLX0,wsp1,&N);
   memcpy(ba->data,wsp1,len*sizeof(complx));
}

void cmv_addto(mv_complx *ab, mv_complx * b)
{
   int len;
   
   len = (b->row)*(b->col);
   zaxpy_(&len,&CPLX1,b->data,&INTONE,ab->data,&INTONE);
}

void cmv_subfrom(mv_complx *ab, mv_complx * b)
{
   int len;
   
   len = (b->row)*(b->col);
   zaxpy_(&len,&CPLXm1,b->data,&INTONE,ab->data,&INTONE);
}

void cm_diag(mv_complx *A, mv_complx *eigs, mv_complx *T)
{
   int N, len, lwsp, info;
   double testwsp, vl[2];
   static mv_complx *wsp1=NULL, *wsp2=NULL;
   complx *tmpA;
   double *rwork;
      
   N = A->row;
   len = N*N;
   if ( (N!=A->col) || (N!=eigs->row) || (N!=T->row) || (N!=T->col) ) {
      fprintf(stderr,"cm_diag error: dimension mismatch\n");
      exit(1);
   }
   wsp1 = cmv_static(wsp1,len+N,1);
   tmpA = wsp1->data;
   rwork = (double*)(tmpA+len);
   memcpy(tmpA, A->data,len*sizeof(complx));
#ifdef LAPACK_NO_WSP_TEST   
   lwsp = 2*N;
#else 
   lwsp = -1;
   zgeev_("N","V",&N,tmpA,&N,eigs->data,vl,&INTONE,T->data,&N,&testwsp,&lwsp,rwork,&info); 
   lwsp = (int)testwsp;
#endif
   wsp2 = cmv_static(wsp2,lwsp,1);
   zgeev_("N","V",&N,tmpA,&N,eigs->data,vl,&INTONE,T->data,&N,wsp2->data,&lwsp,rwork,&info);
   if ( info != 0) {
      fprintf(stderr,"cm_diag error: diagonalization failed\n");
      exit(1);
   }
}

mv_complx * cm_ln(mv_complx * m)
{
   static mv_complx *T=NULL, *wsp1=NULL, *wsp2=NULL;
   double testwsp, *rwork;
   int N, len, lwsp, info;
   mv_complx *res;
   complx *A, *eigs, *vl;
      
   N = m->row;
   len = N*N;
   if ( m->col == 1 ) {
      /* it is a vector, easy job */
      res = complx_matrix_alloc(N,N);
      cmv_zero(res);
      for (lwsp=0; lwsp<N; lwsp++) {
         res->data[lwsp*(N+1)] = Clog(m->data[lwsp]);
      }
      return res;
   } else if ( N != m->col ) {
      /* error */
      fprintf(stderr,"cm_ln error: matrix is not square\n");
      exit(1);
   }
   /* do matrix stuff */
   res = complx_matrix_alloc(N,N);
   T = cmv_static(T,N,N);
   wsp1 = cmv_static(wsp1,len+2*N+1,1);
   A = wsp1->data;
   eigs = A + len;
   vl = eigs + N;
   rwork =(double*)(vl+1);
   memcpy(A,m->data,len*sizeof(complx));
#ifdef LAPACK_NO_WSP_TEST
   lwsp = 2*N;
#else
   lwsp = -1;
   zgeev_("N","V",&N,A,&N,eigs,vl,&INTONE,T->data,&N,&testwsp,&lwsp,rwork,&info); 
   lwsp = (int)testwsp;
#endif
   wsp2 = cmv_static(wsp2,lwsp,1);
   zgeev_("N","V",&N,A,&N,eigs,vl,&INTONE,T->data,&N,wsp2->data,&lwsp,rwork,&info);
   if ( info != 0) {
      fprintf(stderr,"cm_ln error: diagonalization failed\n");
      exit(1);
   }
   cmv_zero(res);
   for (lwsp=0; lwsp<N; lwsp++) {
     res->data[lwsp*(N+1)] = Clog(eigs[lwsp]);
   }
   simtrans(res,T);
   return res;
}

void dm_symdiag(mv_double *A, mv_double *eigs, mv_double *T)
{
   /* diagonalization of a real symmetric matrix (e.g. Hamiltonian) */
   int N,lwsp,info;
   double testwsp;
   static mv_double *mwsp=NULL;
   
   N = A->row;
   if ( (N != eigs->row) || (N != A->col) || (N != T->row) || (N != T->col) ) {
      fprintf(stderr,"dm_symdiag error: dimension missmatch\n");
      exit(1);
   }

   memcpy(T->data,A->data,N*N*sizeof(double));
   lwsp = -1;
   dsyev_("V","L",&N, T->data,&N,eigs->data,&testwsp,&lwsp,&info);
   lwsp = (int)testwsp;
   mwsp = dmv_static(mwsp,lwsp,1);
   dsyev_("V","L",&N, T->data,&N,eigs->data,mwsp->data,&lwsp,&info);
   if ( info != 0) {
      fprintf(stderr,"dm_symdiag error: diagonalization failed\n");
      exit(1);
   }

}



















