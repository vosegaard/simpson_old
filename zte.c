#include <stdio.h>
#include <stdlib.h>
#include "cm_new.h"
#include "pulse.h"

typedef struct _ZTEstruct {
   mv_complx *rho0, *sigma;
   mv_complx *prop;
   mv_int * idx;
} ZTEstruct;

ZTEstruct ZTE;
/*
ZTE.rho0 = NULL;
ZTE.sigma = NULL;
ZTE.prop = NULL;
ZTE.idx = NULL;
*/

/* make visible global variable puls */
extern Pulse* puls;

int tclZTEtestBegin(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
   int len;
   check_pulse();
   _evolve_with_prop();
   _reset_prop();
   /*
   if (ZTE.rho0 != NULL) complx_matrix_free(ZTE.rho0);
   if (ZTE.sigma != NULL) complx_vector_free(ZTE.sigma);
   if (ZTE.prop != NULL) complx_matrix_free(ZTE.prop);
   if (ZTE.idx != NULL) int_vector_free(ZTE.idx);
   */
   ZTE.rho0 = cmv_dup(puls->sigma);
   len = (puls->N)*(puls->N);
   ZTE.idx = int_vector_alloc(len);
   ZTE.sigma = complx_vector_alloc(len); /* these are maximal possible lengths */
   return TCL_OK;
}

int tclZTEtestEnd(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
   int i, j, len, N, n, *idx;
   double d, eps = 1e-5;
   complx *r0, *r1, *sg;
   int a, b, k, kk, l, ll;
      
   check_pulse();

   if (argc == 2) {
      if (Tcl_GetDouble(interp,argv[1],&eps) == TCL_ERROR) 
         return TclError(interp,"zte_test_end: argument 1 must be double <epsilon>");
   }

   _evolve_with_prop();
   r0 = ZTE.rho0->data;
   r1 = puls->sigma->data;
   sg = ZTE.sigma->data;
   N = puls->N;
   len = N*N;

   printf("ZTE report\n==========\n");
   printf("Hilbert space dim   = %d\n",N);
   printf("Liouville space dim = %d\n",len);

   idx = ZTE.idx->data;
   n = 0;
   for (i=0; i<len; i++) {
      d = fabs(r0->re - r1->re) + fabs(r0->im - r1->im);
      if (d > eps) {
         *idx = i;
	 idx++;
	 n++;
	 *sg = *r1;
	 sg++;
      }
      r0++;
      r1++;
   }
   if (n == 0) {
      printf("zte_test_end: there was no evolution, do nothing\n");
      return TCL_OK;
   }
   printf("ZTE space dim       = %d\n",n);
   
   ZTE.sigma->row = ZTE.idx->row = n;
   ZTE.prop = complx_matrix_alloc(n,n);
   cmv_zero(ZTE.prop);
   for (i=0; i<n; i++) {
      a = ZTE.idx->data[i];
      k = a/N;
      l = a-k*N;
      for (j=0; j<n; j++) {
         b = ZTE.idx->data[j];
	 kk = b/N;
	 ll = b-kk*N;
	 ZTE.prop->data[i+j*n] = Cmul((puls->U->data[l+ll*N]),Conj(puls->U->data[k+kk*N]));
      }
   }
   _reset_prop();
   

   return TCL_OK;
}

int tclZTEStore(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
   return TCL_OK;
}

int tclZTEAcq(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
   int np=1;
   int i, j, Nzte;
   complx *ptr, *z1, *z2;
   static mv_complx *y=NULL;
      
   check_pulse();

   if (argc == 2) {
      if (Tcl_GetInt(interp,argv[1],&np) == TCL_ERROR) 
         return TclError(interp,"zte_acq: argument 1 must be integer <np>");
   }
   if ( np < 1) 
      return TclError(interp,"zte_acq error: wrong number of acq points");
   if (puls->curr_nsig + np > LEN(puls->fid)) {
     fprintf(stderr,"zte_acq error: acq overflow in fid points\n");
     exit(1);
   }

   Nzte = ZTE.idx->row;
   y = cmv_static(y,Nzte,1);
   for (i=1; i<=np; i++) {
      ptr = &(puls->fid[++(puls->curr_nsig)]);
      *ptr = CPLX0;
      /* this is not optimal, ignores constant contributions to observable */
      z1 = ZTE.sigma->data;
      for (j=0; j<Nzte; j++) {
         z2 = puls->fdetect->data + ZTE.idx->data[j];
         ptr->re += z1->re*z2->re - z1->im*z2->im;
	 ptr->im += z1->re*z2->im + z1->im*z2->re;
	 z1++;
      }
      cm_mulmv(y,ZTE.prop,ZTE.sigma);
      cmv_copy(ZTE.sigma,y);
   }
   


   return TCL_OK;
}

void tclcmd_zte(Tcl_Interp* interp) {
Tcl_CreateCommand(interp,"zte_test_begin",tclZTEtestBegin,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateCommand(interp,"zte_test_end",tclZTEtestEnd,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateCommand(interp,"zte_store",tclZTEStore,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateCommand(interp,"zte_acq",tclZTEAcq,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

}
