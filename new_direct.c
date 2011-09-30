#include "new_direct.h"

/* make visible global variable pulse */
extern Pulse* puls;

#include <errno.h>
void new_direct_initialize(Tcl_Interp *interp,Sim *s)
{
   Tcl_Obj *data, **listObjElem,**listObjElem2;
   int Nelem, Nacqs, i;
   char buf[512];

   sprintf(buf,"if {[catch {pulseq_parser [info body %s]} __pulseq_parsed]} "
               "{puts \"new_direct_initialize error: $__pulseq_parsed\"\n exit} \n ",s->P->pulsename);
	       
   /* sprintf(buf,"pulseq_parser [info body %s]",s->P->pulsename); */
   /*printf("command = %s\n",buf);*/
   if ( Tcl_EvalEx(interp,buf,-1,TCL_EVAL_GLOBAL) != TCL_OK) {
      TclError(interp, "new_direct_initialize error: Unable to execute command 'pulseq_parser':\n");
      exit(1);
   }
/*printf("ini 1\n");*/
   /* data = Tcl_GetObjResult(interp); */
   data = Tcl_GetVar2Ex(interp,"__pulseq_parsed",NULL,TCL_GLOBAL_ONLY);

   if (Tcl_ListObjGetElements(interp, data, &Nelem, &listObjElem) != TCL_OK) { 
      TclError(interp,"new_direct_initialize error: cannot decompose list in elements");
      exit(1);
   }
/*printf("ini 2 \n");*/
   s->P->Npsq = Nelem-1;
   s->P->psq = listObjElem;
   /*
   s->P->psq = (Tcl_Obj**)Tcl_Alloc((s->P->Npsq)*sizeof(Tcl_Obj*));
   for (i=0; i<s->P->Npsq; i++) {
      s->P->psq[i] = Tcl_DuplicateObj(listObjElem[i]);
      Tcl_IncrRefCount(s->P->psq[i]);
   }
   */

   if (Tcl_ListObjGetElements(interp, listObjElem[Nelem-1], &Nacqs, &listObjElem2) != TCL_OK) { 
      TclError(interp,"new_direct_initialize error: cannot decompose acqs list in elements");
      exit(1);
   }
   s->P->pacqs = int_vector_alloc(Nacqs);
/*printf("ini 3\n");*/
   for (i=0; i<Nacqs; i++) {
      if (Tcl_GetIntFromObj(interp,listObjElem2[i],&Nelem) != TCL_OK) { 
         TclError(interp,"new_direct_initialize error: cannot get int from acqs list elements %i",i);
         exit(1);
      }
      s->P->pacqs->data[i] = Nelem;
   }
/*printf("ini 4\n");*/

   /*Tcl_Free((char*)listObjElem);*/
   //Tcl_Free((char *)listObjElem2);
   /* s->P->cacqs = (Tcl_Obj**)Tcl_Alloc(Nacqs*sizeof(Tcl_Obj*)); */
   s->P->cacqs = (Tcl_Obj**)malloc(Nacqs*sizeof(Tcl_Obj*));
   s->P->cacqs_pos = 0;
   s->P->acqnp = (int*)malloc(Nacqs*sizeof(int));
   s->P->acqph = (double*)malloc(Nacqs*sizeof(double));
   
   /* test output */
   /*
   printf("\nnew_direct_initialize test output\n=================================\n");
   for (i=0; i<s->P->Npsq; i++) printf("El. %i.) = \n '%s'\n",i,Tcl_GetString(s->P->psq[i]));
   printf("pacqs vector: (");
   for (i=0; i<s->P->pacqs->row; i++) printf(" %i",s->P->pacqs->data[i]);
   printf(")\n\n\n");
   */
/*Tcl_Eval(interp,"puts aaaaaHHHHHHHAAAAAAAA ");   
for (i=0; i<s->P->Npsq; i++) {
   printf("EXECUTING El. %i ... START \n",i);
   Tcl_EvalObjEx(interp, s->P->psq[i], TCL_EVAL_GLOBAL);
   sprintf(buf,"%s",Tcl_GetString(s->P->psq[i]));
   printf("%s\n",buf);
   Tcl_Eval(interp,buf);
   printf("--------> El. %i ... STOP \n\n",i);
}
*/
/*printf("ini 5\n");*/
   
   s->dw = 1.0e6/s->sw;
   if (s->sw1 > 0) s->dw1 = 1.0e6/s->sw1;
   if (s->wr > 0.001) {
      s->taur = 2.0e6*M_PI/s->wr;
   } else {
      s->taur = -100;
   }
   if ( (s->imethod == M_GCOMPUTE_NEW) || (s->imethod == M_GCOMPUTE2_NEW) ) {
      s->P->new_gcompute = 1;
   } else {
      s->P->new_gcompute = 0;
   }
   /* other stuff moved from pulse_propagate */
   s->P->H = s->H;
   s->P->fid = s->fid;
}

/* these settings might be valid for NEW_DIRECT only... */
void pulse_new_direct_sim_ini(Sim *s, double t_ini)
{
   Pulse *P;
   int i;
   
   P = s->P;
   if (!(puls->new_gcompute)) {
      P->curr_nsig=0;
      m_zerov(P->fid);  /* don't do this if NEW_GCOMPUTE */
   }
   cmv_copy(P->sigma,P->fstart);
   memset(P->phv,0,(1+P->nchan)*sizeof(double)); 
   memset(P->rfv,0,(1+P->nchan)*sizeof(double)); 

  if (s->H->isdiag) {
     /* limit initial propagator to diagonal, rest of allocated space is not used
        but ready to be used when propagator turns to full form due to pulse    */
     P->U->col = 1;
     cv_d1(P->U);
  } else {
     cm_unit(P->U);
  }
  P->Uisunit=1;
  P->isselectivepulse=0;
  P->cannotbestored=0;
  P->waspulse=0;
  P->brl=P->brlorig;
  if (P->wr == 0.0)
    P->dtmax=1e9; /* dtmax is infinity in the static case. */
  else
    P->dtmax= 1.0e-6; /* dt max is 1 usec per default in case of spining*/
  P->dt= 1e99;
  P->tpulsestart_usec=t_ini*1.0e6;
  P->t_usec=P->tpulsestart_usec;
  P->tpropstart_usec=P->t_usec;

  /* every interaction but the DC offset is turned on */
  ham_turnon(P->H,"all");
  /* ZT: take care of static field inhomogeneity ofsets */
  if ( !(P->inhom_offset) ) {
     /* offset is specifically turned off */
     ham_turnoff(P->H,"offset");
  } else {
     /* add these offsets on all channels */
     memset(P->tmpIz,0,P->N*sizeof(double));
     for (i=1; i<=P->nchan; i++) {
        daxpy_(&(P->N),&(P->inhom_offset[i]),P->chan_Iz[i],&INTONE,P->tmpIz,&INTONE);
	/* printf("chan %i offset %f\n",i,P->inhom_offset[i]); */
     }
     ham_set_offset(P->H,P->tmpIz,2);
  }
}

void new_direct_destroy(Sim *s)
{
   int i;

 /* don't know why this does not work... */   
 /*  for (i=0; i<s->P->pacqs->row; i++) {
      Tcl_Free((char*)(s->P->cacqs[i]));
   }
   Tcl_Free((char*)(s->P->cacqs)); */
   /*free(s->P->cacqs);
   free(s->P->acqnp);
   free(s->P->acqph);
   int_vector_free(s->P->pacqs);
   Tcl_Free((char*)(s->P->psq)); */
}

void new_direct_modnumber(double taur,double dt, int *k, int *l)
{
   double a,b;
   int m,n;
   
   a = taur;
   b = dt;
   m = 1;
   n = 1;
   if (taur>dt) {
      do {
         b += dt;
	 n++;
	 if (b-a > 1e-6) {
	    a += taur;
	    m++;
	 }
      } while ( fabs(a-b) > 1e-6 );
   } else {
      do {
         a += taur;
	 m++;
	 if (a-b > 1e-6) {
	    b += dt;
	    n++;
	 }
      } while ( fabs(a-b) > 1e-6 );
   }
   *k = m;
   *l = n; 
}

void new_direct_acq(Sim *s, int which_acq)
{
   int np, Ndim, i, k, l, m;
   double phase, dw, t0, dt, dt2;
   Tcl_Obj *code;
   complx *fidptr;
   static mv_complx *rho=NULL, *T=NULL, *Ud=NULL, *det=NULL;
   const double dmone=-1.0;
   double *dptr;
   const int itwo=2;
   
   /* note: all reference to 'pulse' can be replaced by 's->P' */
   np = puls->acqnp[which_acq];
   if ( np < 0 ) np = s->np; 
   if (puls->curr_nsig + np > LEN(puls->fid)) {
     fprintf(stderr,"acq_block error: overflow in fid points\n");
     exit(1);
   }
   code = puls->cacqs[which_acq];
   phase = (puls->acqph[which_acq])*DEG2RAD;
   dw = s->dw;
   t0 = puls->t_usec;

   /* acquire first point now */
   _evolve_with_prop();
   _reset_prop();
   fidptr = &(puls->fid[++(puls->curr_nsig)]);
   dptr = (double*)fidptr + 1;
   /* use this trace since it is faster, adjoint back at the end */
   /* warning: works only if dens. matrix is Hermitian           */
   *fidptr = cm_trace_adjoint(puls->fdetect,puls->sigma);

   /* this is to initialize values for check_dwelltime */
   puls->new_direct_mxpos = NEW_DIRECT_MXPOS_INI;
   puls->check_dwelltime_t0 = puls->t_usec;

   /* evaluate acq_block to get its propagator */
   if (Tcl_EvalObjEx(s->P->interp, code, TCL_EVAL_GLOBAL) != TCL_OK) { 
      TclError(s->P->interp,"new_direct_acq error while executing acq_block (%i):\n"
               "'%s'\n", which_acq, Tcl_GetString(code));
      exit(1);
   }
   /*printf("exec code = '%s'\n", Tcl_GetString(code));*/
   dt = puls->t_usec - t0;
   /*printf("=> %f + %f\n",t0, dt);   */
   
   Ndim = s->matdim;
   Ud = cmv_static(Ud, Ndim, 1);
   rho = cmv_static(rho, Ndim, Ndim);
   det = cmv_static(det, Ndim, Ndim);
   
   if (s->taur < 0) {
   /* static case */
      /* at the moment, continue only when acq_block exactly synchronized with dwell time */
      if ( fabs(dt-dw)>1.0e-6 ) {
         fprintf(stderr,"new_direct_acq error: acq_block (%i) not synchronized with dwell time\n",which_acq);
         exit(1);
      }
      /* diagonalize, transform dens. matrix and detection op. */
      if (puls->U->col == 1) {
         cmv_copy(Ud, puls->U);
	 cmv_copy(rho,puls->sigma);
	 cmv_copy(det,puls->fdetect);
      } else {
         T = cmv_static(T, Ndim, Ndim);
         cm_diag(puls->U, Ud, T);
	 cmv_copy(rho,puls->sigma);
	 simtransh_adj(rho,T);
	 cmv_copy(det,puls->fdetect);
	 simtrans_adj(det,T);
      }
      /* evolve and acquire */
      for (i=1; i<np; i++) {
         simtrans_diagprop3(rho, Ud);
	 fidptr++;
         *fidptr = cm_trace_adjoint(det,rho);
      }
      if (puls->acq_adjoint == 0) {
         /* adjoint back the fid */
         dscal_(&np,&dmone,dptr,&itwo);
      }
   } else {
   /* MAS */
      if ( fabs(dw-dt) < 1.0e-6 ) {
         new_direct_modnumber(s->taur,dt,&k,&l);
         /*printf("modulation numbers are %i*tau_r / %i*dt\n",k,l);*/
	 k = l;
      } else if ( dw-dt > 1.0e-6 ) {
         m = (int)floor(dw/dt+1e-6);
         dt2 = dt*(double)m;
         /*printf("---> dw = %f, dt = %f, dt2 = %f, m = dw/dt = %i\n",dw,dt,dt2,m);*/
         if ( fabs(dw-dt2) > 1e-6 ) {
            fprintf(stderr,"new_direct_acq error: acq_block (%i) doesn't fit to dwell time\n",which_acq);
	    exit(1);
         }
         new_direct_modnumber(s->taur,dt2,&k,&l);
         /*printf("modulation numbers are %i*tau_r / %i*%i*dt\n",k,l,m);*/
	 k = l;
	 l *= m;
      } else {
         m = (int)floor(dt/dw+1e-6);
         dt2 = dw*(double)m;
         /*printf("---> dw = %f, dt = %f, dt2 = %f, m = dt/dw = %g = %i\n",dw,dt,dt2,floor(dt/dw),m);*/
         if ( fabs(dt-dt2) > 1e-6 ) {
            fprintf(stderr,"new_direct_acq error: acq_block (%i) doesn't fit to dwell time\n",which_acq);
	    exit(1);
         }
         new_direct_modnumber(s->taur,dt,&k,&l);
         /*printf("modulation numbers are %i*tau_r / %i*dt\n",k,l);*/
	 k = m*l;
      }
      /* finish l-times acq_block, number of propagators is k */
      for (i=1; i<l; i++) {
         tclCheckDwelltime(NULL,s->P->interp, 0, NULL);
	 /*printf("+++> run %i, check_dwelltime passed\n",i);*/
         if (Tcl_EvalObjEx(s->P->interp, code, TCL_EVAL_GLOBAL) != TCL_OK) { 
            TclError(s->P->interp,"new_direct_acq error while executing acq_block (%i), run %i:\n"
                  "'%s'\n", which_acq, i, Tcl_GetString(code));
            exit(1);
         }
	 /*printf("exec code = '%s'\n", Tcl_GetString(code));*/
      }
      tclCheckDwelltime(NULL,s->P->interp, 0, NULL);
      /*printf("last chech_dwelltime passed\n\n\n");*/
      /* redefine k,l: initial and final index of propagators */
      l = NEW_DIRECT_MXPOS_INI + k;
      k = NEW_DIRECT_MXPOS_INI;
      /* test if all propagators were generated */
      for (i=k; i<l; i++) {
         if (puls->STO[i] == NULL) {
	    fprintf(stderr,"new_direct_acq error: found only %i propagators out of %i\n",i-k,l-k);
	    exit(1);
	 }
	 /* printf("prop check: slot %i, start time %g\n",i,puls->STO_tpropstart_usec[i]);
	   cm_print(puls->STO[i],"prop"); */
      }
      l--;
      /* prepare transformed matrices ... */
      if (puls->acq_adjoint == 0) {
         fidptr->im *= -1.0;
	 cmv_copy(det, puls->fdetect);
      } else {
         cm_adjoint(det, puls->fdetect);
      }
      if ( puls->U->col == 1 ) {
         /* ... for diagonal props */
	 cmv_copy(rho,puls->sigma);
	 puls->matrix[k] = cmv_dup(det);
	 for (i=k; i<l; i++) {
	    puls->matrix[i+1] = cmv_dup(det);
	    simtrans_adj_diagprop3(puls->matrix[i+1],puls->STO[i]);
	    complx_matrix_free(puls->STO[i]);
	    puls->STO[i]=NULL;
	 }
         cmv_copy(Ud,puls->STO[l]);;
         complx_matrix_free(puls->STO[l]);
         puls->STO[l]=NULL;
      } else {
         /* ... for full props */
         T = cmv_static(T, Ndim, Ndim);
         cm_diag(puls->STO[l], Ud, T);
         complx_matrix_free(puls->STO[l]);
         puls->STO[l]=NULL;
	 cmv_copy(rho,puls->sigma);
	 simtransh_adj(rho,T);
	 puls->matrix[k] = cmv_dup(det);
	 simtrans_adj(puls->matrix[k],T);
	 for (i=k; i<l; i++) {
	    puls->matrix[i+1] = cmv_dup(det);
	    simtrans_adj(puls->matrix[i+1],puls->STO[i]);
	    simtrans_adj(puls->matrix[i+1],T);
	    complx_matrix_free(puls->STO[i]);
	    puls->STO[i]=NULL;
	 }
      }
      /* run the acquisition */
      m = k+1;
      for (i=1; i<np; i++) {
	 fidptr++;
         *fidptr = cm_trace_adjoint(rho,puls->matrix[m]);
         m++;
	 if ( m > l ) {
	    /* period finished */
	    m = k;
	    simtrans_diagprop3(rho,Ud);
         }
      }
   }

   if ( phase != 0.0 ) {
      complx z;
      
      fidptr = &(puls->fid[puls->curr_nsig]);
      z = Complx(cos(phase),sin(phase));
      zscal_(&np, &z, fidptr, &INTONE);
   }
   
   puls->curr_nsig += np-1;

}

void new_direct_calc(Sim *s, double gamma_add, mv_complx *wfidsum, mv_complx *wfid, int *nfid)
{
   Tcl_Obj *code;
   int ig,ic, c1;
   double t_ini;

   /*
   printf("\nTest acq_block initiated\n");
   code = s->P->psq[s->P->pacqs->data[0]];
   if (Tcl_EvalObjEx(s->P->interp, code, TCL_EVAL_GLOBAL) != TCL_OK) { 
      TclError(s->P->interp,"new_direct_calc error: cannot execute acq_block code");
      exit(1);
   }
   */
      
   /* will work only for 1D */
   if (s->ni != 0) {
      fprintf(stderr,"Sorry, direct_new not supported for 2D simulations yet...\n");
      exit(1);
   }
   /* will work only for Hermitian density matrix */
   if (!cm_ishermit(s->P->fstart)) {
      fprintf(stderr,"Sorry, direct_new not supported for non-Hermitian density matrix...\n");
      exit(1);
   }

   if (puls) {
     fprintf(stderr,"new_direct_calc error: pulseq propagation was called before without destroying\n");
     exit(1);
   }
   puls = s->P;

   for (ig=0;ig<s->ngamma;ig++) {
/*printf("new_direct_calc ig=%d \n",ig);*/
      if ( fabs(s->wr) < 0.001 ) {
        t_ini=0;
      } else {
	t_ini = ig/(double)s->ngamma*2.0*M_PI/s->wr;
      }
      s->H->gamma = 360*ig/(double)s->ngamma + gamma_add;
      pulse_new_direct_sim_ini(s, t_ini);
      c1 = 0;
      for (ic=0; ic<s->P->Npsq; ic++ ) {
/*printf("gugu ic=%d\n",ic);*/
         code = s->P->psq[ic];
/*printf("gDDDgu \n");*/
         if (Tcl_EvalObjEx(s->P->interp, code, TCL_EVAL_GLOBAL) != TCL_OK) { 
            TclError(s->P->interp,"new_direct_calc error while executing block:\n"
	              "'%s'\n",Tcl_GetString(code));
            exit(1);
         }
/*printf("   ic = %d; (%s)\n",ic,Tcl_GetString(code));*/
	 if ( ic == puls->pacqs->data[c1] ) {
	    /* code was acq_block, take care of it */
	    new_direct_acq(s, c1);
	    if (++c1 == s->P->pacqs->row) c1=0;
	 }
      }
      cmv_addto(wfidsum,wfid);
      (*nfid)++;
   }

   puls = NULL;
   
   /*printf("don't STOP here.\n");*/
   /* exit(1); */
}













void new_gcompute_calc(Sim *s,double gamma_add, mv_complx *wfidsum, mv_complx *wfid, int *nfid)
{
   Tcl_Obj *code;
   int ig, k, l, m, c1, np, Ndim, q, ir, im, iq;
   double t_ini, dtg, dt2, phase;
   complx *fidptr, z;
   static mv_complx *T=NULL, *Ud=NULL, *det=NULL;
   
   /*printf("gcompute_new started\n");*/
   
   if (s->ngamma <= 1) {
      fprintf(stderr,"error: gcompute_new called with gamma_angles = %d\n       (must be > 1)\n",s->ngamma);
      exit(1);
   }

   /* will work only for 1D */
   if (s->ni != 0) {
      fprintf(stderr,"Sorry, gcompute_new not supported for 2D simulations yet...\n");
      exit(1);
   }
   /* multiple acquisition blocks not allowed */
   if ( s->P->pacqs->row != 1 ) {
      fprintf(stderr,"error: gcompute_new not compatible with multiple acquisition blocks\n");
      exit(1);
   }
   /* first check on timings */
   dtg = s->taur/s->ngamma;
   m = (int)floor(s->dw/dtg+1e-6);
   dt2 = dtg*(double)m;
   if (fabs(dt2 - s->dw)>1e-6) {
      fprintf(stderr,"new_gcompute_calc error: bad synchronization of acquisition and gamma-averaging\n");
      exit(1);
   }
	 
   if (puls) {
     fprintf(stderr,"new_gcompute_calc error: pulseq propagation was called before without destroying\n");
     exit(1);
   }
   puls = s->P;

   c1 = puls->pacqs->data[0];
   /*printf("c1=%d\n",c1);*/
   cmv_zero(wfid);
   puls->curr_nsig = 1;
   fidptr = &(puls->fid[puls->curr_nsig]);
   Ndim = s->matdim;
   Ud = cmv_static(Ud, Ndim, 1);
   T = cmv_static(T, Ndim, Ndim);
   det = cmv_static(det, Ndim, Ndim);

   k = NEW_DIRECT_MXPOS_INI;
   for (ig=0;ig<s->ngamma;ig++) {
      t_ini = ig/(double)s->ngamma*2.0*M_PI/s->wr;
      s->H->gamma = 360*ig/(double)s->ngamma + gamma_add;
      pulse_new_direct_sim_ini(s, t_ini);
      if ( c1 == 0 ) {
         /* there is no preparation pulse sequence */
	 /*printf("there is no preparation\n");*/
	 if (puls->matrix[k+ig] != NULL) complx_matrix_free(puls->matrix[k+ig]); 
         puls->matrix[k+ig] = cmv_dup(puls->sigma);
      } else {
         /*printf("executing preparation\n");*/
	 /* evaluate preparation pulse sequence */
	 code = puls->psq[0];
         if (Tcl_EvalObjEx(puls->interp, code, TCL_EVAL_GLOBAL) != TCL_OK) { 
            TclError(s->P->interp,"new_gcompute_calc error while executing block:\n"
	              "'%s'\n",Tcl_GetString(code));
            exit(1);
         }
         _evolve_with_prop();
         _reset_prop();
	 if (puls->matrix[k+ig] != NULL) complx_matrix_free(puls->matrix[k+ig]); 
         puls->matrix[k+ig] = cmv_dup(puls->sigma);
      }
      /* evaluate acquisition block */
      if (ig == 0) {
         /*printf("executing acq_block\n");*/
         code = puls->psq[c1];
         if (Tcl_EvalObjEx(puls->interp, code, TCL_EVAL_GLOBAL) != TCL_OK) { 
            TclError(s->P->interp,"new_gcompute_calc error while executing block:\n"
                 "'%s'\n",Tcl_GetString(code));
            exit(1);
         }
	 /* now all info for this check is available */
	 np = puls->acqnp[0];
         /*printf("   np=%d\n",np);*/
         if ( np < 0 ) np = s->np; 
         if ( np > LEN(puls->fid)) {
           fprintf(stderr,"acq_block error: overflow in fid points\n");
           exit(1);
         }
         phase = puls->acqph[0];
         /* it is here since fdetect could possibly be changed in acq_block */
         if (puls->acq_adjoint == 0) {
            cmv_copy(det, puls->fdetect);
         } else {
            cm_adjoint(det, puls->fdetect);
         }
      }
      dt2 = puls->t_usec;
      code = puls->cacqs[0];
      /*printf("executing acq_block code for gamma index %d\n",ig);*/
      if (Tcl_EvalObjEx(s->P->interp, code, TCL_EVAL_GLOBAL) != TCL_OK) { 
         TclError(s->P->interp,"new_gcompute_calc error while executing code from acq_block (%i), run %i:\n"
               "'%s'\n", 0, ig, Tcl_GetString(code));
         exit(1);
      }
      /* final check on timings */
      if (fabs(puls->t_usec - dt2 - dtg) > 1e-6) {
         fprintf(stderr,"new_gcompute_calc error: acq_block not synchronized with gamma-averaging\n");
	 fprintf(stderr,"   acq_block length = %g usec\n",puls->t_usec-dt2);
	 fprintf(stderr,"   gamma time shift = %g usec\n",dtg);
	 fprintf(stderr,"   gamma angle idx  = %d usec\n",ig);
	 exit(1);
      } 
      /* here we acquire first point of fid */
      z = cm_trace_adjoint(puls->matrix[k+ig],det);
      fidptr->re += z.re;
      fidptr->im += z.im;
      /* store cummulative propagators */
      l = k + ig;
      if (puls->STO[l] != NULL) complx_matrix_free(puls->STO[l]);
      if ( ig == 0 ) {
         puls->STO[l] = cmv_dup(puls->U);
      } else {
         puls->STO[l] = complx_matrix_alloc(puls->U->row,puls->U->col);
         if ( puls->U->col == 1) {
	    /* all props are diagonal */
	    cmv_mul_elem(puls->STO[l],puls->U,puls->STO[l-1]);
	 } else {
	    cm_mul(puls->STO[l],puls->U,puls->STO[l-1]);
	 }
      }
      /* finish with pulseq code if there is any */
      if (puls->Npsq > c1+1 ) {
         /*printf("executing last pulseq block\n");*/
	 code = puls->psq[c1+1];
         if (Tcl_EvalObjEx(puls->interp, code, TCL_EVAL_GLOBAL) != TCL_OK) { 
            TclError(s->P->interp,"new_gcompute_calc error while executing block:\n"
                 "'%s'\n",Tcl_GetString(code));
            exit(1);
         }
      
      }

   } /* end of ig loop */
   /* transformations of dens. matrix and detect op. */
   l = k + s->ngamma;
   puls->matrix[l] = cmv_dup(det);
   if ( puls->STO[k]->col == 1 ) {
      /* ... for diagonal props */
      for (ig=1; ig<s->ngamma; ig++) {
         q = k+ig-1;
         simtrans_adj_diagprop3(puls->matrix[k+ig],puls->STO[q]);
         if ( puls->matrix[l+ig] != NULL ) complx_matrix_free(puls->matrix[l+ig]);
	 puls->matrix[l+ig] = cmv_dup(det);
	 simtrans_adj_diagprop3(puls->matrix[l+ig],puls->STO[q]);
	 complx_matrix_free(puls->STO[q]);
	 puls->STO[q] = NULL;
      }
      cmv_copy(Ud,puls->STO[l-1]);
      complx_matrix_free(puls->STO[l-1]);
      puls->STO[l-1] = NULL;
   } else {
      /* ... for full props */
      cm_diag(puls->STO[l-1], Ud, T);
      complx_matrix_free(puls->STO[l-1]);
      puls->STO[l-1] = NULL;
      simtransh_adj(puls->matrix[k],T);
      simtrans_adj(puls->matrix[l],T);
      for (ig=1; ig<s->ngamma; ig++) {
         q = k+ig-1;
         simtransh_adj(puls->matrix[k+ig],puls->STO[q]);
	 simtrans_adj(puls->matrix[k+ig],T);
         if ( puls->matrix[l+ig] != NULL ) complx_matrix_free(puls->matrix[l+ig]);
	 puls->matrix[l+ig] = cmv_dup(det);
	 simtrans_adj(puls->matrix[l+ig],puls->STO[q]);
	 simtrans_adj(puls->matrix[l+ig],T);
	 complx_matrix_free(puls->STO[q]);
	 puls->STO[q] = NULL;
      }
   }
   /* acquisition */
   ir = s->ngamma - 1;
   for (q=1; q<np; q++) {
      fidptr++;
      iq = q % s->ngamma;
      if ( m>1 ) {
         iq *= m;
	 iq = iq % s->ngamma;
      }
      for (im=0; im<m; im++) {
         simtrans_diagprop3(puls->matrix[k+ir],Ud);
         if (--ir < 0) ir = s->ngamma-1;
      }
      for (ig=0; ig<s->ngamma; ig++) {
         z = cm_trace_adjoint(puls->matrix[k+ig],puls->matrix[l+iq]);
         fidptr->re += z.re;
	 fidptr->im += z.im;
	 if (++iq == s->ngamma) iq = 0;
      }
   }
   puls->curr_nsig += np-1;
   *nfid = s->ngamma;

   puls = NULL;

   cmv_copy(wfidsum,wfid);

   if (phase != 0.0 ) {
      phase *= DEG2RAD;
      z = Complx(cos(phase),sin(phase));
      cmv_mulc(wfidsum,z);
   }
   
   
   /*printf("gcompute_new finished\n");*/

}


void new_gcompute2_calc(Sim *s,double gamma_add, mv_complx *wfidsum, mv_complx *wfid, int *nfid)
{
   Tcl_Obj *code;
   int ig, k, l, m, c1, np, Ndim, p, q, ir, im, iq;
   double t_ini, dtg, dt2, phase;
   complx *fidptr, z;
   static mv_complx *T=NULL, *Ud=NULL, *det=NULL;
   int Ng = s->ngamma;
   mv_complx **mx1, **mx2, **mx3;
   
   //printf("gcompute2_new started\n");
   
   if (Ng <= 1) {
      fprintf(stderr,"error: gcompute_new called with gamma_angles = %d\n       (must be > 1)\n",Ng);
      exit(1);
   }

   /* will work only for 1D */
   if (s->ni != 0) {
      fprintf(stderr,"Sorry, gcompute_new not supported for 2D simulations yet...\n");
      exit(1);
   }
   /* multiple acquisition blocks not allowed */
   if ( s->P->pacqs->row != 1 ) {
      fprintf(stderr,"error: gcompute_new not compatible with multiple acquisition blocks\n");
      exit(1);
   }
   /* first check on timings */
   dtg = s->taur/Ng;
   m = (int)floor(s->dw/dtg+1e-6);
   dt2 = dtg*(double)m;
   if (fabs(dt2 - s->dw)>1e-6) {
      fprintf(stderr,"new_gcompute_calc error: bad synchronization of acquisition and gamma-averaging\n");
      exit(1);
   }
	 
   if (puls) {
     fprintf(stderr,"new_gcompute_calc error: pulseq propagation was called before without destroying\n");
     exit(1);
   }
   puls = s->P;

   c1 = puls->pacqs->data[0];
   /*printf("c1=%d\n",c1);*/
   cmv_zero(wfid);
   puls->curr_nsig = 1;
   /* fidptr = &(puls->fid[puls->curr_nsig]); */
   fidptr = wfidsum->data;
   Ndim = s->matdim;
   Ud = cmv_static(Ud, Ndim, 1);
   T = cmv_static(T, Ndim, Ndim);
   det = cmv_static(det, Ndim, Ndim);

   mx1 = &(puls->matrix[NEW_DIRECT_MXPOS_INI]);
   mx2 = mx3 = &(puls->STO[NEW_DIRECT_MXPOS_INI]);
   for (ig=0;ig<Ng;ig++) {
      t_ini = ig/(double)Ng*2.0*M_PI/s->wr;
      s->H->gamma = 360*ig/(double)Ng + gamma_add;
      pulse_new_direct_sim_ini(s, t_ini);
      if ( c1 == 1 ) {
         /*printf("executing preparation\n");*/
	 /* evaluate preparation pulse sequence */
	 code = puls->psq[0];
         if (Tcl_EvalObjEx(puls->interp, code, TCL_EVAL_GLOBAL) != TCL_OK) { 
            TclError(s->P->interp,"new_gcompute_calc error while executing block:\n"
	              "'%s'\n",Tcl_GetString(code));
            exit(1);
         }
         _evolve_with_prop();
         _reset_prop();
      }
      /*printf("pointer to matrix pointer %p\n",mx1);*/
      if (*mx1 != NULL) complx_matrix_free(*mx1); 
      *mx1 = cmv_dup(puls->sigma);
      /* evaluate acquisition block */
      if (ig == 0) {
         /*printf("executing acq_block\n");*/
         code = puls->psq[c1];
         if (Tcl_EvalObjEx(puls->interp, code, TCL_EVAL_GLOBAL) != TCL_OK) { 
            TclError(s->P->interp,"new_gcompute_calc error while executing block:\n"
                 "'%s'\n",Tcl_GetString(code));
            exit(1);
         }
	 /* now all info for this check is available */
	 np = puls->acqnp[0];
         /*printf("   np=%d\n",np);*/
         if ( np < 0 ) np = s->np; 
         if ( np > LEN(puls->fid)) {
           fprintf(stderr,"acq_block error: overflow in fid points\n");
           exit(1);
         }
         phase = puls->acqph[0];
         /* it is here since fdetect could possibly be changed in acq_block */
         if (puls->acq_adjoint == 0) {
            cmv_copy(det, puls->fdetect);
         } else {
            cm_adjoint(det, puls->fdetect);
         }
      }
      dt2 = puls->t_usec;
      code = puls->cacqs[0];
      /*printf("executing acq_block code for gamma index %d\n",ig);*/
      if (Tcl_EvalObjEx(s->P->interp, code, TCL_EVAL_GLOBAL) != TCL_OK) { 
         TclError(s->P->interp,"new_gcompute_calc error while executing code from acq_block (%i), run %i:\n"
               "'%s'\n", 0, ig, Tcl_GetString(code));
         exit(1);
      }
      /* final check on timings */
      if (fabs(puls->t_usec - dt2 - dtg) > 1e-6) {
         fprintf(stderr,"new_gcompute_calc error: acq_block not synchronized with gamma-averaging\n");
	 fprintf(stderr,"   acq_block length = %g usec\n",puls->t_usec-dt2);
	 fprintf(stderr,"   gamma time shift = %g usec\n",dtg);
	 fprintf(stderr,"   gamma angle idx  = %d usec\n",ig);
	 exit(1);
      } 
      /* here we acquire first point of fid */
      z = cm_trace_adjoint(*mx1,det);
      fidptr->re += z.re;
      fidptr->im += z.im;
      mx1++;
      /* store cummulative propagators */
      if (*mx2 != NULL) complx_matrix_free(*mx2);
      if ( ig == 0 ) {
         *mx2 = cmv_dup(puls->U);
	 mx2++;
      } else {
         *mx2 = complx_matrix_alloc(puls->U->row,puls->U->col);
         if ( puls->U->col == 1) {
	    /* all props are diagonal */
	    cmv_mul_elem(*mx2,puls->U,*mx3);
	 } else {
	    cm_mul(*mx2,puls->U,*mx3);
	 }
	 mx2++;
	 mx3++;
      }
      /* finish with pulseq code if there is any */
      if (puls->Npsq > c1+1 ) {
         /*printf("executing last pulseq block\n");*/
	 code = puls->psq[c1+1];
         if (Tcl_EvalObjEx(puls->interp, code, TCL_EVAL_GLOBAL) != TCL_OK) { 
            TclError(s->P->interp,"new_gcompute_calc error while executing block:\n"
                 "'%s'\n",Tcl_GetString(code));
            exit(1);
         }
      
      }

   } /* end of ig loop */
   
   /* transformations of dens. matrix and detect op. */
   mx1 = &(puls->matrix[NEW_DIRECT_MXPOS_INI]);
   mx2 = mx1 + Ng;
   if (*mx2 != NULL) complx_matrix_free(*mx2);
   *mx2 = cmv_dup(det);
   mx3 = &(puls->STO[NEW_DIRECT_MXPOS_INI]);
   if ( (*mx3)->col == 1 ) {
      /* ... for diagonal props */
      for (ig=1; ig<Ng; ig++) {
         mx1++;
	 mx2++;
         q = k+ig-1;
         simtrans_adj_diagprop3(*mx1,*mx3);
         if ( *mx2 != NULL ) complx_matrix_free(*mx2);
	 *mx2 = cmv_dup(det);
	 simtrans_adj_diagprop3(*mx2,*mx3);
	 complx_matrix_free(*mx3);
	 *mx3 = NULL;
	 mx3++;
      }
      cmv_copy(Ud,*mx3);
      complx_matrix_free(*mx3);
      *mx3 = NULL;
   } else {
      /* ... for full props */
      l = NEW_DIRECT_MXPOS_INI + Ng -1;
      cm_diag(puls->STO[l], Ud, T);
      complx_matrix_free(puls->STO[l]);
      puls->STO[l] = NULL;
      simtransh_adj(*mx1,T);
      simtrans_adj(*mx2,T);
      for (ig=1; ig<Ng; ig++) {
	 mx1++;
	 mx2++;
         simtransh_adj(*mx1,*mx3);
	 simtrans_adj(*mx1,T);
         if ( *mx2 != NULL ) complx_matrix_free(*mx2);
	 *mx2 = cmv_dup(det);
	 simtrans_adj(*mx2,*mx3);
	 simtrans_adj(*mx2,T);
	 complx_matrix_free(*mx3);
	 *mx3 = NULL;
	 mx3++;
      }
   }
   /*printf("transformations done\n");*/
   
   /* averaging over gamma angles */
   mx1 = &(puls->STO[NEW_DIRECT_MXPOS_INI]); /* result of averaging */
   complx *z1, *z2, *z3, *z4, *z5;
   const int len = Ndim*Ndim;
   mx2 = &(puls->matrix[NEW_DIRECT_MXPOS_INI + Ng]);
   for (l=0; l<Ng; l++) {
      *mx1 = complx_matrix_alloc(Ndim,Ndim);
      cmv_zero(*mx1);
      if (l>0) {
         simtrans_adj_diagprop3(*mx2,Ud);
	 mx2++;
      }
      for (ig=0; ig<Ng; ig++){
         z1 = (*mx1)->data;
	 z2 = z4 = puls->matrix[NEW_DIRECT_MXPOS_INI + ig]->data;
	 z5 = z2 + len;
	 z3 = puls->matrix[NEW_DIRECT_MXPOS_INI + Ng + ((ig+l)%Ng)]->data;
         for (k=0; k<len; k++) {
	    z1->re += z2->re*z3->re - z2->im*z3->im;
	    z1->im += z2->re*z3->im + z2->im*z3->re;
	    z1++;
	    z3++;
	    z2 += Ndim;
	    if (z2 >= z5) z2 = ++z4;
	 }
      }
      mx1++;
   }
   /*printf("averaging done\n");*/

   /* acquisition */
   cmv_zero(T);
   zgerc_(&Ndim,&Ndim,&CPLX1,Ud->data,&INTONE,Ud->data,&INTONE,T->data,&Ndim);
   mx1 = &(puls->STO[NEW_DIRECT_MXPOS_INI]);
   simtrans_adj_diagprop3(*mx1,Ud);
   for (l=1; l<np; l++) {
      z1 = T->data;
      q = l % Ng;
      mx2 = mx1 + q; 
      z2 = (*mx2)->data;
      fidptr++;
      for (k=0; k<len; k++) {
         double d1 = z2->re;
	 double d2 = z2->im;
	 fidptr->re += d1;
	 fidptr->im += d2;
	 z2->re = d1*z1->re + d2*z1->im;
	 z2->im = -d1*z1->im + d2*z1->re;
	 z1++;
	 z2++;
      }
   }


   puls->curr_nsig += np-1;
   *nfid = Ng;

   puls = NULL;

   if (phase != 0.0 ) {
      phase *= DEG2RAD;
      z = Complx(cos(phase),sin(phase));
      cmv_mulc(wfidsum,z);
   }
   
   
   /*printf("gcompute2_new finished\n");*/

}





int tclAcqBlock(ClientData data,Tcl_Interp* interp, int objc, Tcl_Obj *objv[])
{
   int i, np=-1;
   double phase=0.0;
   Tcl_Obj *cacq;
   char *s;

   /* do nothing when all acq_blocks were already called once */
   if (puls->cacqs_pos >= puls->pacqs->row) return TCL_OK;

   /*
   for (i=1; i<objc; i++) printf("acq_block par %i = '%s'\n",i,Tcl_GetString(objv[i]));
   */
   
   check_pulse();
   
   if ( objc%2 != 0 ) 
      return TclError(interp,"acq_block argument count error. usage: \n\t"
                             "acq_block -np X -ph x { pulse code }");
   for (i=1; i<objc-1; i++) {
      s = Tcl_GetString(objv[i]);
/* printf("   analyzing '%s'\n",s); */
      i++;
      if (!strncmp(s,"-np",3)) {
         if ( Tcl_GetIntFromObj(interp,objv[i],&np) != TCL_OK ) 
	    return TclError(interp,"acq_block error: can not get int for -np parameter");
      } else if (!strncmp(s,"-ph",3)) {
         if ( Tcl_GetDoubleFromObj(interp,objv[i],&phase) != TCL_OK ) {
	    /* perhaps not a number... */
	    s = Tcl_GetString(objv[i]);
            if (!strcasecmp(s,"-X")) {
               phase=180;
            } else if (!strcasecmp(s,"-Y")) {
               phase=270;    
            } else if (!strcasecmp(s,"X")) {
               phase=0;
            } else if (!strcasecmp(s,"Y")) {
               phase=90;
            } else {
               return TclError(interp,"acq_block error: can not get phase for -ph parameter");
            }
	 }
      } else {
         return TclError(interp,"acq_block error: unknown option %s\n",s);
      }
   }
   
   puls->cacqs[puls->cacqs_pos] = Tcl_DuplicateObj(objv[objc-1]);
   /* this was needed to prevent overwritting of the object with other Tcl stuff */
   Tcl_IncrRefCount(puls->cacqs[puls->cacqs_pos]);
   puls->acqnp[puls->cacqs_pos] = np;
   puls->acqph[(puls->cacqs_pos)++] = phase;
   
   /* test output */
   /*i = puls->cacqs_pos;
   printf("acq_block (%i) parameters summary\n================================\n", i-1);
   printf("Number of points = %d\n",np);
   printf("Phase = %f\n",phase);
   printf("Acq code = \n'%s'\n",Tcl_GetString(puls->cacqs[i-1]));
   */
   return TCL_OK;
}

int tclCheckDwelltime(ClientData data,Tcl_Interp* interp, int objc, Tcl_Obj *objv[])
{
   double dt;
   const double dw=1.0e6*puls->dwellt;

   if (puls->new_gcompute) {
      /*printf("check_dwelltime in gcompute_new does nothing\n");*/
      return TCL_OK;
   }

   dt = puls->t_usec - puls->check_dwelltime_t0;
   /*printf("check_dwelltime report: t0 = %g; t = %g; dt = %g\n",puls->check_dwelltime_t0,puls->t_usec,dt);*/

   if (fabs(dt) < 1e-6) {
      /*printf("check_dwelltime report: t = %g - there was no action\n",puls->t_usec);*/
      return TCL_OK;
   }
   if ( dw-dt > 1e-6 ) {
      /*printf("check_dwelltime report: t = %g - step less than dw (%g < %g) \n",puls->t_usec,dt,dw);*/
      return TCL_OK;
   }
   if ( dt-dw > 1e-6) {
      /*printf("check_dwelltime report: t = %g - step more than dw (%g > %g) \n",puls->t_usec,dt,dw);*/
      /* TclError is not properly displayed, using the hard way */
      /*return TclError(interp,"acq_block error: some events not synchronized with dwell time");*/
      fprintf(stderr,"acq_block error: some events not synchronized with dwell time");
      exit(1);
   }
   
   if (puls->wr < 0.001) {
       /* static case: do nothing, more will be added later */
   } else {
       /* MAS */
      /*printf("check_dwelltime report: t = %g - storing prop to slot %i \n",puls->t_usec,puls->new_direct_mxpos);*/
       _store(puls->new_direct_mxpos);
       (puls->new_direct_mxpos)++;
       puls->check_dwelltime_t0 = puls->t_usec;
      /*printf("check_dwelltime report: reset t0 = %g \n",t0);*/
   }

   return TCL_OK;
}


void tclcmd_new_direct(Tcl_Interp* interp) {

Tcl_CreateObjCommand(interp,"acq_block",(Tcl_ObjCmdProc*)tclAcqBlock,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateObjCommand(interp,"check_dwelltime",(Tcl_ObjCmdProc*)tclCheckDwelltime,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

}
