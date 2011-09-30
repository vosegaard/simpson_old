#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "cm_new.h"
#include "tclutil.h"
#include "cryst.h"
#include "OCroutines.h"
#include "rfshapes.h"
#include "pulse.h"
#include "iodata.h"
#include "tclutil.h"
#include "B0inhom.h"

/* global variable holding all OC parameters */
OCoptPars OCpar;

/* make visible also global variable pulse */
extern Pulse* puls;

/****
 * initialize global variable OCpar
 ****/ 
void OCpar_initialize(void) 
{
  int i;

  /* matrices for OC propagators */
  for (i=0;i<MAXOCPROPS;i++) {
    OCpar.prop[i]=NULL;
    OCpar.dens[i]=NULL;
    OCpar.mx_code[i]=NULL;
  }
  OCpar.mx_pos=0;
  OCpar.gradmode=0;
  OCpar.gradmodeprop=0;
  OCpar.propstatus=0;
  OCpar.var_shapes=NULL;
  OCpar.var_shapes_min=NULL;
  OCpar.var_shapes_max=NULL;
  OCpar.var_shapes_rmsmax=NULL;
  OCpar.grad_shapes=NULL;
  OCpar.isinit=1;
}

/****
 * clear global variable OCpar
 ****/
void OCpar_destroy(void)
{
  int i;
  
  /* matrices for OC propagators */
  for (i=0;i<MAXOCPROPS;i++) {
    if (OCpar.prop[i]) {
      complx_matrix_free(OCpar.prop[i]);
      OCpar.prop[i]=NULL;
    }
    if (OCpar.dens[i]) {
      complx_matrix_free(OCpar.dens[i]);
      OCpar.dens[i]=NULL;
    }
    if (OCpar.mx_code[i]) {
      free( (char*)(OCpar.mx_code[i]) );
      OCpar.mx_code[i]=NULL;
      /* printf("Freeing OCpar.mx_code[%d] \n",i); */
    }
  }
  
  OCpar.mx_pos=0; 
  OCpar.gradmode=0;
  OCpar.gradmodeprop=0;
  
  if (OCpar.var_shapes) {
    free_int_vector(OCpar.var_shapes);
    OCpar.var_shapes=NULL;
  }
  if (OCpar.var_shapes_min) {
    free_double_vector(OCpar.var_shapes_min);
    OCpar.var_shapes_min=NULL;
  }
  if (OCpar.var_shapes_max) {
    free_double_vector(OCpar.var_shapes_max);
    OCpar.var_shapes_max=NULL;
  }
  if (OCpar.var_shapes_rmsmax) {
    free_double_vector(OCpar.var_shapes_rmsmax);
    OCpar.var_shapes_rmsmax=NULL;
  }
  if (OCpar.grad_shapes) {
    free_int_vector(OCpar.grad_shapes);
    OCpar.grad_shapes=NULL;
  }
  OCpar.isinit=0;
}

/****
 * increase matrix pointer
 ****/
void incr_OCmx_pos(void)
{
  if (OCpar.mx_pos+1>MAXOCPROPS) {
     fprintf(stderr,"error: there is no more space for OC propagators\n");
     exit(1);
  }
  
  (OCpar.mx_pos)++;
}

/****
 * store current propagator to OCpar->prop[i]
 ****/
void store_OCprop(void)
{
  int i=OCpar.mx_pos;
  
  /* debug print */
  /* printf("Storing OC propagator to slot %d\n",i); */
   
  if (OCpar.prop[i] == NULL)
     OCpar.prop[i] = complx_matrix_alloc(puls->N,puls->N);
  if (puls->U->col == 1)
     cm_copydiagmv(OCpar.prop[i],puls->U);
  else 
     cmv_copy(OCpar.prop[i],puls->U);

}

/****
 * store current density matrix to OCpar->dens[i]
 ****/
void store_OCdens(void)
{
  int i=OCpar.mx_pos;
  
  /* debug print */
   /*  printf("Storing OC matrix density to slot %d\n",i); */
   
  if (OCpar.dens[i] == NULL)
     OCpar.dens[i]=complx_matrix_alloc(puls->N,puls->N);
  
  cmv_copy(OCpar.dens[i],puls->sigma);

}

/****
 * this serves debugging purposes only
 ****/
void test_print_codes(void)
{
  int i;
  for (i=1; i<=OCpar.mx_pos; i++) {
     if (OCpar.mx_code[i]) {
        printf("/%d/ code is '%s'\n",i,OCpar.mx_code[i]);
     }
  }
}

/****
 * store code to matrix to OCpar.mx_code[i]
 ****/
void set_OCmx_code(char *c)
{
  int i=OCpar.mx_pos;
  
  if (OCpar.mx_code[i] == NULL) {
     OCpar.mx_code[i] = malloc(strlen(c)+1);
     /* printf("Alocating OCpar.mx_code[%d] \n",i); */
  } else {
     free( (char*)(OCpar.mx_code[i]) );
     /* printf("Realocating OCpar.mx_code[%d] \n",i); */
     OCpar.mx_code[i] = malloc(strlen(c)+1);
  }
  strcpy(OCpar.mx_code[i],c);
  /* printf("%3d) -%s-\n",OCpar.mx_pos,OCpar.mx_code[OCpar.mx_pos]); */
  
}

/****
 * store filter matrix as OCprop
 ****/
void store_OCfilter(int num) 
{
  int i=OCpar.mx_pos;
  
  /* debug print */
  /*   printf("Storing OC filter matrix to slot %d\n",i); */
   
  if (OCpar.prop[i] == NULL)
     OCpar.prop[i]=complx_matrix_alloc(puls->N,puls->N);
  
  cmv_copy(OCpar.prop[i],puls->matrix[num]);
}

/****
 * check and set switch for optimization of propagators
 *  (in that mode cummulative propagators are calculated)
 ****/
void test_pulseq_for_acqOC_prop(Tcl_Interp *interp)
{
  char *res, *lin, *dum;
  char buf[256],pulseq[128];
  int l=0;
  
  TclGetString(interp,pulseq,"par","pulse_sequence",0,"pulseq");
  strcpy(buf,"info body ");
  strcat(buf,pulseq);
  if ( Tcl_Eval(interp,buf) != TCL_OK ) {
     fprintf(stderr,"test_pulseq_for_acqOC_prop can not execute %s :\n",buf);
     fprintf(stderr,"%s",interp->result);
     exit(1);
  }
     
  res = Tcl_GetStringResult(interp);

  lin = strtok(res,"\n");
  while ( lin != NULL) {
     l++;
     dum = strstr(lin,"oc_acq_prop");
     if ( dum ) {
        /* found the command, now check if it is commented out */
	printf("found oc_acq_prop on line %d ",l);
	char *comt;
	comt = strchr(lin,'#');
	if ( comt ) {
	   /* comment sign found, check if it is before the command */
	   if ( comt < dum ) {
	      printf("but it is within a comment\n");
	   } else {
	      printf(" and it seems to be active command\n");
	      OCpar.gradmodeprop++;
	   }
        } else {
	   printf(" and it seems to be active command\n");
	   OCpar.gradmodeprop++;
	}
     }  
     lin = strtok(NULL,"\n");
  }
  
  if (OCpar.gradmodeprop > 1) {
     fprintf(stderr,"test_pulseq_for_acqOC_prop error: the command 'oc_acq_prop' can be used\n   only once in the pulse sequence!!! Please do not use reserved word\n   'oc_acq_prop' neither as a text argument for other commands.\n");
     exit(1);
  }

}



/****
 * filter command in gradient mode 
 ****/
void _filterOC(int num)
{
  mv_complx *m;
  int i,j,N;
  complx *sl1, *sl2;

  N = puls->N;
  m=puls->matrix[num];

  if (!m) {
    fprintf(stderr,"error: illegal filter number %d\n",num);
    exit(1);
  }
  if (m->row != N || m->col != N) {
    fprintf(stderr,"error: filter: matrix %d must be a %dx%d matrix\n",num,N,N);
    exit(1);
  }

  if (puls->Uisunit != 1) {
     incr_OCmx_pos();
     store_OCprop();
     _evolve_with_prop();
     _reset_prop();
     store_OCdens();
     set_OCmx_code("P");
  }
  puls->cannotbestored=1;

  sl1 = puls->sigma->data;
  sl2 = m->data;
  for (i=0;i<N*N;i++) {
    if (sl2->re == 0.0 && sl2->im == 0.0) *sl1 = CPLX0;
    sl1++;
    sl2++;
  }

  incr_OCmx_pos();
  store_OCfilter(num);
  store_OCdens();
  set_OCmx_code("F");
}

/****
 * filter command in backevolution when calculating gradients 
 ****/
void lambda_filterOC(int num, mv_complx *lam)
{
  int i, j, N;
  mv_complx *m;
  complx *sl1, *sl2;
  
  N = lam->row;
  m = OCpar.prop[num];
  
  sl1 = lam->data;
  sl2 = m->data;
  for (i=0;i<N*N;i++) {
    if (sl2->re == 0.0 && sl2->im == 0.0) *sl1 = CPLX0;
    sl1++;
    sl2++;
  }
}


/****
 * helper function for pulse_shaped in gradient mode and state to state optimization
 *        (here step by step propagators are saved together with step-evolved
 *         density matrices)
 ****/
 void _pulse_shapedOC(char *code, int Nch, int Nelem, int *mask, double steptime)
{
  int i, j;
  char cd[128];

   if (puls->Uisunit != 1) {
      /* there is some pending propagator, store it but set code for not take gradient */
      incr_OCmx_pos();
      store_OCprop();
      _evolve_with_prop();
      store_OCdens();
      _reset_prop();
      set_OCmx_code("P");
   }
   /* do pulsing and evolving, storing the results */
   for (j=1; j<=Nelem; j++) {
      for (i=1; i<=puls->nchan; i++) {
         if (mask[i] == -1) {
            _rf(i,0.0);
            _ph(i,0.0);
         } else {
	    _rf(i,RFshapes[mask[i]][j].ampl);
	    _ph(i,RFshapes[mask[i]][j].phase);
         }
      }
      _pulse(steptime);
      incr_OCmx_pos();
      store_OCprop();
      _evolve_with_prop();
      store_OCdens();
      _reset_prop();
      if (j==1)
         set_OCmx_code("G");
   }
   /* set the code only for the last memory slot */ 
   sprintf(cd,"G%dE%d",Nelem,Nch); 
   strcat(cd,code);
   set_OCmx_code(cd);

}


/****
 * helper function for pulse_shaped in gradient mode and propagator optimization
 *        (here cummulative propagators are stored and no density matrices)
 ****/
void _pulse_shapedOCprops(char *code, int Nch, int Nelem, int *mask, double steptime)
{
  int i, j;
  char cd[128];

   if ( (puls->Uisunit != 1) && (OCpar.propstatus == 1) ) {
      /* there is some pending propagator, store it but set code for not take gradient */
      incr_OCmx_pos();
      store_OCprop();
      set_OCmx_code("P");
   }
   /* do pulsing and storing cummulative propagators */
   for (j=1; j<=Nelem; j++) {
      for (i=1; i<=puls->nchan; i++) {
         if (mask[i] == -1) {
            _rf(i,0.0);
            _ph(i,0.0);
         } else {
	    _rf(i,RFshapes[mask[i]][j].ampl);
	    _ph(i,RFshapes[mask[i]][j].phase);
         }
      }
      _pulse(steptime);
      incr_OCmx_pos();
      store_OCprop();
      if (j==1) {
         /* set the code only for the first memory slot */ 
         sprintf(cd,"G%dE%d",Nelem,Nch); 
         strcat(cd,code);
         set_OCmx_code(cd);
      }
         
   }
   /* this is just to mark last slot with variable shape propagator*/
   set_OCmx_code("G");
   /* mark current propagator as saved */
   OCpar.propstatus = 0;
}





/****
 * helper function for pulse_and_zgrad_shaped in gradient mode and state to state optimization
 *        (here step by step propagators are saved together with step-evolved
 *         density matrices)
 ****/
 void _pulse_and_zgrad_shapedOC(char *code, int Nch, int Nelem, int *mask, int zgrslot, double steptime)
{
  int i, j, used, Nchan, N;
  char cd[128];
  double *mx, *Hst1;
  double *ovals, zcoor, offs;

  /* prepare for offset term */
  N = puls->N;
  Hst1 = puls->H->Hstatic[1];
  mx = (double*)malloc(N*sizeof(double));
  used = puls->H->nstatic_used[1];
  if (used) {
    /* there already is some offset, remember it */
       memcpy(mx,Hst1,N*sizeof(double));
    } else {
       memset(mx,0,N*sizeof(double));
       memset(Hst1,0,N*sizeof(double));
  }
  Nchan = puls->nchan;
  ovals = double_vector(Nchan);
  zcoor = puls->zcoor;

   if (puls->Uisunit != 1) {
      /* there is some pending propagator, store it but set code for not take gradient */
      incr_OCmx_pos();
      store_OCprop();
      _evolve_with_prop();
      store_OCdens();
      _reset_prop();
      set_OCmx_code("P");
   }
   /* do pulsing and evolving, storing the results */
   for (j=1; j<=Nelem; j++) {
     if (j != 1) memcpy(Hst1,mx,N*sizeof(double));
      get_chan_offset_ratios(puls->ss, 2.0*M_PI*ZgradShapes[zgrslot][j], ovals);
      for (i=1; i<=Nchan; i++) {
         /* finalize offset term */
         offs = zcoor*ovals[i];
         daxpy_(&N,&offs,puls->chan_Iz[i],&INTONE,Hst1,&INTONE);
         /* prepare rf term */
         if (mask[i] == -1) {
            _rf(i,0.0);
            _ph(i,0.0);
         } else {
	    _rf(i,RFshapes[mask[i]][j].ampl);
	    _ph(i,RFshapes[mask[i]][j].phase);
         }
      }
      /* do step offset */
       puls->H->nstatic_used[1] = 1;
      /* do step pulse */
      _pulse(steptime);
      incr_OCmx_pos();
      store_OCprop();
      _evolve_with_prop();
      store_OCdens();
      _reset_prop();
      if (j==1)
         set_OCmx_code("G");
   }
   /* set the code only for the last memory slot */ 
   sprintf(cd,"G%dE%d",Nelem,Nch); 
   strcat(cd,code);
   set_OCmx_code(cd);
   /* clean up after z grad offsets */
   ham_set_offset(puls->H,mx,used);
   free((char*)mx);
   free_double_vector(ovals);  
}


/****
 * helper function for pulse_and_zgrad_shaped in gradient mode and propagator optimization
 *        (here cummulative propagators are stored and no density matrices)
 ****/
void _pulse_and_zgrad_shapedOCprops(char *code, int Nch, int Nelem, int *mask, int zgrslot, double steptime)
{
  int i, j, used, Nchan, N;
  char cd[128];
  double *mx, *Hst1;
  double *ovals, zcoor, offs;

  /* prepare for offset term */
  N = puls->N;
  Hst1 = puls->H->Hstatic[1];
  mx = (double*)malloc(N*sizeof(double));
  used = puls->H->nstatic_used[1];
  if (used) {
    /* there already is some offset, remember it */
       memcpy(mx,Hst1,N*sizeof(double));
    } else {
       memset(mx,0,N*sizeof(double));
       memset(Hst1,0,N*sizeof(double));
  }
  Nchan = puls->nchan;
  ovals = double_vector(Nchan);
  zcoor = puls->zcoor;

   if ( (puls->Uisunit != 1) && (OCpar.propstatus == 1) ) {
      /* there is some pending propagator, store it but set code for not take gradient */
      incr_OCmx_pos();
      store_OCprop();
      set_OCmx_code("P");
   }
   /* do pulsing and storing cummulative propagators */
   for (j=1; j<=Nelem; j++) {
     if (j != 1) memcpy(Hst1,mx,N*sizeof(double));
      get_chan_offset_ratios(puls->ss, 2.0*M_PI*ZgradShapes[zgrslot][j], ovals);
      for (i=1; i<=puls->nchan; i++) {
         /* finalize offset term */
         offs = zcoor*ovals[i];
         daxpy_(&N,&offs,puls->chan_Iz[i],&INTONE,Hst1,&INTONE);
         /* prepare rf term */
         if (mask[i] == -1) {
            _rf(i,0.0);
            _ph(i,0.0);
         } else {
	    _rf(i,RFshapes[mask[i]][j].ampl);
	    _ph(i,RFshapes[mask[i]][j].phase);
         }
      }
      /* do step offset */
      puls->H->nstatic_used[1] = 1;
      /* do step pulse */
      _pulse(steptime);
      incr_OCmx_pos();
      store_OCprop();
      if (j==1) {
         /* set the code only for the first memory slot */ 
         sprintf(cd,"G%dE%d",Nelem,Nch); 
         strcat(cd,code);
         set_OCmx_code(cd);
      }
         
   }
   /* this is just to mark last slot with variable shape propagator*/
   set_OCmx_code("G");
   /* mark current propagator as saved */
   OCpar.propstatus = 0;
   /* clean up after z grad offsets */
   ham_set_offset(puls->H,mx,used);
   free((char*)mx);
   free_double_vector(ovals);  
}



/****
 * function calculating gradients for Hermitian state to state transfer
 ****/
 void gradOC_hermit(void)
{
  int i,j,k,Nelem,Nch,idx,chan,dum,N;
  mv_complx *lam;
  complx cc,gr;
  int *gridx, *chnl, *slt;
  char code, *dumstr, mxcd[128];
  
  N = puls->N;
  /* printf("in gradOC_hermit\n"); */

  /* check if there is pending propagator and save and evolve with it */
  if ( puls->Uisunit != 1) {
      incr_OCmx_pos();
      store_OCprop();
      _evolve_with_prop();
      store_OCdens();
      _reset_prop();
      set_OCmx_code("P");
   }
   
  /* initialize info where to store gradients */
  Nelem = OCpar.grad_shapes[0];
  gridx = int_vector(Nelem);
  dum=0;
  for (j=1;j<=Nelem;j++) {
     dum += RFshapes_len(OCpar.grad_shapes[j]);
     gridx[j]=dum; 
  }
  /* check the size of fid */
  if ( dum > LEN(puls->fid) ) {
    fprintf(stderr,"error: gradient function detected overflow in fid points\n");
    exit(1);
  }
  /* optimisticaly set this and hope all in this function works... */
  puls->curr_nsig = dum;
  
  /* final lambda is */
  lam = cmv_dup(puls->fdetect);
  /* use this as workspace of proper size */
  puls->tmpU->col = puls->N;
  
  i=OCpar.mx_pos;
  while (i>0) {
     strcpy(mxcd,OCpar.mx_code[i]);
     code = mxcd[0];
     /* printf("%d/ code is '%c'\n",i,code); */
     if (code == 'G') {
        /* calculate gradients */
        /* printf("%d is G: %s\n",i,OCpar.mx_code[i]); */
	/* disassemble the code */
	dumstr = strtok(mxcd," ");
	if (sscanf(dumstr,"G%dE%d", &Nelem, &Nch) != 2) {
	   fprintf(stderr,"gradOC_hermit can't decompose first part of code\n");
	   exit(1);
	}
	/* printf("   -> %d elements",Nelem); */
	chnl = int_vector(Nch);
	slt = int_vector(Nch);
	for (j=1; j<=Nch; j++) {
	   dumstr = strtok(NULL, " ");
	   if ( sscanf(dumstr,"I%dC%d",&idx, &chan) != 2 ) {
	      fprintf(stderr,"gradOC_hermit - can't read grad code number %d\n",j);
	      exit(1);
	   }
	   chnl[j] = chan;
	   slt[j] = idx;
	   /* printf(", index %d on channel %d",idx, chan); */
	}
	/* printf("\n"); */
	
	for (j=Nelem; j>0; j--) {
	   /* printf("     doing elem %d\n",i); */
	   if (j == Nelem) {
	      if (i != OCpar.mx_pos) {
                 if (strncmp(OCpar.mx_code[i+1],"F",1) == 0) {
                    /* filter */
                    /* printf("     execute filter with %d\n",i+1); */
		    lambda_filterOC(i+1,lam);
                 } else {
                    /* backevolve */
                    /* printf("     backevolve with %d\n",i+1); */
                    simtrans_adj(lam,OCpar.prop[i+1]);
                 }
              }
	   } else {
	      /* printf("     backevolve with %d\n",i+1); */
              simtrans_adj(lam,OCpar.prop[i+1]);
	   }
	   for (k=1; k<=Nch; k++) {
	      /* printf("     gradient calc. for shape idx %d on channel %d",slt[k],chnl[k]); */
	      zlacrm_(&N,&N,lam->data,&N,puls->chan_Ix[chnl[k]],&N,puls->tmpU->data,&N,(double*)(puls->dU->data));
	      cc = cm_trace(puls->tmpU,OCpar.dens[i]);
	      gr.re = 2.0*cc.im;
              zgemm_("N","N",&N,&N,&N,&CPLX1,lam->data,&N,puls->chan_Iy[chnl[k]],&N,&CPLX0,puls->tmpU->data,&N);
	      cc = cm_trace(puls->tmpU,OCpar.dens[i]);
	      gr.im = 2.0*cc.im;
	      /* printf(", stored in fid[%d]\n",gridx[slt[k]]); */
	      puls->fid[gridx[slt[k]]] = gr;
	      (gridx[slt[k]])--;
	   }
	   i--;
	}
	free_int_vector(chnl);
	free_int_vector(slt);

     } else {
        if (i != OCpar.mx_pos) {
           /* code is either P or F, just evolve or filter */
           /* printf("%d is not G: %s\n",i,OCpar.mx_code[i]); */
           if (strncmp(OCpar.mx_code[i+1],"F",1) == 0) {
              /* filter */
              /* printf("     execute filter with %d\n",i+1); */
	      lambda_filterOC(i+1,lam);
           } else {
              /* backevolve */
              /* printf("     backevolve with %d\n",i+1); */
              simtrans_adj(lam,OCpar.prop[i+1]);
           }
	}
	i--;
     }

  } /*end while */
  
  complx_matrix_free(lam);
  free_int_vector(gridx);
}

/****
 * function calculating gradients for non-Hermitian state to state transfer
 ****/
 void gradOC_nonhermit(void)
{
  int i,j,k,Nelem,Nch,idx,chan,dum,N;
  mv_complx *lam;
  mv_double *wrd;
  complx cc,gr, cc1;
  int *gridx, *chnl, *slt;
  char code, *dumstr, mxcd[128];
  
  /* printf("in gradOC_nonhermit\n"); */
  N = puls->N;
   
  /* check if there is pending propagator and save and evolve with it */
  if ( puls->Uisunit != 1) {
      incr_OCmx_pos();
      store_OCprop();
      _evolve_with_prop();
      store_OCdens();
      _reset_prop();
      set_OCmx_code("P");
   }
   
  /* initialize info where to store gradients */
  Nelem = OCpar.grad_shapes[0];
  gridx = int_vector(Nelem);
  dum=0;
  for (j=1;j<=Nelem;j++) {
     dum += RFshapes_len(OCpar.grad_shapes[j]);
     gridx[j] = dum;
  }
  /* check the size of fid */
  if ( dum > LEN(puls->fid) ) {
    fprintf(stderr,"error: gradient function detected overflow in fid points\n");
    exit(1);
  }
  /* optimisticaly set this and hope all in this function works... */
  puls->curr_nsig = dum;
  
  /* here, lam is holding adjoint of lambda throughout the calculations */
  lam = complx_matrix_alloc(N,N);
  cm_adjoint(lam,puls->fdetect);
  /* be some workspace */
  puls->tmpU->col = N;
  puls->dU->col = N;
  /* prepare constant complex term Tr[rho_N+ lam_N] */
  cc1 = cm_trace_adjoint(OCpar.dens[OCpar.mx_pos],puls->fdetect);
   
  i=OCpar.mx_pos;
  while (i>0) {
     strcpy(mxcd,OCpar.mx_code[i]);
     code = mxcd[0];
     /* printf("%d/ code is '%c'\n",i,code); */
     if (code == 'G') {
        /* calculate gradients */
        /* printf("%d is G: %s\n",i,OCpar.mx_code[i]); */
	/* disassemble the code */
	dumstr = strtok(mxcd," ");
	if (sscanf(dumstr,"G%dE%d", &Nelem, &Nch) != 2) {
	   fprintf(stderr,"gradOC_nonhermit can't decompose first part of code\n");
	   exit(1);
	}
	/* printf("   -> %d elements",Nelem); */
	chnl = int_vector(Nch);
	slt = int_vector(Nch);
	for (j=1; j<=Nch; j++) {
	   dumstr = strtok(NULL, " ");
	   if ( sscanf(dumstr,"I%dC%d",&idx, &chan) != 2 ) {
	      fprintf(stderr,"gradOC_hermit - can't read grad code number %d\n",j);
	      exit(1);
	   }
	   chnl[j] = chan;
	   slt[j] = idx;
	   /* printf(", index %d on channel %d",idx, chan); */
	}
	/* printf("\n"); */
	
	
	/* start loop over elements in shapes */
	wrd = (mv_double*)malloc(sizeof(mv_double));
	wrd->row = wrd->col = N;
	for (j=Nelem; j>0; j--) {
	   /* printf("     doing elem %d\n",i); */
	   if (j == Nelem) {
	      if (i != OCpar.mx_pos) {
                 if (strncmp(OCpar.mx_code[i+1],"F",1) == 0) {
                    /* filter */
                    /* printf("     execute filter with %d\n",i+1); */
		    lambda_filterOC(i+1,lam);
                 } else {
                    /* backevolve */
                    /* printf("     backevolve with %d\n",i+1); */
                    simtrans_adj(lam,OCpar.prop[i+1]);
                 }
              }
	   } else {
	      /* printf("     backevolve with %d\n",i+1); */
              simtrans_adj(lam,OCpar.prop[i+1]);
	   }
	   /* prepare constant matrix [rho_j, lam_j+], note that lam is already adjoint! */
	   cm_commutator(puls->tmpU, OCpar.dens[i], lam);
	   /* loop over channels */
	   for (k=1; k<=Nch; k++) {
	      /* printf("     gradient calc. for shape idx %d on channel %d",slt[k],chnl[k]); */
	      /* cc = Cmul( m_trace(puls->chan_Ix[chnl[k]], puls->tmp), cc1); */
	      wrd->data = puls->chan_Ix[chnl[k]];
	      dmv_complx(wrd,puls->dU);
	      cc = Cmul( cm_trace_adjoint(puls->dU,puls->tmpU), cc1);
	         /* I use this trace since it is faster, Ix is hermitian so adjoint does not change it */
	      gr.re = 2.0*cc.im;
	      /* cc = Cmul( m_trace(puls->chan_Iy[chnl[k]], puls->tmp), cc1); */
	      memcpy(puls->dU->data,puls->chan_Iy[chnl[k]],N*N*sizeof(complx));
	      cc = Cmul( cm_trace_adjoint(puls->dU,puls->tmpU), cc1);
	         /* I use this trace since it is faster, Iy is hermitian so adjoint does not change it */
	      gr.im = 2.0*cc.im;
	      /* printf(", stored in fid[%d]\n",gridx[slt[k]]); */
	      puls->fid[gridx[slt[k]]] = gr;
	      (gridx[slt[k]])--;
	   }
	   i--;
	}
	free_int_vector(chnl);
	free_int_vector(slt);
	free((char*)wrd);

     } else {
        if (i != OCpar.mx_pos) {
           /* code is either P or F, just evolve or filter */
           /* printf("%d is not G: %s\n",i,OCpar.mx_code[i]); */
           if (strncmp(OCpar.mx_code[i+1],"F",1) == 0) {
              /* filter */
              /* printf("     execute filter with %d\n",i+1); */
	      lambda_filterOC(i+1,lam);
           } else {
              /* backevolve */
              /* printf("     backevolve with %d\n",i+1); */
              simtrans_adj(lam,OCpar.prop[i+1]);
           }
	}
	i--;
     }

  } /*end while */
  
  complx_matrix_free(lam);
  free_int_vector(gridx);
}

/****
 * function for calculation gradients for propagator optimization
 ****/
void gradOC_prop(int ud)
{
  int i,j,k,Nelem,Nch,idx,chan,dum,N;
  mv_complx *pr;
  mv_double *wrd;
  complx cc,gr, cc1;
  int *gridx, *chnl, *slt;
  char code, *dumstr, mxcd[128];
  
  /* printf("in gradOC_prop\n"); */
  N = puls->N;
  
  if ( (puls->Uisunit != 1) && (OCpar.propstatus == 1) ) {
     /* there is some pending propagator, store it but set code for not take gradient */
     incr_OCmx_pos();
     store_OCprop();
     set_OCmx_code("P");
  }
 
  /* initialize info where to store gradients */
  Nelem = OCpar.grad_shapes[0];
  gridx = int_vector(Nelem);
  dum=1;
  for (j=1;j<=Nelem;j++) {
     gridx[j]=dum;
     dum += RFshapes_len(OCpar.grad_shapes[j]);
  }
  /* check the size of fid */
  if ( dum-1 > LEN(puls->fid) ) {
    fprintf(stderr,"error: gradient function detected overflow in fid points\n");
    exit(1);
  }
  /* optimisticaly set this and hope all in this function works... */
  puls->curr_nsig = dum-1;

  /* pr is holding U_d+ U_final */
  pr = complx_matrix_alloc(N,N);
  cm_adjoint(pr,puls->STO[ud]);
  cm_multo(pr,OCpar.prop[OCpar.mx_pos]);
  /* cc1 is constant Tr { U_final+ U_d } */
  cc1 = cm_trace_adjoint( OCpar.prop[OCpar.mx_pos], puls->STO[ud]);
  /* workspace */
  puls->tmpU->col = N;
  
  i=1;
  while (i <= OCpar.mx_pos) {
     strcpy(mxcd,OCpar.mx_code[i]);
     code = mxcd[0];
     /* printf("%d/ code is '%c'\n",i,code); */
     if (code == 'G') {
        /* calculate gradients */
        /* printf("%d is G: %s\n",i,OCpar.mx_code[i]); */
	/* disassemble the code */
	dumstr = strtok(mxcd," ");
	if (sscanf(dumstr,"G%dE%d", &Nelem, &Nch) != 2) {
	   fprintf(stderr,"gradOC_nonhermit can't decompose first part of code\n");
	   exit(1);
	}
	/* printf("   -> %d elements",Nelem); */
	chnl = int_vector(Nch);
	slt = int_vector(Nch);
	for (j=1; j<=Nch; j++) {
	   dumstr = strtok(NULL, " ");
	   if ( sscanf(dumstr,"I%dC%d",&idx, &chan) != 2 ) {
	      fprintf(stderr,"gradOC_hermit - can't read grad code number %d\n",j);
	      exit(1);
	   }
	   chnl[j] = chan;
	   slt[j] = idx;
	   /* printf(", index %d on channel %d",idx, chan); */
	}
	/* printf("\n"); */
	
	
	/* start loop over elements in shapes */
	wrd = (mv_double*)malloc(sizeof(mv_double));
	wrd->row = wrd->col = N;
	for (j=1; j<=Nelem; j++) {
	   /* printf("     doing elem %d\n",i); */
	   /* loop over channels */
	   for (k=1; k<=Nch; k++) {
	      /* printf("     gradient calc. for shape idx %d on channel %d",slt[k],chnl[k]); */
	      wrd->data = puls->chan_Ix[chnl[k]];
	      dmv_complx(wrd,puls->tmpU);
	      simtrans_adj(puls->tmpU,OCpar.prop[i]);
	      cc = Cmul( cm_trace(pr, puls->tmpU), cc1);
	      gr.re = 2.0*cc.im;
	      memcpy(puls->tmpU->data,puls->chan_Iy[chnl[k]],N*N*sizeof(complx));
	      simtrans_adj(puls->tmpU,OCpar.prop[i]);
	      cc = Cmul( cm_trace(pr, puls->tmpU), cc1);
	      gr.im = 2.0*cc.im;
	      /* printf(", stored in fid[%d]\n",gridx[slt[k]]); */
	      puls->fid[gridx[slt[k]]] = gr;
	      (gridx[slt[k]])++;
	   }
	   i++;
	}
	free_int_vector(chnl);
	free_int_vector(slt);
	free((char*)wrd);

     } else {
        i++;
     }

  } /*end while */
  
  complx_matrix_free(pr);
  free_int_vector(gridx);

}




/**************************************************************
 * ZT: functions for Optimal Control by Conjugated Gradients  *
 **************************************************************/

/****
 * ZT: function to evaluate target function with current RFshapes
 ****/
 double evaluate_target_function(Tcl_Interp* interp)
{
  double tfval;
  Tcl_Obj* objptr;  

  /* reset grad mode to calculate target function */
  OCpar.gradmode = 0;
  
  if (Tcl_EvalEx( (Tcl_Interp*)interp, "target_function", -1, TCL_EVAL_GLOBAL) != TCL_OK) {
    TclError(interp, "error: Unable to execute command 'target_function' ");
    exit(1);
  }
  objptr = Tcl_GetObjResult( (Tcl_Interp*)interp);
  if ( Tcl_GetDoubleFromObj( (Tcl_Interp*)interp, objptr, &tfval) != TCL_OK) {
    TclError(interp, "error in evaluate_target_function: Unable to convert result into double");
    exit(1);
  }

  return tfval;

}

/****
 * ZT: evaluate gradient procedure with current RFshapes
 ****/
 int evaluate_gradient(Tcl_Interp* interp)
{
  int gr_slot;
  Tcl_Obj* objptr;

  /* set gradient mode for calculating gradients and reset matrix counter */
  OCpar.gradmode = 1;
  OCpar.mx_pos = 0;
  
  if (Tcl_EvalEx((Tcl_Interp*)interp,"gradient",-1,TCL_EVAL_GLOBAL) != TCL_OK) {
    TclError(interp, "error: Unable to execute command 'gradient' ");
    exit(1);
  }
  objptr = Tcl_GetObjResult( (Tcl_Interp*)interp);
  if ( Tcl_GetIntFromObj( (Tcl_Interp*)interp, objptr, &gr_slot) != TCL_OK) {
    TclError(interp, "error in evaluate_gradient: Unable to convert result into integer");
    exit(1);
  }
    
  OCpar.gradmode = 0;
  
  return gr_slot;
}

/****
 * ZT: function to update RFshapes with data moved by a step along dir
 *     IMPORTANT: dir is in x,y format!!!
 ****/
 void update_RFshapes_by_step(double step, double *dir, RFelem *CURRshapes[])
{
  double x, y, am, ph, rms;
  int i,j,ii,Nii;
  int Nshapes = OCpar.var_shapes[0];

  if (step == 0) {
     /* just copy CURRshapes to RFshapes */
     for (i=1; i<=Nshapes; i++) {
        for (ii=1; ii<=RFshapes_len(OCpar.var_shapes[i]); ii++) {
           RFshapes[OCpar.var_shapes[i]][ii].ampl = CURRshapes[i-1][ii].ampl;
	   RFshapes[OCpar.var_shapes[i]][ii].phase = CURRshapes[i-1][ii].phase;
	}
     }
  } else {
     j=0;
     for (i=1; i<=Nshapes; i++) {
        rms = 0.0;
        Nii = RFshapes_len(OCpar.var_shapes[i]);
        for (ii=1; ii<=Nii; ii++) {
	   /* get original x,y components for this shape */
           x = CURRshapes[i-1][ii].ampl*cos(CURRshapes[i-1][ii].phase*DEG2RAD);
           y = CURRshapes[i-1][ii].ampl*sin(CURRshapes[i-1][ii].phase*DEG2RAD);
	   /* move them to new value */
           j++;
	   x += step*dir[j];
	   j++;
	   y += step*dir[j];
	   /* transfer them to ampl, phase */
	   am = sqrt(x*x+y*y);
	   ph = RAD2DEG*atan2(y,x);
           rms += am*am;
	   /* store them in relevant RFshape */
	   RFshapes[OCpar.var_shapes[i]][ii].ampl = am;
	   if (am > OCpar.var_shapes_max[i]) {
	      RFshapes[OCpar.var_shapes[i]][ii].ampl = OCpar.var_shapes_max[i];
	   } 
	   if (am < OCpar.var_shapes_min[i]) {
	      RFshapes[OCpar.var_shapes[i]][ii].ampl = OCpar.var_shapes_min[i];
	   }  
	   RFshapes[OCpar.var_shapes[i]][ii].phase = ph;
        }
        rms = sqrt(rms/((double)Nii));
        if (rms > OCpar.var_shapes_rmsmax[i]) {
           double kkk;
           kkk = OCpar.var_shapes_rmsmax[i]/rms;
           for (ii=1; ii<=Nii; ii++) {
              RFshapes[OCpar.var_shapes[i]][ii].ampl *= kkk;
           }
        }
     }
  } /* end of else */

}


/****
 * ZT: function to minimize during linesearch
 ****/
 double func_to_min(double step, double *dir, RFelem *CURRshapes[], Tcl_Interp *interp)
{
  double fval;
  
  /* update RFshapes. Assumes RFshapes in ampl/phase and gradients in x/y */
  update_RFshapes_by_step(step,dir,CURRshapes);
  /* call target function */
  fval = -evaluate_target_function(interp);
  
  return fval;
}

/****
 * ZT: bracketing minimum with 3 values
 ****/
 int bracketminCG(double *a, double *b, double *c, double *fa, double *fb, double *fc, double *dir, RFelem *CURRshapes[], Tcl_Interp *interp)
{
   double golden = 0.381966;
   double EPS = 1.0e-10 ; /*insignificant difference between two numbers*/
   double f_l, f_r, f_m, x_l, x_r, x_m, p, q, u, umax, v;
   int Neval = 0;
   int status=0;

   x_l = *a;
   x_r = *b;
   f_l = *fa; /* passed from caller, save one evaluation! */
   f_r = func_to_min(x_r, dir, CURRshapes, interp);
   Neval++;
   
   if (f_r >= f_l) {
      /* go into interval x_l,x_r */
      x_m = (x_r -x_l)*golden + x_l;
      f_m = func_to_min(x_m, dir, CURRshapes, interp);
      Neval++;
   } else {
      /* go further to the right */
      x_m = x_r;
      f_m = f_r;
      x_r = (x_m - x_l)/golden + x_l;
      f_r = func_to_min(x_r, dir, CURRshapes, interp);
      Neval++;
   }
   
   while ( (Neval<OCpar.max_brack_eval) && (fabs(x_l-x_r)>EPS) ) {
      if (f_m < f_l) {
         if (f_m < f_r) {
	    *a = x_l;
	    *b = x_m;
	    *c = x_r;
	    *fa = f_l;
	    *fb = f_m;
	    *fc = f_r;
	    if ( OCpar.verb ) {
               printf("bracketing reached after %d iterations\n",Neval);
            }
	    return status;
	 } else {
	    /* golden expansion limit */
	    v = (x_r-x_m)/golden+x_m;
	    /* use parabolic function */
	    p = (x_m-x_l)*(f_m-f_r);
	    q = (x_m-x_r)*(f_m-f_l);
            u = ( (x_m+x_r)*q - (x_m+x_l)*p )/(q-p)/2;
	    if ( OCpar.verb ) { printf("-1-"); }
	    if ( u > v ) {
	       /* accept parabolic fit only when it goes further than golden expansion */
	    /* if ( (u-x_r) > EPS ) {  */
	       /* u is further on the right = GOOD, but don't go too far */
	       umax = x_l+100.0*(x_r-x_l);
	       x_l = x_m;
               f_l = f_m;
               x_m = x_r;
               f_m = f_r;
               x_r = (u < umax) ? u : umax;
	       f_r = func_to_min(x_r, dir, CURRshapes, interp);
	       Neval++;
	       if ( OCpar.verb ) { printf("-2-"); }
	    } else if ( ((u-x_m)>EPS) && ((x_r-u)>EPS) ) {
	       /* u is between x_m, x_r = GOOD */
	       x_l = x_m;
	       f_l = f_m;
	       x_m = u;
	       f_m = func_to_min(x_m, dir, CURRshapes, interp);
	       Neval++;
	       if ( OCpar.verb ) { printf("-3-"); }
	    } else {
	       /* parabolic function was useless, take golden expansion */
	       x_l = x_m;
	       f_l = f_m;
	       x_m = x_r;
	       f_m = f_r;
	       x_r = (x_m-x_l)/golden+x_l; 
	       f_r = func_to_min(x_r, dir, CURRshapes, interp);
	       Neval++;
	       if ( OCpar.verb ) { printf("-4-"); }
	    }
	 }
      } else {
         if ( f_m < f_r ) {
	    /* try parabola fist */
	    p = (x_m-x_l)*(f_m-f_r);
	    q = (x_m-x_r)*(f_m-f_l);
            u = ( (x_m+x_r)*q - (x_m+x_l)*p )/(q-p)/2;
	    if ( OCpar.verb ) { printf("-5-"); }
	    if ( ((u-x_l)>EPS) && ((x_m-u)>EPS) ) {
	       /* u is between x_l, x_m = GOOD */
	       x_r = x_m;
	       f_r = f_m;
	       x_m = u;
	       f_m = func_to_min(x_m, dir, CURRshapes, interp);
	       Neval++;
	       if ( OCpar.verb ) { printf("-6-"); }
	    } else {
	       /* parabola was vaste of time, but min must be somewhere between x_l, x_m */
	       x_r = x_m;
	       f_r = f_m;
	       x_m =(x_r-x_l)*golden+x_l;
	       f_m = func_to_min(x_m, dir, CURRshapes, interp);
	       Neval++;
	       if ( OCpar.verb ) { printf("-7-"); }
	    }
	 } else {
	    /* also here we insist that min is right to x_l and step to right was too big */
            x_r = x_m;
            f_r = f_m;
            x_m =(x_r-x_l)*golden+x_l;
            f_m = func_to_min(x_m, dir, CURRshapes, interp);
            Neval++;
	    if ( OCpar.verb ) { printf("-8-"); }
	 }
      }
      
      if ( OCpar.verb ) { 
         printf("bracketing status:x=(%f, %f, %f), F=(%.15g, %.15g, %.15g)\n",x_l,x_m,x_r,f_l,f_m,f_r);   
      }

   } /* end of while */
   
   *a = x_l;
   *b = x_m;
   *c = x_r;
   *fa = f_l;
   *fb = f_m;
   *fc = f_r;
   printf("Bracketmin: Failed to bracket, maximum function evaluation reached or limits closer than machine precision\n");
   printf("ERROR!!! --> Aborting optimization...\n");
   /* exit(1); */
   status = 1;
   return status;
   
}

/****
 * ZT: find minimum with Brent method
 ****/
 void brentCG(double *a, double *b, double *c, double *fa, double *fb, double *fc, double *dir, RFelem *CURRshapes[], Tcl_Interp *interp)
{
   double golden = 0.3819660;
   double EPS = 1.4901161193847656e-08; /* sqrt(double precision) */
   double f_l, f_r, f_m, x_l, x_r, x_m, v, w, d, e, f_v, f_w, z, f_z;
   double w_lo, w_up, tolerance, p, q, r, midpoint, u, f_u;
   double stoptol;
   int Neval = 0;
   
   x_l = *a;
   f_l = *fa;
   x_m = *b;
   f_m = *fb;
   x_r = *c;
   f_r = *fc;
   
   /* initialization */
   stoptol = ( x_l*x_r > 0.0 ) ? (OCpar.tol*x_m) : OCpar.stepmin; 
   v = x_l + golden*(x_r - x_l);
   w = v;
   d = 0.0;
   e = 0.0;
   f_v = func_to_min(v, dir, CURRshapes, interp);
   Neval++;
   f_w = f_v;
   
   while ( (Neval<OCpar.max_brent_eval) && (fabs(x_r-x_l)>stoptol) ) {
    /* printf("Brent status1: x=<%f, %f>, f=<%f, %f>\n",x_l,x_r,f_l,f_r); */
      z = x_m;
      f_z = f_m;
      d = e;
      e = d;
      w_lo = z -x_l;
      w_up = x_r -z;
      tolerance = EPS*abs(z);
      p = 0.0;
      q = 0.0;
      r = 0.0;
      midpoint = 0.5*(x_l+x_r);
      
      if ( fabs(e)>tolerance) {
         /* fit parabola */
	 r = (z-w)*(f_z-f_v);
	 q = (z-v)*(f_z-f_w);
	 p = (z-v)*q - (z-w)*r;
	 q = 2*(q-r);
	 if (q>0) {
	    p = -p;
	 } else {
	    q = -q;
	 }
	 r = e;
	 e = d;
      }
      
      if  ( (fabs(p)<fabs(0.5*q*r)) && (p<q*w_lo) && (p<q*w_up) ) {
         double t2 = 2.0*tolerance;
	 d = p/q;
	 u = z + d;
	 if ( (u-x_l<t2) || (x_r-u<t2) )
	    d = (z<midpoint) ? tolerance : -tolerance;
      } else {
         e = ( (z<midpoint) ? x_r : x_l ) - z;
	 d = golden*e;
      }
      
      if ( fabs(d) >= tolerance ) {
         u = z + d;
      } else {
         u = ( (d>0.0) ? tolerance : -tolerance) + z;
      }
      
      f_u = func_to_min(u, dir, CURRshapes, interp);
      Neval++;
      
      if ( f_u<f_z ) {
         if ( u<z) {
	    x_r = z;
	    f_r = f_z;
	 } else {
	    x_l = z;
	    f_l = f_z;
	 }
	 v = w;
	 f_v = f_w;
	 w = z;
	 f_w = f_z;
	 x_m = u;
	 f_m = f_u;
      } else {
         if ( u<z) {
	    x_l = u;
	    f_l = f_u;
	 } else {
	    x_r = u;
	    f_r = f_u;
	 }
	 if ( (f_u <= f_w) || (w == z) ) {
	    v = w;
	    f_v = f_w;
	    w = u;
	    f_w = f_u;
	 } else if ( (f_u<=f_v) || (v==z) || (v==w) ) {
	    v = u;
	    f_v = f_u;
	 }
      }
      if ( OCpar.verb ) { 
         printf("Brent status2: x=<%f, %f>, f=<%f, %f>\n",x_l,x_r,f_l,f_r);
      }
      stoptol = ( x_l*x_r > 0.0 ) ? (OCpar.tol*x_m) : OCpar.stepmin; 
   } /* end of 'while' iteration loop */
   
   if ( OCpar.verb ) {
      printf("Brent used %d iterations\n",Neval);
   }
   *b = x_m;
   *fb = f_m;
   /* other brackets are not changed on exit (since they are not relevant) */
   /* but I added print-out so I rather change them, it looks nicer */
   *a = x_l;
   *fa = f_l;
   *c = x_r;
   *fc = f_r;

}

/****
 * ZT: linesearch envelope
 ****/
 double linesearchCG(double* fcurr, double *dir, RFelem *CURRshapes[], Tcl_Interp* interp)
{
  double stepsize, a, b, c, fa, fb, fc;
  int bkstat;
  
  /* find triplet enclosing minimum */
  a = 0.0;                /* we start here with bracketing */
  b = (double)OCpar.mnbkstep; /* this defines initial step size for bracketing */
  c = 2.0*b;                /* this is not realy used on call, just on exit */
  fa = *fcurr;              /* passed from caller, save one evaluation! */
  bkstat = bracketminCG(&a,&b,&c,&fa,&fb,&fc,dir,CURRshapes,interp);
  if ( OCpar.verb ) {
     printf("BRACKETING minimum results:\n %10.5f %10.5f\n %10.5f %10.5f\n %10.5f %10.5f\n",a,fa,b,fb,c,fc);
  }
  if ( bkstat == 0) {
     /* find local minimum using Brent method */
     brentCG(&a,&b,&c,&fa,&fb,&fc,dir,CURRshapes,interp);
     if ( OCpar.verb ) {
        printf("BRENT minimization results:\n %10.5f %10.5f\n %10.5f %10.5f\n %10.5f %10.5f\n",a,fa,b,fb,c,fc);
     }
     stepsize = b;
     *fcurr = fb;
  } else {
     stepsize = 0;
     /* *fcurr is not changed when bracketing failed... */
  }

  
  /******** D E B U G G I N G ****** R E M O V E    T H I S ******/
  /*    double g1=-500.0, g2=1000.0, gs=2.0;
      double dum;
      printf("----PROFILE START-----\n");
      while (g1<=g2) {
          dum = func_to_min(g1,dir,CURRshapes,interp);
          printf("   %.15g   %.15g\n",g1, dum);
	  fflush(stdout);
          g1 += gs; 
      }
      printf("----PROFILE END------\n");
      fflush(stdout);                                            */
  /******** D E B U G G I N G ****** R E M O V E    T H I S ******/
  
  
  /* make sure that current RFshapes hold this new position */
  update_RFshapes_by_step(stepsize,dir,CURRshapes);

  return stepsize;
}

/****
 * ZT: function to store RFshapes, option enabling monitoring of changes to shapes
 *     during optimization
 ****/
 void store_OCshapes(Tcl_Interp *interp)
{
  /* not used when 'none' was specified */
  if ( strcmp(OCpar.VarSaveProc,"none") ) {
     if (Tcl_EvalEx( (Tcl_Interp*)interp, OCpar.VarSaveProc, -1, TCL_EVAL_GLOBAL) != TCL_OK) {
        TclError(interp, "error: Unable to execute command '%s' for storing RF shapes during oc_optimize\n",OCpar.VarSaveProc);
        exit(1);
     }
  }
  
}


/****
 * helper function for freeing fid from memory (when gradients are no longer needed)
 ****/
void free_grad_fid(int gr_slot, Tcl_Interp *interp)
{ 
  char buf1[64];
  extern int funload(int fidN);
  
  if (funload(gr_slot) == -1) {
    fprintf(stderr,"Error in free_grad_fid\n");
    exit(1);
  }
  sprintf(buf1,"%d",gr_slot);
  Tcl_UnsetVar2(interp,"FD_Internal::f",buf1,0);
}

/****
 * ZT: conjugated gradients
 ****/
 double OptimizeCG(Tcl_Interp* interp)
{
 double fold, fcurr, dum, g, beta;
 double step, tfval;
 double *r, *rN, *dir;
 int gr_slot, count, iter, i, j, N;
 int Nshapes = OCpar.var_shapes[0];
 int Nvars=0;
 RFelem *CURRshapes[Nshapes];
 extern FD** fd;
 
 
 /* initialize function value and gradient */
 tfval = evaluate_target_function(interp); /* target function should be maximzed, here we minimize */
 fold = -tfval;
 printf("     Initial target function: %.10f \n", tfval);
 gr_slot = evaluate_gradient(interp);

 /* store current RF shapes at safe place */
 for (i=0; i<Nshapes; i++) {
    N = RFshapes_len(OCpar.var_shapes[i+1]);
    Nvars += N;
    CURRshapes[i]=RFshapes_alloc(N);
    for (j=1; j<=N; j++) {
       CURRshapes[i][j].ampl = RFshapes[OCpar.var_shapes[i+1]][j].ampl;
       CURRshapes[i][j].phase = RFshapes[OCpar.var_shapes[i+1]][j].phase;
    }
 }
 
 /* check length of gradients with total number of variables */
 if (Nvars != (fd[gr_slot]->np)) {
    fprintf(stderr,"oc_optimize error: number of gradients does not match total number of variables in RF shapes (%d versus %d)\n",Nvars,(fd[gr_slot]->np) );
    exit(1);
  }
 /* Nvars so far counts just elements in RF shapes. But vector dimension is double (aml, ph) */
 Nvars *=2;
 
 /* initialize variables for conjugated gradients */
 dir = double_vector(Nvars); 
 r = double_vector(Nvars);
 rN = double_vector(Nvars);
 i = 0;
 for (j=1; j<=fd[gr_slot]->np; j++) {
     i++;
     dum = (fd[gr_slot]->data)[j].re;
     r[i] = dum; /* gradient calculated for maximizing, here we minimize. r[i]=-dum */
     dir[i] = r[i];
     i++;
     dum = (fd[gr_slot]->data)[j].im;
     r[i] = dum; /* gradient calculated for maximizing, here we minimize. r[i]=-dum */
     dir[i] = r[i];
 }
 /* here all work with gadients in FID is done, free it */
 free_grad_fid(gr_slot,interp);
 
 count = 0;
 fcurr = fold;
 
 /* start iteration loop */
 for (iter=1; iter<=OCpar.nIterations; iter++) {
    /* do linesearch */
    step = linesearchCG(&fcurr, dir, CURRshapes, interp);
    tfval = -fcurr;
    /* store new point RFshapes to safe place (update CURRshapes) */
     for (i=1; i<=Nshapes; i++) {
        for (j=1; j<=RFshapes_len(OCpar.var_shapes[i]); j++) {
           CURRshapes[i-1][j].ampl = RFshapes[OCpar.var_shapes[i]][j].ampl;
	   CURRshapes[i-1][j].phase = RFshapes[OCpar.var_shapes[i]][j].phase;
	}
     }
     count++;
    
    /* print out results */
    printf("    Iter %d: tf=%.10f (step=%g)\n",iter, tfval, step);
    
    /* store sequences every nreport steps */
    if ( (iter % OCpar.nreport) < 1 ) {
       printf("   Storing data for iteration %d\n",iter);
       store_OCshapes(interp);
    }
    
    /* check for done */
    if (fabs(fcurr-fold)<=OCpar.eps) {
       printf("Done by reaching target function tolerance limit.\n\n");
       break;
    }
    if( (tfval<OCpar.cut) && (iter>OCpar.ncut) ) {
       printf("Done %d iterations, but target function below cut-off limit.\n\n",iter);
       break;
    }
    if ( fabs(step)<OCpar.stepmin ) {
       printf("Done by reaching minimal step in line-search.\n\n");
       break;
    }              
    
    /* calculate new gradient */
    gr_slot = evaluate_gradient(interp);
    
    /* prepare parameters for conjugation */
    g = 0.0;
    beta = 0.0;
    i = 0;
    for (j=1; j<=fd[gr_slot]->np; j++) {
       i++;
       dum = (fd[gr_slot]->data)[j].re;
       rN[i] = dum; /* gradient calculated for maximizing, here we minimize. rN[i]=-dum */
       g += r[i]*r[i];
       beta += (rN[i]-r[i])*rN[i]; /* this is for Polak-Ribiere */
       i++;
       dum = (fd[gr_slot]->data)[j].im;
       rN[i] = dum; /* gradient calculated for maximizing, here we minimize. rN[i]=-dum */
       g += r[i]*r[i];
       beta += (rN[i]-r[i])*rN[i]; /* this is for Polak-Ribiere */
    }
    /* here all work with gadients in FID is done, free it */
    free_grad_fid(gr_slot,interp);
    
    /* check for exact reaching of extremum (nul-gradient) */
    if (g==0.0) {
       printf("Done by reaching extremum (nul gradient).\n\n");
       break;
    }
    beta = beta/g;
    /* reset search direction if beta<0 or all independent directions were already taken */
    if ( (beta < 0.0) || (count == Nvars) ) {
       printf("(reseting search directions)\n");
       beta = 0.0;
       count = 0;
    }
    
    /* generate new direction by conjugation and prepare for next iteration */
    for (i=1; i<=Nvars; i++) {
       dir[i] = rN[i]+beta*dir[i]; 
       r[i]=rN[i];
    }
    fold = fcurr;
        
 } /* end of iteration loop */
 
 if (iter >= OCpar.nIterations) 
    printf("Done by reaching maximal number of iterations.\n\n");
 
 
 free_double_vector(dir);
 free_double_vector(r);
 free_double_vector(rN);
 for (i=0; i<Nshapes; i++) {
    free((char *)CURRshapes[i]);
 }
 /* make sure that correct final target function value is returned */
 tfval = -fcurr;
 
 return tfval;
 
}


/****
 * ZT: implementation of oc_optimize
 ****/
 int tclOCoptimize2(ClientData data,Tcl_Interp* interp, int argc, char *argv[])
{
  double tfval, dumfloat;
  Tcl_Obj* objptr;
  int i, slot, ndim, nsh, ish;
  char rxbuf[32];

  /* disable when relaxation is ON */
  TclGetString(interp,rxbuf,"par","relax",0,"off");
  printf("OC detects: relax is %s\n",rxbuf);
  if (!strncmp(rxbuf,"on",2)) {
     fprintf(stderr,"oc_optimize error: relaxation not supported during optimization.\n");
     exit(1);
  }
  
  /* read input parameters */
  if ( argc == 1)
    return TclError(interp,"Usage: <double> oc_optimize <RFshape> ?-min XXX -max XXX -rmsmax XXX? ?<RFshape> ...?");
   
  /* go through argv and count shapes */
  nsh=0;
  for (i=1; i<argc; i++) {
     if (strncmp(argv[i],"-",1)) {
        /* this seems to be a shape */
        nsh++;
     } else {
        /* it is a switch, skip over its value */
        i++;
     }
  }
  if (nsh == 0)
     return TclError(interp,"oc_optimize has not detected any variable shape among arguments!");
     
  ndim = 0;
  ish = 0;
  OCpar_initialize();
  OCpar.var_shapes = int_vector(nsh); 
  OCpar.var_shapes_min = double_vector(nsh); 
  OCpar.var_shapes_max = double_vector(nsh); 
  OCpar.var_shapes_rmsmax = double_vector(nsh);
  OCpar.grad_shapes = int_vector(nsh);
  for (i=1; i<argc; i++) {
     if (strncmp(argv[i],"-",1)) {
        /* this seems to be a shape */
        if (Tcl_GetInt(interp,argv[i],&slot) != TCL_OK) {
           OCpar_destroy();
           return TclError(interp,"oc_optimize: argument %d must be integer <RFshape>",i);
           }
        if (!RFshapes[slot]) {
           OCpar_destroy();
           return TclError(interp,"oc_optimize: argument %d -> RFshape does not exist", i);
        }
        ish++;
        OCpar.var_shapes[ish] = slot;
        OCpar.var_shapes_min[ish] = 0.0;
        OCpar.var_shapes_max[ish] = 1.0e6;
        OCpar.var_shapes_rmsmax[ish] = 1.0e6;
        OCpar.grad_shapes[ish] = slot;
        ndim += RFshapes_len(slot);
     } else {
        /*  it seems to be limits */
        if (ish == 0) {
           OCpar_destroy();
           return TclError(interp,"oc_optimize error reading arguments: limits preceed shape!");
        }
        if ( !strcmp(argv[i],"-min") ) {
           if ( Tcl_GetDouble( (Tcl_Interp*)interp,argv[i+1],&dumfloat) != TCL_OK ) {
              OCpar_destroy();
              return TclError(interp,"oc_optimize error: %d. shape, switch -min: can not get double from %s",ish,argv[i+1]);
           } 
           OCpar.var_shapes_min[ish] = dumfloat;
           i++;
        } else if ( !strcmp(argv[i],"-max") ) {
           if ( Tcl_GetDouble( (Tcl_Interp*)interp,argv[i+1],&dumfloat) != TCL_OK ) {
              OCpar_destroy();
              return TclError(interp,"oc_optimize error: %d. shape, switch -max: can not get double from %s",ish,argv[i+1]);
           } 
           OCpar.var_shapes_max[ish] = dumfloat;
           i++;
        } else if ( !strcmp(argv[i],"-rmsmax") ) {
           if ( Tcl_GetDouble( (Tcl_Interp*)interp,argv[i+1],&dumfloat) != TCL_OK ) {
              OCpar_destroy();
              return TclError(interp,"oc_optimize error: %d. shape, switch -rmsmax: can not get double from %s",ish,argv[i+1]);
           } 
           OCpar.var_shapes_rmsmax[ish] = dumfloat;
           i++;
        } else {
           OCpar_destroy();
           return TclError(interp,"oc_optimize error reading %d.argument '%s' ",i,argv[i]);
        }
     }
  }
  /* debug output */
  for (i=1; i<=nsh; i++) {
     printf("Argument %d: shape slot %d, min=%g, max=%g, RMSmax=%g\n",i,OCpar.var_shapes[i],OCpar.var_shapes_min[i], OCpar.var_shapes_max[i], OCpar.var_shapes_rmsmax[i]);
  }
  /* reading of input parameters is complete */
 
  test_pulseq_for_acqOC_prop(interp);
  OCpar.ndim = ndim;

  /* read in all parameters from par(OC_?) or set their default values */
  objptr = Tcl_GetVar2Ex( (Tcl_Interp*)interp, "par", "oc_tol_cg", TCL_GLOBAL_ONLY);
  if (!objptr) {
     OCpar.eps = 1.0e-6;
  } else {
     if ( Tcl_GetDoubleFromObj( (Tcl_Interp*)interp,objptr,&OCpar.eps) != TCL_OK ) {
       OCpar.eps = 1.0e-6;
       printf("Warning! Could'n read par(oc_tol_cg), using default value %f\n",OCpar.eps);
     } 
  }
  objptr = Tcl_GetVar2Ex( (Tcl_Interp*)interp, "par", "oc_tol_ls", TCL_GLOBAL_ONLY);
  if (!objptr) {
     OCpar.tol = 1.0e-3;
  } else {
     if ( Tcl_GetDoubleFromObj( (Tcl_Interp*)interp,objptr,&OCpar.tol) != TCL_OK ) {
       OCpar.tol = 1.0e-3;
       printf("Warning! Could'n read par(oc_tol_ls), using default value %f\n",OCpar.tol);
     } 
  }
  objptr = Tcl_GetVar2Ex( (Tcl_Interp*)interp, "par", "oc_mnbrak_step", TCL_GLOBAL_ONLY);
  if (!objptr) {
     OCpar.mnbkstep = 10.0;
  } else {
     if ( Tcl_GetDoubleFromObj( (Tcl_Interp*)interp,objptr,&OCpar.mnbkstep) != TCL_OK ) {
       OCpar.mnbkstep = 10.0;
       printf("Warning! Could'n read par(oc_mnbrak_step), using default value %f\n",OCpar.mnbkstep);
     } 
  }
  objptr = Tcl_GetVar2Ex( (Tcl_Interp*)interp, "par", "oc_verbose", TCL_GLOBAL_ONLY);
  if (!objptr) {
     OCpar.verb = 0;
  } else {
     if ( Tcl_GetIntFromObj( (Tcl_Interp*)interp,objptr,&OCpar.verb) != TCL_OK ) {
       OCpar.verb = 0;
       printf("Warning! Could'n read par(oc_verbose), using default value %d\n",OCpar.verb);
     } 
  }
  objptr = Tcl_GetVar2Ex( (Tcl_Interp*)interp, "par", "oc_max_iter", TCL_GLOBAL_ONLY);
  if (!objptr) {
     OCpar.nIterations = 1000;
  } else {
     if ( Tcl_GetIntFromObj( (Tcl_Interp*)interp,objptr,&OCpar.nIterations) != TCL_OK ) {
       OCpar.nIterations = 1000;
       printf("Warning! Could'n read par(oc_max_iter), using default value %d\n",OCpar.nIterations);
     } 
  }
  objptr = Tcl_GetVar2Ex( (Tcl_Interp*)interp, "par", "oc_cutoff", TCL_GLOBAL_ONLY);
  if (!objptr) {
     OCpar.cut = 0.0;
  } else {
     if ( Tcl_GetDoubleFromObj( (Tcl_Interp*)interp,objptr,&OCpar.cut) != TCL_OK ) {
       OCpar.cut = 0.0;
       printf("Warning! Could'n read par(oc_cutoff), using default value %f\n",OCpar.cut);
     } 
  }
  objptr = Tcl_GetVar2Ex( (Tcl_Interp*)interp, "par", "oc_cutoff_iter", TCL_GLOBAL_ONLY);
  if (!objptr) {
     OCpar.ncut = 1000;
  } else {
     if ( Tcl_GetIntFromObj( (Tcl_Interp*)interp,objptr,&OCpar.ncut) != TCL_OK ) {
       OCpar.ncut = 1000;
       printf("Warning! Could'n read par(oc_cutoff_iter), using default value %d\n",OCpar.ncut);
     } 
  }
  objptr = Tcl_GetVar2Ex( (Tcl_Interp*)interp, "par", "oc_var_save_iter", TCL_GLOBAL_ONLY);
  if (!objptr) {
     OCpar.nreport = 1001;
  } else {
     if ( Tcl_GetIntFromObj( (Tcl_Interp*)interp,objptr,&OCpar.nreport) != TCL_OK ) {
       OCpar.nreport = 1001;
       printf("Warning! Could'n read par(oc_var_save_iter), using default value %d\n",OCpar.nreport);
     } 
  } 
  objptr = Tcl_GetVar2Ex( (Tcl_Interp*)interp, "par", "oc_cg_min_step", TCL_GLOBAL_ONLY);
  if (!objptr) {
     OCpar.stepmin = 1.0e-3;
  } else {
     if ( Tcl_GetDoubleFromObj( (Tcl_Interp*)interp,objptr,&OCpar.stepmin) != TCL_OK ) {
       OCpar.stepmin = 1.0e-3;
       printf("Warning! Could'n read par(oc_cg_min_step), using default value %f\n",OCpar.stepmin);
     } 
  }
 
  objptr = Tcl_GetVar2Ex( (Tcl_Interp*)interp, "par", "oc_var_save_proc", TCL_GLOBAL_ONLY);
  if (!objptr) {
     strcpy(OCpar.VarSaveProc,"none");
     objptr = Tcl_NewStringObj(OCpar.VarSaveProc,strlen(OCpar.VarSaveProc));
     objptr = Tcl_SetVar2Ex( (Tcl_Interp*)interp,"par", "oc_var_save_proc", objptr, TCL_GLOBAL_ONLY);
  } else {
     strcpy(OCpar.VarSaveProc, Tcl_GetString(objptr));
     if ( OCpar.VarSaveProc  == NULL ) {
       strcpy(OCpar.VarSaveProc,"none");
       objptr = Tcl_NewStringObj(OCpar.VarSaveProc,strlen(OCpar.VarSaveProc));
       objptr = Tcl_SetVar2Ex( (Tcl_Interp*)interp,"par", "oc_var_save_proc", objptr, TCL_GLOBAL_ONLY);
       printf("Warning! Could'n read par(oc_var_save_proc), using default value %s\n",OCpar.VarSaveProc);
     } 
  }
  objptr = Tcl_GetVar2Ex( (Tcl_Interp*)interp, "par", "oc_max_brack_eval", TCL_GLOBAL_ONLY);
  if (!objptr) {
     OCpar.max_brack_eval = 100;
  } else {
     if ( Tcl_GetIntFromObj( (Tcl_Interp*)interp,objptr,&OCpar.max_brack_eval) != TCL_OK ) {
       OCpar.max_brack_eval = 100;
       printf("Warning! Could'n read par(oc_max_brack_eval), using default value %d\n",OCpar.max_brack_eval);
     } 
  }
  objptr = Tcl_GetVar2Ex( (Tcl_Interp*)interp, "par", "oc_max_brent_eval", TCL_GLOBAL_ONLY);
  if (!objptr) {
     OCpar.max_brent_eval = 100;
  } else {
     if ( Tcl_GetIntFromObj( (Tcl_Interp*)interp,objptr,&OCpar.max_brent_eval) != TCL_OK ) {
       OCpar.max_brent_eval = 100;
       printf("Warning! Could'n read par(oc_max_brent_eval), using default value %d\n",OCpar.max_brent_eval);
     } 
  }
 
 
  printf("Number of variables is %d\n",OCpar.ndim);   
  printf("Global tolerance on target function is %g\n",OCpar.eps);
  printf("Maximal number of iterations is %d\n",OCpar.nIterations);
  printf("Optimization is terminated after %d iterations if target function is not higher than %g\n",OCpar.ncut,OCpar.cut);
  printf("Every %d iterations procedure '%s' (for storing variables) is executed\n",OCpar.nreport, OCpar.VarSaveProc);
  printf("Minimal step size along conjugated gradient is %g\n",OCpar.stepmin);
  printf("Tolerance for line-search (Brent) is %g\n",OCpar.tol);
  printf("Line-search will terminate after reaching %d evaluations\n",OCpar.max_brent_eval);
  printf("Initial step for bracketing minimum is %g\n",OCpar.mnbkstep);
  printf("Bracketing will fail after conducting %d evaluations\n\n",OCpar.max_brack_eval);
 
 
  tfval = OptimizeCG(interp);
 
  sprintf(interp->result,"%lf",tfval); 
  
 
  OCpar_destroy();
 
  return TCL_OK;
} 





/****
 * grad_shapes is used to define RFshapes for gradient calculation with respect to them
 *   can be called only from inside of 'gradient' procedure (here is just simple
 *    test on existence of OCpar.gradmode)
 ****/
 int tclGradShapes(ClientData data,Tcl_Interp* interp, int argc, char *argv[])
{
  int *shapes;
  int i, slot;

  /* test for existence of gradmode */
  if (!OCpar.gradmode) 
     return TclError(interp,"ERROR: gradient mode not activated! oc_grad_shapes can be called only within 'gradient' procedure");
     
  /* read input parameters (only RFshapes) */
  if ( argc == 1)
    return TclError(interp,"Usage: oc_grad_shapes <RFshape> ?<RFshape> ...?");
    
  shapes = int_vector(argc-1);
  for (i=1; i<argc; i++) {
    if (Tcl_GetInt(interp,argv[i],&slot) != TCL_OK) 
      return TclError(interp,"oc_grad_shapes: argument %d must be integer <RFshape>",i);
    if (!RFshapes[slot])
      return TclError(interp,"oc_grad_shapes: argument %d -> RFshape does not exist", i);
    shapes[i] = slot;
  }
  
  i = shapes[0];
  if (OCpar.grad_shapes != NULL) {
     /* printf("oc_grad_shapes is realocating OCpar.grad_shapes\n"); */
     free_int_vector(OCpar.grad_shapes);
  } 
  /* OCpar.grad_shapes = int_vector(i);
  for (slot=1; slot<=i; slot++) {
     OCpar.grad_shapes[slot] = shapes[slot];
  }
  free_int_vector(shapes);  
  */
  OCpar.grad_shapes = shapes;
  
  return TCL_OK;
}

/****
 *  this specifically turns gradient mode on or off
 *  can be used only within 'gradient' procedure, the test is insufficient...
 *  (purpose: to enable target function evaluation from 'gradient' procedure)
 ****/
 int tclGradModeSwitch(ClientData data,Tcl_Interp* interp, int argc, char *argv[])
{
  /* test for existence of OCpar */
  if (OCpar.isinit != 1)
     return TclError(interp,"ERROR: gradient mode not activated! oc_gradmode can be called only within 'gradient' procedure");
     
  /* read input parameters (only RFshapes) */
  if ( argc != 2)
    return TclError(interp,"Usage: oc_gradmode on | off");

  if ( !strcmp(argv[1],"on") ) {
    OCpar.gradmode = 1;
    printf("GradMode turned  ON\n");
  } else if ( !strcmp(argv[1],"off") ) {
    OCpar.gradmode = 0;
    printf("GradMode turned OFF\n");
  }
  
  return TCL_OK;
}

/****
 * this command serves debugging gradient routine in input file
 * Procedure calculating gradients needs to be called either from oc_optimize 
 * or by help of this command. It makes sure that all necessary switches are set
 ****/ 
 int tclOCevaluate_grad(ClientData data,Tcl_Interp* interp, int argc, char *argv[])
{

  if (argc != 2)
    return TclError(interp,"usage: oc_evaluate_grad <name of procedure>");

  OCpar_initialize();
  test_pulseq_for_acqOC_prop(interp);
  OCpar.gradmode = 1;

  if (Tcl_Eval(interp,argv[1]) != TCL_OK) {
    fprintf(stderr,"oc_evaluate_grad error in evaluation of: %s\n",interp->result);
    OCpar_destroy();
    exit(1);
  }
  /* on succesfull evaluation, interp->result is already set by return value
     of 'gradient' procedure */  

  /* debugging savestate: */   
  move_OCpar_to_tcl(interp);
  
  OCpar_destroy();

  return TCL_OK;
}

/****
 * C implementation of oc_acq_hermit
 *   ( result saved in current FID point )
 ****/
void acqOC_hermit(void)
{
  complx c;
  
  if (puls->fid == NULL) {
    fprintf(stderr,"error: the 'oc_acq_hermit' command can only be used when the computing method is 'direct'\n");
    exit(1);  
  }
  
  if (puls->curr_nsig + 1 > LEN(puls->fid)) {
    fprintf(stderr,"error: oc_acq_hermit overflow in fid points\n");
    exit(1);
  }

  _evolve_with_prop();
  _reset_prop();

  c = cm_trace(puls->sigma,puls->fdetect);
  c.im = 0.0;
  
  puls->fid[++(puls->curr_nsig)] = c;
}

/****
 * C implementation of oc_acq_nonhermit
 *   ( result saved in current FID point )
 *   normalization removed...
 ****/
void acqOC_nonhermit(void)
{
  complx c;
  double r, n1, n2;
  
  if (puls->fid == NULL) {
    fprintf(stderr,"error: the 'oc_acq_nonhermit' command can only be used when the computing method is 'direct'\n");
    exit(1);  
  }
  
  if (puls->curr_nsig + 1 > LEN(puls->fid)) {
    fprintf(stderr,"error: oc_acq_nonhermit overflow in fid points\n");
    exit(1);
  }

  _evolve_with_prop();
  _reset_prop();

  /* normalize with sizes of start and detect operators -should we do?- */
  /* m_adjoint(puls->tmp,puls->fstart);
  c = m_trace(puls->tmp,puls->fstart);
  n1 = c.re;
  m_adjoint(puls->tmp,puls->fdetect);
  c = m_trace(puls->tmp, puls->fdetect);
  n2 = c.re;
  */

  c = cm_trace_adjoint(puls->fdetect,puls->sigma);
  r = c.re*c.re+c.im*c.im;
  /* r /= sqrt(n1);
     r /= sqrt(n2);
  */

  puls->fid[++(puls->curr_nsig)] = Complx(r,0.0);
}

/****
 * C implementation of oc_acq_prop
 *   ( result saved in current FID point )
 ****/
void acqOC_prop(int Ud)
{
  complx c;
  double r;
  
  if (puls->fid == NULL) {
    fprintf(stderr,"error: the 'oc_acq_prop' command can only be used when the computing method is 'direct'\n");
    exit(1);  
  }
  
  if (puls->curr_nsig + 1 > LEN(puls->fid)) {
    fprintf(stderr,"error: oc_acq_prop overflow in fid points\n");
    exit(1);
  }

  /* make the adjoint of desired propagator  */
  /* multiply it with optimized propagator and get the trace */
  c = cm_trace_adjoint(puls->STO[Ud], puls->U);

  r = c.re*c.re+c.im*c.im;
  
  puls->fid[++(puls->curr_nsig)] = Complx(r,0.0);
}


/****
 * Tcl implementation of oc_acq_hermit
 ****/
 int tclAcqOCHermit(ClientData data,Tcl_Interp* interp, int argc, char *argv[])
{

  check_pulse();
  /* acqOC can be called only once in pulseq, check it */
  if (puls->curr_nsig != 0)
     return TclError(interp,"oc_acq_hermit can be used only once in pulse sequence!!!");
  
  /* this does not accept any arguments */
  if (argc != 1) 
     return TclError(interp,"usage: oc_acq_hermit\n (oc_acq_hermit does not accept any parameters)");
  
  /* test whether to calculate target function or gradient */
  if (OCpar.gradmode) {
     /* printf("do gradients\n"); */
     gradOC_hermit();
  } else {
     /* printf("do target function\n"); */
     acqOC_hermit();
  }
    
  return TCL_OK;

}

/****
 * Tcl implementation of oc_acq_nonhermit
 ****/
 int tclAcqOCnonHermit(ClientData data,Tcl_Interp* interp, int argc, char *argv[])
{
  check_pulse();
  /* acqOC can be called only once in pulseq, check it */
  if (puls->curr_nsig != 0)
     return TclError(interp,"oc_acq_nonhermit can be used only once in pulse sequence!!!");
  
  /* this does not accept any arguments */
  if (argc != 1)
     return TclError(interp,"usage: oc_acq_nonhermit\n (oc_acq_nonhermit does not accept any parameters)");
  
  /* test whether to calculate target function or gradient */
  if (OCpar.gradmode) {
     /* printf("do gradients\n"); */
     gradOC_nonhermit();
  } else {
     /* printf("do target function\n"); */
     acqOC_nonhermit();
  }
  
  return TCL_OK;

}

/****
 * Tcl implementation of oc_acq_prop
 ****/
 int tclAcqOCProp(ClientData data,Tcl_Interp* interp, int argc, char *argv[])
{
  int ud;
  
  check_pulse();
  /* acqOC can be called only once in pulseq, check it */
  if (puls->curr_nsig != 0)
     return TclError(interp,"oc_acq_prop can be used only once in pulse sequence!!!");
  
  /* check arguments */
  if (argc != 2)
     return TclError(interp,"usage: oc_acq_prop <desired propagator>");
  if (Tcl_GetInt(interp,argv[1],&ud) != TCL_OK)
       return TclError(interp,"oc_acq_prop: can't read first argument, must be integer");
  if (puls->STO[ud] == NULL)
     return TclError(interp,"oc_acq_prop: desired propagator seems not to exist");
     
  /* test whether to calculate target function or gradient */
  if (OCpar.gradmode) {
     /* printf("do gradients\n"); */
     gradOC_prop(ud);
  } else {
     /* printf("do target function\n"); */
     acqOC_prop(ud);
  }

  return TCL_OK;

}

/****
 * Tcl implementation of oc_grad_add_energy_penalty [ units rad.s-1]
 ****/
 int tclGradPenalty(ClientData data,Tcl_Interp* interp, int argc, char *argv[])
{
  int fidN, slot, ndim, Nslot, i, ii, idx;
  int *slotvec;
  double *lamvec, lam, xx, yy;
  extern FD** fd;
  extern int nfd;
  double2 *gr;
  
  if (argc < 2 || ((argc % 2) != 0))
     return TclError(interp,"Usage: oc_grad_add_energy_penalty <fid> <RFshape> <lambda> ?<RFshape> <lambda>? ...");
     
  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR) 
     return TclError(interp,"oc_grad_add_energy_penalty: argument 1 must be integer <data set>");

  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL) 
    return TclError(interp,"oc_grad_add_energy_penalty: data set %d was not previously loaded\n",fidN);

  Nslot = argc/2-1;
  slotvec = int_vector(Nslot);
  lamvec = double_vector(Nslot);
  ii = 1;
  ndim = 0;
  for (i=1; i<=Nslot; i++) {
     ii++;
     if (Tcl_GetInt(interp,argv[ii],&slot) == TCL_ERROR) {
        free_int_vector(slotvec);
	free_double_vector(lamvec);
	return TclError(interp,"oc_grad_add_energy_penalty: argument %d must be integer <RFshape>",ii);
     }
     if (!RFshapes[slot]) {
        free_int_vector(slotvec);
	free_double_vector(lamvec);
	return TclError(interp,"oc_grad_add_energy_penalty: RFshape %d does not exist",slot);
     }
     ii++;
     if (Tcl_GetDouble(interp,argv[ii],&lam) == TCL_ERROR) {
        free_int_vector(slotvec);
	free_double_vector(lamvec);
	return TclError(interp,"oc_grad_add_energy_penalty: unable to get double <lambda> from %s",argv[ii]);
     }
     slotvec[i] = slot;
     lamvec[i] = lam;
     ndim += RFshapes_len(slot);
  }

  if ( ndim != fd[fidN]->np ) {
     free_int_vector(slotvec);
     free_double_vector(lamvec);
     return TclError(interp,"oc_grad_add_energy_penalty: mismatch in lengths of gradients and RF shapes");
  }
  
  gr = (double2*)(fd[fidN]->data);
  idx = 0;
  /* IMPORTANT: this assumes that gradients in fid are {gr_x, gr_y} !!! */
  for (i=1; i<=Nslot; i++) {
     slot = slotvec[i];
     for (ii=1; ii<=RFshapes_len(slot); ii++) {
        xx = (RFshapes[slot][ii].ampl)*cos( (RFshapes[slot][ii].phase)*DEG2RAD );
        yy = (RFshapes[slot][ii].ampl)*sin( (RFshapes[slot][ii].phase)*DEG2RAD );
        idx++;
	gr[idx].re += lamvec[i]*4.0*xx*M_PI;
	gr[idx].im += lamvec[i]*4.0*yy*M_PI;
     }
  }
  
  free_int_vector(slotvec);
  free_double_vector(lamvec);
  
  return TCL_OK;

}


/****
 * Tcl implementation of TF_line
 *      ( this is ment for debugging mainly!!! )
 ****/
 int tclTFLine(ClientData data,Tcl_Interp* interp, int argc, char *argv[])
{
  int gr_slot, i, j, N, Nsh, slot;
  int Nvars=0;
  double g1, g2, gs, dum;
  double *dir;
  extern FD** fd;
  
  /* check arguments */
  if (argc < 5)
     return TclError(interp,"usage: TF_line <from> <to> <step> <RFshape> ?<RFshape> ...?");
  if (Tcl_GetDouble(interp,argv[1],&g1) != TCL_OK)
       return TclError(interp,"TF_line: can't read second argument, must be double <from>");
  if (Tcl_GetDouble(interp,argv[2],&g2) != TCL_OK)
       return TclError(interp,"TF_line: can't read third argument, must be double <to>");
  if (Tcl_GetDouble(interp,argv[3],&gs) != TCL_OK)
       return TclError(interp,"TF_line: can't read fourth argument, must be double <step>");
  if (gs<=0.0)
     return TclError(interp,"TF_line: wrong step size!");

  Nsh = argc-4;
  if (Nsh<=0)
     return TclError(interp,"TF_line: there are no variable shapes!");
  OCpar_initialize();
  OCpar.var_shapes = int_vector(Nsh);
  OCpar.var_shapes_min = double_vector(Nsh);
  OCpar.var_shapes_max = double_vector(Nsh);
  OCpar.var_shapes_rmsmax = double_vector(Nsh);
  OCpar.grad_shapes = int_vector(Nsh);
  j = 0; 
  for (i=4; i<argc; i++) {
     if (Tcl_GetInt(interp,argv[i],&slot) != TCL_OK) {
          OCpar_destroy();
          return TclError(interp,"TF_line: can't read argument %d, must be integer <RFshape>",i);
     }
     if (!RFshapes[slot]) {
	  OCpar_destroy();
          return TclError(interp,"TF_line: argument %d -> RFshape does not exist", i);
     } 
     j++;
     OCpar.var_shapes[j] = slot;
     OCpar.var_shapes_min[j] = 0.0;
     OCpar.var_shapes_max[j] = 1.0e6;
     OCpar.var_shapes_rmsmax[j] = 1.0e6;
     OCpar.grad_shapes[j] = slot; 
     printf("reading RFshape slot %d\n",j);
  }

  test_pulseq_for_acqOC_prop(interp);
  
  /* evaluate gradient */
  gr_slot = evaluate_gradient(interp);
 printf("gradient done\n"); 
  /* store current RF shapes at safe place */
  RFelem *CURRshapes[Nsh];
  for (i=0; i<Nsh; i++) {
     N = RFshapes_len(OCpar.var_shapes[i+1]);
     Nvars += N;
     CURRshapes[i]=RFshapes_alloc(N);
     for (j=1; j<=N; j++) {
        CURRshapes[i][j].ampl = RFshapes[OCpar.var_shapes[i+1]][j].ampl;
        CURRshapes[i][j].phase = RFshapes[OCpar.var_shapes[i+1]][j].phase;
     }
  }
 printf("RFshapes stored separately\n");
 /* check length of gradients with total number of variables */
 if (Nvars != (fd[gr_slot]->np)) {
    fprintf(stderr,"TF_line error: number of gradients does not match total number of variables in RF shapes (%d versus %d)\n",Nvars,(fd[gr_slot]->np) );
    exit(1);
  }
 /* Nvars so far counts just elements in RF shapes. But vector dimension is double (aml, ph) */
 Nvars *=2;
printf("FID length and variables checked\n"); 
 /* initialize variable with gradient direction */
 dir = double_vector(Nvars);
 i = 0;
 for (j=1; j<=fd[gr_slot]->np; j++) {
     i++;
     dum = (fd[gr_slot]->data)[j].re;
     /* gradient calculated for maximizing, here we minimize. r[i]=-dum */
     dir[i] = dum;
     i++;
     dum = (fd[gr_slot]->data)[j].im;
     /* gradient calculated for maximizing, here we minimize. r[i]=-dum */
     dir[i] = dum;
 }
printf("gradient moved to direction\n");
 Tcl_ResetResult(interp);
 while (g1<=g2) {
    dum = func_to_min(g1,dir,CURRshapes,interp);
    TclAppendResult(interp,"%.15g",dum); /* doesn't work fo unknown reasons... */
    printf("   %.15g   %.15g\n",g1, dum);
    g1 += gs; 
 }
printf("calculations done\n");

  free_grad_fid(gr_slot,interp);
  free_double_vector(dir);
  for (i=0; i<Nsh; i++) {
    free((char *)CURRshapes[i]);
  } 
  OCpar_destroy();
  
  return TCL_OK;

}

/****
 *  function to replicate global C variable OCPar into Tcl environment
 *      (used in cluster calculation mode)
 ****/
void move_OCpar_to_tcl(Tcl_Interp* interp)
{
  char buf[8];
  int i;

  sprintf(buf,"%d",OCpar.isinit);
  /* printf("OCpar.isinit = %s\n",buf); */
  Tcl_SetVar2(interp,"__OCPAR","isinit",buf,TCL_GLOBAL_ONLY);
  
  sprintf(buf,"%d",OCpar.gradmode);
  /* printf("OCpar.gradmode = %s\n",buf); */
  Tcl_SetVar2(interp,"__OCPAR","gradmode",buf,TCL_GLOBAL_ONLY);

  sprintf(buf,"%d",OCpar.gradmodeprop);
  /* printf("OCpar.gradmodeprop = %s\n",buf); */
  Tcl_SetVar2(interp,"__OCPAR","gradmodeprop",buf,TCL_GLOBAL_ONLY);

  sprintf(buf,"%d",OCpar.propstatus);
  /* printf("OCpar.propstatus = '%s'\n",buf); */
  Tcl_SetVar2(interp,"__OCPAR","propstatus",buf,TCL_GLOBAL_ONLY); 

  sprintf(buf,"%d",OCpar.verb);
  /* printf("OCpar.verb = %s\n",buf); */
  Tcl_SetVar2(interp,"__OCPAR","verb",buf,TCL_GLOBAL_ONLY);

  /* OCpar.var_shapes */
  if (OCpar.var_shapes != NULL) {
     /* printf("OCpar.var_shapes = "); */
     if (Tcl_GetVar2Ex(interp, "__OCPAR", "var_shapes", TCL_GLOBAL_ONLY) != NULL) Tcl_UnsetVar2(interp,"__OCPAR","var_shapes",TCL_GLOBAL_ONLY);
     for (i=1; i<=LEN(OCpar.var_shapes); i++) {
        sprintf(buf,"%d",OCpar.var_shapes[i]);
	/* printf(" %s",buf); */
        Tcl_SetVar2(interp, "__OCPAR", "var_shapes", buf, TCL_GLOBAL_ONLY|TCL_LIST_ELEMENT|TCL_APPEND_VALUE);
     }
     /* printf("\n"); */
  }
       
  /* OCPar.grad_shapes */
  if (OCpar.grad_shapes != NULL) {
     if (Tcl_GetVar2Ex(interp, "__OCPAR", "grad_shapes", TCL_GLOBAL_ONLY) != NULL) Tcl_UnsetVar2(interp,"__OCPAR","grad_shapes",TCL_GLOBAL_ONLY);
     for (i=1; i<=LEN(OCpar.grad_shapes); i++) {
        sprintf(buf,"%d",OCpar.grad_shapes[i]);
        Tcl_SetVar2(interp, "__OCPAR", "grad_shapes", buf, TCL_GLOBAL_ONLY|TCL_LIST_ELEMENT|TCL_APPEND_VALUE);
     }
  }     
  
}

/****
 *  function to replicate Tcl variable __OCPAR into global C
 *      (used in cluster calculation mode)
 ****/
void create_OCpar_from_tcl(Tcl_Interp* interp)
{
  int c,i,p;
  char *src;
  Tcl_Obj *lptr, **elptr;
    
  if ((src=Tcl_GetVar2(interp,"__OCPAR","isinit",TCL_GLOBAL_ONLY)) != NULL) {
     if (Tcl_GetInt(interp,src,&c) != TCL_OK) {
        TclError(interp,"cannot get int from __OCPAR(isinit)");
	exit(0);
     }
     OCpar.isinit = c;
  }
  if ((src=Tcl_GetVar2(interp,"__OCPAR","gradmode",TCL_GLOBAL_ONLY)) != NULL) {
     if (Tcl_GetInt(interp,src,&c) != TCL_OK) {
        TclError(interp,"cannot get int from __OCPAR(gradmode)");
	exit(0);
     }
     OCpar.gradmode = c;
  }
  if ((src=Tcl_GetVar2(interp,"__OCPAR","gradmodeprop",TCL_GLOBAL_ONLY)) != NULL) {
     if (Tcl_GetInt(interp,src,&c) != TCL_OK) {
        TclError(interp,"cannot get int from __OCPAR(gradmodeprop)");
	exit(0);
     }
     OCpar.gradmodeprop = c;
  }
  if ((src=Tcl_GetVar2(interp,"__OCPAR","propstatus",TCL_GLOBAL_ONLY)) != NULL) {
     if (Tcl_GetInt(interp,src,&c) != TCL_OK) {
        TclError(interp,"cannot get int from __OCPAR(propstatus)");
	exit(0);
     }
     OCpar.propstatus = c;
  }
  if ((src=Tcl_GetVar2(interp,"__OCPAR","verb",TCL_GLOBAL_ONLY)) != NULL) {
     if (Tcl_GetInt(interp,src,&c) != TCL_OK) {
        TclError(interp,"cannot get int from __OCPAR(verb)");
	exit(0);
     }
     OCpar.verb = c;
  }
  if ((lptr=Tcl_GetVar2Ex(interp,"__OCPAR","var_shapes",TCL_GLOBAL_ONLY)) != NULL) {
     if (Tcl_ListObjGetElements(interp, lptr, &c, &elptr) != TCL_OK) {
	TclError(interp,"cannot decompose list in __OCPAR(var_shapes)");
	exit(0);
     }
     if (OCpar.var_shapes != NULL) {
        printf("move_OCoptPars_from_tcl is realocating OCpar.var_shapes\n");
        free_int_vector(OCpar.var_shapes);
     } 
     OCpar.var_shapes = int_vector(c);
     for (i=1;i<=c;i++) {
        if (Tcl_GetIntFromObj(interp,elptr[i-1],&p) != TCL_OK ) {
	   TclError(interp,"cannot get integer within __OCPAR(var_shapes)");
	   exit(0);
	}
	OCpar.var_shapes[i] = p;
     }
  }
  if ((lptr=Tcl_GetVar2Ex(interp,"__OCPAR","grad_shapes",TCL_GLOBAL_ONLY)) != NULL) {
     if (Tcl_ListObjGetElements(interp, lptr, &c, &elptr) != TCL_OK) {
	TclError(interp,"cannot decompose list in __OCPAR(grad_shapes)");
	exit(0);
     }
     if (OCpar.grad_shapes != NULL) {
        printf("move_OCoptPars_from_tcl is realocating OCpar.grad_shapes\n");
        free_int_vector(OCpar.grad_shapes);
     } 
     OCpar.grad_shapes = int_vector(c);
     for (i=1;i<=c;i++) {
        if (Tcl_GetIntFromObj(interp,elptr[i-1],&p) != TCL_OK ) {
	   TclError(interp,"cannot get integer within __OCPAR(grad_shapes)");
	   exit(0);
	}
	OCpar.grad_shapes[i] = p;
     }
  }

}

/****
 * Tcl implementation of test_move_OCpar
 *      ( this is ment for debugging mainly!!! )
 ****/
 int tclTest_move_OCpar(ClientData data,Tcl_Interp* interp, int argc, char *argv[])
{
  move_OCpar_to_tcl(interp);
  
  return TCL_OK;
}

/****
 * Tcl implementation of test_move_OCpar2
 *      ( this is ment for debugging mainly!!! )
 ****/
 int tclTest_move_OCpar2(ClientData data,Tcl_Interp* interp, int argc, char *argv[])
{
  int i;
  
  create_OCpar_from_tcl(interp);
  
  printf("OCpar.gradmode = %d\n",OCpar.gradmode);
  printf("OCpar.gradmodeprop = %d\n",OCpar.gradmodeprop);
  printf("OCpar.propstatus = %d\n",OCpar.propstatus);
  if (OCpar.var_shapes != NULL) {
     printf("OCpar.var_shapes = ");
     for (i=1; i<=LEN(OCpar.var_shapes); i++) {
        printf(" %d",OCpar.var_shapes[i]);
     }
     printf("\n");
  }
  if (OCpar.grad_shapes != NULL) {
     printf("OCpar.grad_shapes = ");
     for (i=1; i<=LEN(OCpar.grad_shapes); i++) {
        printf(" %d",OCpar.grad_shapes[i]);
     }
     printf("\n");
  }

  return TCL_OK;
}


void tclcmd_OCroutines(Tcl_Interp* interp) {

Tcl_CreateCommand(interp,"oc_optimize",tclOCoptimize2,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateCommand(interp,"oc_acq_hermit",tclAcqOCHermit,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateCommand(interp,"oc_acq_nonhermit",tclAcqOCnonHermit,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateCommand(interp,"oc_acq_prop",tclAcqOCProp,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateCommand(interp,"oc_evaluate_grad",tclOCevaluate_grad,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateCommand(interp,"oc_grad_shapes",tclGradShapes,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateCommand(interp,"oc_grad_add_energy_penalty",tclGradPenalty,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateCommand(interp,"oc_gradmode",tclGradModeSwitch,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

Tcl_CreateCommand(interp,"TF_line",tclTFLine,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateCommand(interp,"test_move_OCpar",tclTest_move_OCpar,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateCommand(interp,"test_move_OCpar2",tclTest_move_OCpar2,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

}
