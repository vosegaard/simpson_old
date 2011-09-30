/*
    Pulse propagation
    Copyright (C) 1999 Mads Bak, Jimmy T. Rasmussen
     Zdenek Tosner & Niels Chr. Nielsen: changes for optimal control
     
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
    
    Performs propagation with the pulses given in the 'pulseq' procedure
    given in the input file. Makes available several commands that
    can be read about in the manual.
    
    Called from sim.c that setup and perform the simulation.    
*/

#ifdef __APPLE__
#  include <malloc/malloc.h>
#else
#  include <malloc.h>
#endif
#include <stdlib.h>
#include <string.h>
#include "matrix_new.h"
#include "cm_new.h"
#include "pulse.h"
#include "tclutil.h"
#include "defs.h"
#include "rfshapes.h"
#include "iodata.h"
#include "OCroutines.h"
#include "B0inhom.h"
#include "relax.h"

/*
  tpropstart_usec : the time when the propagator was reset last time
  tpulsestart_usec : the time when the pulse sequence was started
  t_usec : the current time.
  tstartincr_usec : the time that will be added to the current time after resetting.
  t_proplength : the length of a propagator
 
  In this file t*_usec is in microseonds, the other time variables
  in seconds. t_usec alias 't' can be changed from the tcl script and must therefore be in microseconds
*/


void pulse_setpulsename(Pulse* P,char* pulsename)
{
  strcpy(P->pulsename,pulsename);
}

void pulse_initialize(Pulse* P,Tcl_Interp* interp,SpinSys* ss,double sw,double wr,double brl,int dor, double wr1,double brl1,double wr2,double brl2)
{
  int i,j;
  int N,NN;
  mv_complx *tmp, *tmpIy[100];
  const double done=1.0;
  
  strcpy(P->pulsename,"pulseq");
  P->ss=ss;
  P->interp=interp;
  P->wr=wr;
  P->sw=sw;
  P->dwellt=1.0/sw;
  P->brl=brl;
  P->brlorig=brl;
  P->gamma_add=0;
  P->dor=dor;
  P->wr1=wr1;
  P->brl1=brl1;
  P->wr2=wr2;
  P->brl2=brl2;
  P->nspins=ss->nspins;
  N=P->N=ss->matdim;
  P->nchan=ss->nchan;
  P->acq_adjoint=0; 

  P->t_usec=0.0;
  Tcl_UnlinkVar(interp,"t");
  if (Tcl_LinkVar(interp,"t",(char*)&(P->t_usec),TCL_LINK_DOUBLE | TCL_LINK_READ_ONLY) != TCL_OK) {
    fprintf(stderr,"error: unable to link variable t: %s\n",P->interp->result);
    exit(1);
  }

  P->tpropstart_usec=0.0;
  P->tpulsestart_usec=0.0;

  P->sigma = complx_matrix_alloc(N,N);
  P->fstart = complx_matrix_alloc(N,N);
  P->fdetect = complx_matrix_alloc(N,N);
  P->dU = complx_matrix_alloc(N,N);
  P->U = complx_matrix_alloc(N,N);
  P->tmpU = complx_matrix_alloc(N,N);
  P->Htot = double_matrix_alloc(N,N);
  NN = N*N;
  P->sumHrf = (double*)malloc(NN*sizeof(double));
  P->sumUph = (double*)malloc(N*sizeof(double));
  P->rfscalefactors = double_vector(ss->nchan); /* ZT: for rf inhomogeneity loop, old style */
  
  for (i=1;i<=P->nspins;i++) { 
    tmp = Iz(ss,i);
    P->Iz[i] = get_real_diag(tmp);
    complx_matrix_free(tmp);
    tmp = Ix(ss,i);
    P->Ix[i] = get_real(tmp);
    complx_matrix_free(tmp);
    tmpIy[i] = Iy(ss,i); 
  }

  for (i=0;i<MAXMATRIX;i++) {
    P->matrix[i]=NULL;
  }
  for (i=0;i<MAXSTO;i++) {
    P->STO[i]=NULL;
  }

  for (i=1;i<=P->nchan;i++) {
    P->chan_Ix[i] = (double*)malloc(NN*sizeof(double));
    memset(P->chan_Ix[i],0,NN*sizeof(double));
    P->chan_Iy[i] = (complx*)malloc(NN*sizeof(complx)); 
    memset(P->chan_Iy[i],0,NN*sizeof(complx));          
    P->chan_Iz[i] = (double*)malloc(N*sizeof(double));
    memset(P->chan_Iz[i],0,N*sizeof(double));

    for (j=1;j<=ss->nchanelem[i];j++) {
      daxpy_(&NN,&done,P->Ix[ss->chan[i][j]],&INTONE,P->chan_Ix[i],&INTONE);
      zaxpy_(&NN,&CPLX1,tmpIy[ss->chan[i][j]]->data,&INTONE,P->chan_Iy[i],&INTONE);
      daxpy_(&N,&done,P->Iz[ss->chan[i][j]],&INTONE,P->chan_Iz[i],&INTONE);
    }
  }
  P->tmpIz = (double*)malloc(N*sizeof(double));
  
  for (i=1;i<=P->nspins;i++) { 
    complx_matrix_free(tmpIy[i]);
  }
  /* ZT: offsets due to static field inhomogeneity */
  P->inhom_offset = NULL;
}

void pulse_destroy(Pulse* P)
{
  int i;
  for (i=0;i<MAXMATRIX;i++) {
    if (P->matrix[i]) {
      complx_matrix_free(P->matrix[i]);
    }
  }
  for (i=0;i<MAXSTO;i++) {
    if (P->STO[i]) {
      complx_matrix_free(P->STO[i]);
    }
  }
  
  double_matrix_free(P->Htot);
  complx_matrix_free(P->tmpU);
  complx_matrix_free(P->dU);
  complx_matrix_free(P->U);
  complx_matrix_free(P->sigma);
  complx_matrix_free(P->fstart);
  complx_matrix_free(P->fdetect);
  free((char*)P->sumHrf);
  free((char*)P->sumUph);
  free_double_vector(P->rfscalefactors);    /* ZT: for rf inhomogeneity loop, old style */
  
  for (i=1;i<=P->nspins;i++) {    
    free((char*)P->Iz[i]);
    free((char*)P->Ix[i]);
  }
  for (i=1;i<=P->nchan;i++) {
    free((char*)P->chan_Ix[i]);
    free((char*)P->chan_Iy[i]);   
    free((char*)P->chan_Iz[i]);
  }
  free((char*)P->tmpIz);
  
  /* ZT: offsets due to static field inhomogeneity */
  if (P->inhom_offset) free((char*)P->inhom_offset);

}

void pulse_setprop(Pulse* P,mv_complx *fstart,mv_complx *fdetect)
{
  cmv_copy(P->fstart,fstart);
  cmv_copy(P->fdetect,fdetect);  
}

Pulse* puls = NULL;

void pulse_propagate(Pulse* P,Hamilton* H,mv_complx *U,double dt,double t,complx* fid)
{
  int i;

  if (fid) {
    P->curr_nsig=0;
    P->fid=fid;
    m_zerov(P->fid);
    cmv_copy(P->sigma,P->fstart);
  }

  for (i=1;i<=P->nchan;i++) {
    P->phv[i]=0.0;
    P->rfv[i]=0.0;
  }

  P->H = H;
  if (H->isdiag) {
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
  
/*  P->dtmax= dt/2.0; */

  if (P->wr == 0.0)
    P->dtmax=1e9; /* dtmax is infinity in the static case. */
  else
    P->dtmax= 1.0e-6; /* dt max is 1 usec per default in case of spining*/
  P->dt= dt;
  P->tpulsestart_usec=t*1.0e6;
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

  if (puls) {
    fprintf(stderr,"error: pulse_propagate() was called before without destroying\n");
    exit(1);
  }
  puls=P;
  if (Tcl_Eval(P->interp,P->pulsename) != TCL_OK) {
    fprintf(stderr,"error: in evaluation of pulse sequence: %s\n",P->interp->result);
    exit(1);
  }
  if (P->U->col == 1) 
     cm_copydiagmv(U,P->U);
  else
     cmv_copy(U,P->U);
  P->H = NULL;
  puls = NULL;
  
  /* ZT: when calculating gradients, reset mx_pos when crystallite calculation finishes */
  if (OCpar.gradmode) OCpar.mx_pos=0;

}


/* Internal routines to be called from the TCL routines
   using 'puls' as a global Pulse structure points
 */

/* ZT: these functions moved here from bottom of this file */
int string2matrix(Tcl_Interp* interp,mv_complx** mptr, char* str)
{
  int i,j,nrows,ncols,nelem;
  char **rows,**cols,**elem;
  mv_complx * m = NULL;
  complx *data;

  *mptr=NULL;  

  Tcl_ResetResult(interp);
  if (Tcl_SplitList(interp,str,&nrows,&rows) != TCL_OK)
    return TCL_ERROR;

  for (i=0;i<nrows;i++) {
    if (Tcl_SplitList(interp,rows[i],&ncols,&cols) != TCL_OK) {
      free(rows);
      return TCL_ERROR;
    }
    if (!m) {
      if (nrows != ncols) {
         free(rows);
         free(cols);
         return TclError(interp,"matrix: matrix must be quadratic");
      }        
      m = complx_matrix_alloc(nrows,ncols);
      data = m->data;
    }
    for (j=0;j<ncols;j++) {
      if (Tcl_SplitList(interp,cols[j],&nelem,&elem) != TCL_OK) {
        free(rows);
        free(cols);
        return TCL_ERROR;
      }
      if (nelem == 1 || nelem == 2) {
        if (Tcl_GetDouble(interp,elem[0],&(data[i+j*nrows].re)) != TCL_OK) {
          free(rows);
          free(cols); 
          free(elem); 
          return TCL_ERROR;
        }
      } else {      
         free(rows);
         free(cols);
         free(elem);
         return TclError(interp,"matrix: matrix element must contain one or two real numbers");          
      }
      if (nelem == 2) {
        if (Tcl_GetDouble(interp,elem[1],&(data[i+j*nrows].im)) != TCL_OK) {
          free(rows);
          free(cols); 
          free(elem); 
          return TCL_ERROR;
        }
      } else {
        data[i+j*nrows].im=0.0;      
      }
      free(elem);
    }
    free(cols);
  }
  free(rows);
  *mptr = m;
  return TCL_OK; 
}


mv_complx * get_matrix_1(Tcl_Interp* interp,char* name)
{
  int num, dupl=1;
  mv_complx * m;
  
  if (Tcl_GetInt(interp,name,&num) == TCL_OK) {    
    if (num < 1 || num > MAXMATRIX) {
       TclError(interp,"matrix: called with number larger than %d or smaller than 1\n",MAXMATRIX);
       return NULL;
    }
    if (!puls->matrix[num]) {
       TclError(interp,"matrix number %d is undefined\n",num);
       return NULL;
    }
    m=puls->matrix[num];
  } else if (!strcmp(name,"start")) {
    m=puls->fstart; 
  } else if (!strcmp(name,"detect")) {
    m=puls->fdetect; 
  } else if (!strcmp(name,"density")) {
    m=puls->sigma; 
  } else if (!strcmp(name,"propagator")) {
    if (puls->U->col == 1) {
       m = complx_matrix_alloc(puls->N,puls->N);
       cm_copydiagmv(m,puls->U);
       dupl = 0;
    } else {
       m=puls->U;
    }
  } else if (!strcmp(name,"hamiltonian")) {
    ham_rotate2(puls->H,puls->wr*puls->t_usec*1.0e-6*(180.0/M_PI)+puls->gamma_add,puls->brl);
    ham_hamilton(puls->H);
    m = complx_matrix_alloc(puls->N,puls->N);
    if (puls->H->isdiag) {
       mv_complx * dum;
       dum = complx_vector_alloc(puls->N);
       dmv_complx(puls->H->H,dum);
       cm_copydiagmv(m,dum);
       complx_vector_free(dum);
    } else {
       dmv_complx(puls->H->H,m);
    }
    dupl = 0;
  } else if (!strcmp(name,"avgham")) {
    double dt = puls->t_usec - puls->tpropstart_usec;
    m = complx_matrix_alloc(puls->N,puls->N);
    if (dt < 1e-6) {
      cmv_zero(m);
    } else {
      m = cm_ln(puls->U);
      cmv_mulc(m, Complx(0.0, 1.0e6/dt));
/*      printf("dt: %g\n", dt); */
    }
    dupl = 0;
  } else if (!strcmp(name,"avghamtv")) {
    m=puls->H->H;
    
  } else {
    TclError(interp,"error: matrix: argument must be 'start', 'detect', 'density', 'propagator', 'avgham', or "
                    "'hamiltonian'\n");
    return NULL;
  }
  Tcl_ResetResult(interp);
  
  if (dupl) {
    return cmv_dup(m);
  } else {
    return m;
  }
}

mv_complx* get_matrix_2(Tcl_Interp* interp,char* name,char* list)
{
  mv_complx * m;
  int N;
  
  if (!strcmp(name,"operator")) {
    return ss_readoper(puls->ss,list);
  }
  N=puls->N;
  m = complx_matrix_alloc(N,N);
  cmv_zero(m);

  if (!strcmp(name,"notelements") || !strcmp(name,"elements")) {
      char **list1,**list2;
      int i,j,nlist1,nlist2,val1,val2;
      
      if (Tcl_SplitList(interp,list,&nlist1,&list1) != TCL_OK)
        return NULL;

      for (i=0;i<nlist1;i++) {
        if (Tcl_SplitList(interp,list1[i],&nlist2,&list2) != TCL_OK) {
          Tcl_Free((char *) list1);
          return NULL;
        }
        if (nlist2 != 2) {
          Tcl_Free((char *) list1);
          Tcl_Free((char *) list2);
          TclError(interp,"matrix: expecting two elements (matrix row and column) in list");
          return NULL;
        }
        if (Tcl_GetInt(interp,list2[0],&val1) != TCL_OK) {
          Tcl_Free((char *) list1);
          Tcl_Free((char *) list2); 
          return NULL;
        }
        if (Tcl_GetInt(interp,list2[1],&val2) != TCL_OK) {
          Tcl_Free((char *) list1);
          Tcl_Free((char *) list2); 
          return NULL;
        }
        if ( val1 < 1 || val1 > N || val2 < 1 || val2 > N) {
          Tcl_Free((char *) list1);
          Tcl_Free((char *) list2); 
          TclError(interp,"matrix: value out of range in list");
          return NULL;
        }
        m->data[val1-1+N*(val2-1)].re=1.0;
        Tcl_Free((char *) list2);        
      }
      if (!strcmp(name,"notelements")) {
	for (i=0; i<N*N; i++) {
	   if (m->data[i].re != 0.0)
	      m->data[i].re = 0.0;
	   else
	      m->data[i].re = 1.0;
	}
      }
      Tcl_Free((char *) list1);
  } else if (!strcmp(name,"totalcoherence")) {
    int coh[MAXSPINS+1];
    char **list1;
    int nlist1,i,j,k;    
    mv_complx * Q;
    
    if (Tcl_SplitList(interp,list,&nlist1,&list1) != TCL_OK)
      return NULL;
    if (nlist1 > MAXSPINS) {      
      TclError(interp,"matrix: internal array overflow\n");
      return NULL;
    }

    for (i=0;i<nlist1;i++) {
      if (Tcl_GetInt(interp,list1[i],&coh[i]) != TCL_OK) {
        Tcl_Free((char *) list1);
        return NULL;
      }
    }
    Q=Iqdelta(puls->ss);

    for (i=0;i<N;i++) {
      for (j=0;j<N;j++) {
        int qdelta = (int)( Q->data[i+j*N].re );
        for (k=0;k<nlist1;k++) {
          if (coh[k] == qdelta) {
             m->data[i+j*N].re = 1.0;
             break;
          }
        }
      }
    }
    Tcl_Free((char *) list1);
    complx_matrix_free(Q);
  } else if (!strcmp(name,"coherence")) {
    char **list1,**list2;
    int i,j,nlist1,nlist2;      
    double coh[MAXSPINS+1];
    mv_complx *tmp;


    if (Tcl_SplitList(interp,list,&nlist1,&list1) != TCL_OK)
      return NULL;

    for (i=0;i<nlist1;i++) {
      if (Tcl_SplitList(interp,list1[i],&nlist2,&list2) != TCL_OK) {
        Tcl_Free((char *) list1);
        return NULL;
      }
      if (nlist2 != puls->nspins) {
        Tcl_Free((char *) list1);
        Tcl_Free((char *) list2);
        TclError(interp,"matrix: expecting %d elements (number of spins)",puls->nspins);
        return NULL;
      }
      for (j=0;j<nlist2;j++) {
        if (Tcl_GetDouble(interp,list2[j],&coh[j+1]) != TCL_OK) {
          Tcl_Free((char *) list1);
          Tcl_Free((char *) list2);
          return NULL;
        }
      }
      tmp=Icoherence(puls->ss,coh);
      cm_or(m,m,tmp);      
      complx_matrix_free(tmp);  
      Tcl_Free((char *) list2);
    }
    Tcl_Free((char *) list1);
    
  } else if (!strcmp(name,"list")) {
    char **list1,**list2,**list3;
    int i,j,nlist1,nlist2,nlist3;

    if (Tcl_SplitList(interp,list,&nlist1,&list1) != TCL_OK) {
      TclError(interp,"matrix: specify list");
      return NULL;
    }
    if (nlist1 != N) {
      TclError(interp,"matrix: expecting %d elements (matrix dimension)",
	  N);
      return NULL;
    }

    for (i=0;i<nlist1;i++) {
      if (Tcl_SplitList(interp,list1[i],&nlist2,&list2) != TCL_OK) {
        Tcl_Free((char *) list1);
        TclError(interp,"matrix: specify list (in list)");
        return NULL;
      }
      if (nlist2 != N) {
        Tcl_Free((char *) list1);
        Tcl_Free((char *) list2);
        TclError(interp,"matrix: expecting %d elements (matrix dimension)",
	  N);
        return NULL;
      }
      for (j=0;j<nlist2;j++) {
        if (Tcl_SplitList(interp,list2[j],&nlist3,&list3) != TCL_OK) {
	  Tcl_Free((char *) list1);
	  Tcl_Free((char *) list2);
          TclError(interp,"matrix: specify {re im} for each element");
	  return NULL;
	}
	if (nlist3 != 1 && nlist3 != 2) {
	  Tcl_Free((char *) list1);
	  Tcl_Free((char *) list2);
	  Tcl_Free((char *) list3);
          TclError(interp,"matrix: specify {re im} for each element");
	  return NULL;
	}
        if (Tcl_GetDouble(interp,list3[0],&(m->data[i+j*N].re)) != TCL_OK) {
          Tcl_Free((char *) list1);
          Tcl_Free((char *) list2);
	  Tcl_Free((char *) list3);
          TclError(interp,"matrix: real element (%d,%d) does not appear to be a double",
	    i,j);
          return NULL;
        }
        if (nlist3 == 2) {
	  if (Tcl_GetDouble(interp,list3[1],&(m->data[i+j*N].im)) != TCL_OK) {
            Tcl_Free((char *) list1);
            Tcl_Free((char *) list2);
  	    Tcl_Free((char *) list3);
            TclError(interp,"matrix: imaginary element (%d,%d) does not appear to be a double",
	      i,j);
            return NULL;
          }
	}
	Tcl_Free((char *) list3);
      }
      Tcl_Free((char *) list2);
    }
    Tcl_Free((char *) list1);
  } else {
    TclError(interp, "Usage: matrix: arguments must be"
                            "  elements {{i j} {i j} ..}\n"
                            "  notelements {{i j} {i j} ..}\n"                            
                            "  totalcoherence {dm1 dm2 ..}\n"
                            "  coherence {dm1 dm2 dmN} {dm1 dm2 dmN} ..., where N is number of nuclei\n"
                            "  list {{c11 c12 ..} {c21 c22 ..} ..}, where c is a real number or {re im}\n"
			    );
    return NULL;
  }
  return m;
}
/* ZT: these functions (above) moved here from bottom of this file */

void _offset(int nchan,double* offset)
{
  int i,used;
  double dof;

  if (nchan < 1 || nchan > puls->nchan) {
    fprintf(stderr,"error: illegal number of channels '%d' in offset\n",nchan);
    exit(1);
  }

  /* ZT: here we need to take care of possibly existing static field inhom. offsets */
  if (!puls->inhom_offset) {
     /* do the old stuf since field is perfect */
     used=0;  
     memset(puls->tmpIz,0,puls->N*sizeof(double));
     for (i=1;i<=nchan;i++) {
       if (offset[i] != 0.0) {
         used=1;
	 dof = offset[i]*M_PI*2.0;
	 daxpy_(&(puls->N),&dof,puls->chan_Iz[i],&INTONE,puls->tmpIz,&INTONE);
       }
     }
     ham_set_offset(puls->H,puls->tmpIz,used);
  } else {
     /* add given offsets to static field inhomogeneity offsets */
     used=2; 
     memset(puls->tmpIz,0,puls->N*sizeof(double));
     for (i=1;i<=nchan;i++) {
        if (offset[i] != 0.0) {
	  used = 1;
	}
	dof = offset[i]*M_PI*2.0+puls->inhom_offset[i];
	daxpy_(&(puls->N),&dof,puls->chan_Iz[i],&INTONE,puls->tmpIz,&INTONE);
     }
     ham_set_offset(puls->H,puls->tmpIz,used);
  }
}

/****
 * ZT: modification in scalling rf power by rf inhomogeneity factor
 ****/
void _rf(int channel,double rffield)
{
  if (channel < 1 || channel > puls->nchan) {
    fprintf(stderr,"error: channel %d out of range in rf()\n",channel);
    exit(1);
  }
  puls->rfv[channel]=rffield*M_PI*2.0*(puls->rfscalefactors[channel]);
}

void _ph(int channel,double phase)
{
  if (channel < 1 || channel > puls->nchan) {
    fprintf(stderr,"error: channel %d out of range in rf()\n",channel);
    exit(1);
  }
  puls->phv[channel]=phase*(M_PI/180.0);
}

void _reset_prop()
{
  puls->tpropstart_usec= puls->t_usec;

  if (puls->H->isdiag) {
     /* limit initial propagator to diagonal, rest of allocated space is not used
        but ready to be used when propagator turns to full form due to pulse    */
     puls->U->col = 1;
     cv_d1(puls->U);
  } else {
     puls->U->col = puls->U->row;
     cm_unit(puls->U);
  }
  puls->Uisunit=1;
  puls->waspulse=0;
  puls->isselectivepulse=0;
  puls->cannotbestored=0;
}

void _evolve_with_prop()
{
  if (!puls->Uisunit) {
    /* if no pulse was called and the hamiltonian is diagonal we can use this command */
    /* if (!puls->waspulse && puls->H->isdiag)
         m_simtrans_diagprop(puls->sigma,puls->U);
       else
         m_simtrans(puls->sigma,puls->U);
     */ 
     if (puls->U->col == 1)
        simtrans_diagprop3(puls->sigma,puls->U);
     else
        simtrans(puls->sigma,puls->U);
  }
}

void _reset(double tstartincr_usec)
{
  cmv_copy(puls->sigma,puls->fstart);
  puls->t_usec=puls->tpulsestart_usec+tstartincr_usec;
  _reset_prop();
  /* ZT: here we need to take care of possibly existing static field inhom. offsets */
  if (!puls->inhom_offset) {
     /* do the old stuf since field is perfect */
     ham_set_offset(puls->H,NULL,0); 
     ham_turnon(puls->H,"all");
  } else {
     int i;
     ham_turnon(puls->H,"all");
     /* set static field inhomogeneity offsets */
     memset(puls->tmpIz,0,puls->N*sizeof(double));
     for (i=1;i<=puls->nchan;i++) {
        daxpy_(&(puls->N),&(puls->inhom_offset[i]),puls->chan_Iz[i],&INTONE,puls->tmpIz,&INTONE);
     }
     ham_set_offset(puls->H,puls->tmpIz,2);
  }
}


void _filter(int num)
{
  mv_complx * m;
  int i,j,N;
  complx *sl1, *sl2;

  N = puls->N;
  m = puls->matrix[num];

  if (!m) {
    fprintf(stderr,"error: illegal filter number %d\n",num);
    exit(1);
  }
  if ( m->row != N || m->col != N) {
    fprintf(stderr,"error: filter: matrix %d must be a %dx%d matrix\n",num,N,N);
    exit(1);
  }

  _evolve_with_prop();
  _reset_prop();
  puls->cannotbestored=1;
  
  sl1 = puls->sigma->data;
  sl2 = m->data;
  for (i=0;i<N*N;i++) {
    if (sl2->re == 0.0 && sl2->im == 0.0) *sl1 = CPLX0;
    sl1++;
    sl2++;
  }
}

void _store(int num)
{
  if (puls->cannotbestored) {
    fprintf(stderr,"error: the command 'filter' cannot be stored\n");
    exit(-1);
  }

  if (puls->STO[num] != NULL)
     complx_matrix_free(puls->STO[num]);

  puls->STO[num] = cmv_dup(puls->U);

  puls->STO_tproplength_usec[num] = puls->t_usec - puls->tpropstart_usec;   
  puls->STO_tpropstart_usec[num] = puls->tpropstart_usec;

  /* save whether a pulse was applied */
  puls->STO_waspulse[num]=puls->waspulse;
  puls->STO_brl[num]=puls->brl;
}


void update_prop(mv_complx * dU)
{
   int N1, N2;
   static mv_complx * wsp=NULL;
   
   if ( puls->Uisunit) {
      puls->U->col = dU->col;
      cmv_copy(puls->U,dU);
      return;
   }
   N1 = puls->U->col;
   N2 = dU->col;
   if ( N1 == 1 && N2 == 1) {
      /* both diagonal */
      wsp = cmv_static(wsp,dU->row,1);
      cmv_mul_elem(wsp,puls->U,dU);
      cmv_copy(puls->U,wsp);
      return;
   }
   if ( N1 == 1) {
      /* U is diagonal, dU is full */
      wsp = cmv_static(wsp,dU->row,N2);
      cm_mul_fd(wsp,dU,puls->U);
      puls->U->col = puls->N;
      cmv_copy(puls->U,wsp);
      return;
   }
   if ( N2 == 1) {
      /* U is full, dU is diagonal */
      wsp = cmv_static(wsp,puls->U->row,N1);
      cm_mul_df(wsp,dU,puls->U);
      cmv_copy(puls->U,wsp);
      return;
   }
   /* both are full */
   wsp = cmv_static(wsp,puls->U->row,N1);
   cm_mul(wsp,dU,puls->U);
   cmv_copy(puls->U,wsp);
   return;
}

void _prop(int num,int ntimes)
{
  int i;
  mv_complx *m1, *m2;

  if (!puls->STO[num]) {
    fprintf(stderr,"error: illegal prop number %d\n",num);
    exit(1);
  }
  if (ham_ischanged(puls->H)) {
    fprintf(stderr,"error: 'store', 'turnon' and 'turnoff' will not have any"
                   " effect\n on a stored propagator\n"
                   "Call 'reset' before using stored propagators\n");
                   
    exit(1);
  }

  /* ZT: added condition about non-negativity of propagator start time
   *     negative start time is used with propagators created from matrices,lists...
   *     for which I don't want to check start time
   */
  if ( (puls->wr != 0.0) && (puls->STO_tpropstart_usec[num] >= 0.0) && !(VARIOUS_NOCHECKPROPTIME & various)) {  
     /* check if the time the propagator is applied is the
        time it was calculated. Time is relative to the
        rotor period because the Hamiltonian is periodic
        with the rotor period. */

     double tr,tdiff,t,tpropstart;

     tr=M_PI*2.0e6/puls->wr;
     tpropstart=puls->STO_tpropstart_usec[num];
     t=puls->t_usec;

     tdiff=t-tpropstart;

     if (fabs(tdiff-floor(tdiff/tr+0.5)*tr) > 1e-5) {
       fprintf(stderr,"error: a propagator was calculated at time %.10gusec\nrelative to the"
       " start of the rotor period but reused at time %.10gusec.\n",
        tpropstart-floor(tpropstart/tr)*tr,
        t-floor(t/tr)*tr );
        exit(1);
     }
     if (fabs(puls->STO_brl[num]-puls->brl) > 1e-5) {
       fprintf(stderr,"error: a propagator was calculated using a rotor angle of %g\nand reused with a rotor angle of %g\n",
         puls->STO_brl[num],puls->brl);
       exit(1);
     }
  }

  /* does the current propagator contain any pulses ? */  
  puls->waspulse = puls->STO_waspulse[num] || puls->waspulse; 

  if ( puls->STO[num]->col == 1 ) {
     puls->dU->col = 1;
     puls->tmpU->col = 1;
     cmv_copy(puls->dU,puls->STO[num]);
     for (i=1; i<ntimes; i++) {
        if (i%2 == 0) {
           m1 = puls->tmpU;
           m2 = puls->dU;
	} else {
           m1 = puls->dU;
           m2 = puls->tmpU;
	}
        cmv_mul_elem(m2, puls->STO[num], m1);
     }
  } else {
     puls->dU->col = puls->STO[num]->col;
     puls->tmpU->col = puls->STO[num]->col;
     cmv_copy(puls->dU,puls->STO[num]);
     for (i=1; i<ntimes; i++) {
        if (i%2 == 0) {
           m1 = puls->tmpU;
           m2 = puls->dU;
	} else {
           m1 = puls->dU;
           m2 = puls->tmpU;
	}
        cm_mul(m2, puls->STO[num], m1);
     }
  }
  if ( ntimes%2 == 0 ) 
     update_prop(puls->tmpU);
  else 
     update_prop(puls->dU);

  puls->t_usec += puls->STO_tproplength_usec[num] * ntimes;
  puls->Uisunit=0;  

}


void _acq(double phase)
{
  complx *ptr;
  
  if (puls->fid == NULL) {
    fprintf(stderr,"error: the 'acq' command can only be used when the computing method is 'direct'\n");
    exit(1);  
  }
  
  if (puls->curr_nsig + 1 > LEN(puls->fid)) {
    fprintf(stderr,"error: acq overflow in fid points\n");
    exit(1);
  }

  /* in relaxation mode rho is already evolved */
  if (!(puls->is_relax)) { 
     _evolve_with_prop();
  }
  _reset_prop();

  /* puls->fid[++(puls->curr_nsig)] = m_trace(puls->sigma,puls->fdetect); */
  ptr = &(puls->fid[++(puls->curr_nsig)]);
  if (puls->acq_adjoint == 0) {
    *ptr = cm_trace(puls->fdetect,puls->sigma);
  } else {
    *ptr = cm_trace_adjoint(puls->fdetect,puls->sigma);
  }

  if (phase != 0.0) {
      double cosph,sinph,re;

      phase *= DEG2RAD;
      cosph=cos(phase);
      sinph=sin(phase);
      re=ptr->re;
      ptr->re=cosph*re+sinph*ptr->im;
      ptr->im=-sinph*re+cosph*ptr->im;
  }
}

/***
 * ZT: acq with adjoint of detect operator and normalization
 ***/
void _acq2(double phase)
{
  complx acqnorm, *ptr;
  double dnorm;

  if (puls->fid == NULL) {
    fprintf(stderr,"error: the 'acq' command can only be used when the computing method is 'direct'\n");
    exit(1);  
  }
  
  if (puls->curr_nsig + 1 > LEN(puls->fid)) {
    fprintf(stderr,"error: acq overflow in fid points\n");
    exit(1);
  }

  /* in relaxation mode sigma is already evolved */
  if (!(puls->is_relax)) { 
     _evolve_with_prop();
  }
  _reset_prop();

  /* normlize according to starting as well as detect operators */
  acqnorm = cm_trace_adjoint(puls->fstart,puls->fstart);
  dnorm = sqrt(acqnorm.re);
  acqnorm = cm_trace_adjoint(puls->fdetect, puls->fdetect);
  dnorm *= sqrt(acqnorm.re);
  ptr = &(puls->fid[++(puls->curr_nsig)]);
  *ptr = CRmul(cm_trace_adjoint(puls->fdetect,puls->sigma),1.0/dnorm);
  if (phase != 0.0) {
      double cosph,sinph,re;

      phase *= DEG2RAD;
      cosph=cos(phase);
      sinph=sin(phase);
      re=ptr->re;
      ptr->re=cosph*re+sinph*ptr->im;
      ptr->im=-sinph*re+cosph*ptr->im;
  }
}


void _nacq(int np,int num,double phase)
{
  int k;
  double cosph=0,sinph=0,re;
  complx *ptr;


  if (puls->fid == NULL) {
    fprintf(stderr,"error: the 'acq' command can only be used when the computing method is 'direct'\n");
    exit(1);  
  }
  if (!puls->STO[num]) {
    fprintf(stderr,"error: illegal prop number %d\n",num);
    exit(1);
  }
  if (puls->curr_nsig + np > LEN(puls->fid)) {
    fprintf(stderr,"error: acq overflow in fid points\n");
    exit(1);
  }

  if (ham_ischanged(puls->H)) {
    fprintf(stderr,"error: 'store', 'turnon' and 'turnoff' will not have any"
                   " effect\n on a stored propagator (indirectly used"
                   " by 'acq <prop> <ntimes>' in this case)\n"
                   "Call 'reset' before using stored propagators\n");
    exit(1);
  }
  _acq(phase);

  puls->waspulse = puls->STO_waspulse[num]; 
  /*cmv_copy(puls->U,puls->STO[num]); */
  update_prop(puls->STO[num]);

  if (phase != 0.0) {
    phase *= DEG2RAD;
    cosph=cos(phase);
    sinph=sin(phase);
  }
  ptr = &(puls->fid[(puls->curr_nsig)]);
  for (k=2;k<=np;k++) {
     /* ZT: added condition about non-negativity of propagator start time
      *     negative start time is used with propagators created from matrices,lists...
      *     for which I don't want to check start time
      */
      if ( (puls->wr != 0.0) && (puls->STO_tpropstart_usec[num] >= 0.0) && !(VARIOUS_NOCHECKPROPTIME & various)) {  

         double tr,tdiff,t,tpropstart;

         tr=M_PI*2.0e6/puls->wr;
         tpropstart=puls->STO_tpropstart_usec[num];
         t=puls->t_usec;

         tdiff=t-tpropstart;

         if (fabs(tdiff-floor(tdiff/tr+0.5)*tr)  > 1e-5) {
           fprintf(stderr,"error: a propagator was calculated at time %.10gusec\nrelative to the"
           " start of the rotor period but reused at time %.10gusec.\n",
            tpropstart-floor(tpropstart/tr)*tr,
            t-floor(t/tr)*tr );
            exit(1);
         }
      }

      puls->t_usec += puls->STO_tproplength_usec[num];
      puls->Uisunit=0;  


      if ( puls->U->col == 1) {
        simtrans_diagprop3(puls->sigma,puls->U);
      } else {
        simtrans(puls->sigma,puls->U);
      }

      ptr++; (puls->curr_nsig)++;
      if (puls->acq_adjoint == 0) {
        *ptr = cm_trace(puls->fdetect,puls->sigma);
      } else {
        *ptr = cm_trace_adjoint(puls->fdetect,puls->sigma);
      }
  
      if (phase != 0.0) {
          re=ptr->re;
          ptr->re=cosph*re+sinph*ptr->im;
          ptr->im=-sinph*re+cosph*ptr->im;
      }
  }
  _reset_prop();
}

/***
 * ZT: acq with adjoint of detect operator and normalization
 ***/
void _nacq2(int np,int num,double phase)
{
  int k;
  double cosph=0,sinph=0,re, dnorm;
  complx acqnorm, *ptr;

  if (puls->fid == NULL) {
    fprintf(stderr,"error: the 'acq' command can only be used when the computing method is 'direct'\n");
    exit(1);  
  }
  if (!puls->STO[num]) {
    fprintf(stderr,"error: illegal prop number %d\n",num);
    exit(1);
  }
  if (puls->curr_nsig + np > LEN(puls->fid)) {
    fprintf(stderr,"error: acq overflow in fid points\n");
    exit(1);
  }

  if (ham_ischanged(puls->H)) {
    fprintf(stderr,"error: 'store', 'turnon' and 'turnoff' will not have any"
                   " effect\n on a stored propagator (indirectly used"
                   " by 'acq <prop> <ntimes>' in this case)\n"
                   "Call 'reset' before using stored propagators\n");
    exit(1);
  }
  
  _acq2(phase);

  puls->waspulse = puls->STO_waspulse[num]; 

  cmv_copy(puls->U,puls->STO[num]);

  if (phase != 0.0) {
    phase *= DEG2RAD;
    cosph=cos(phase);
    sinph=sin(phase);
  }
  ptr = &(puls->fid[(puls->curr_nsig)]);

  /* normlize according to starting as well as detect operators */
  acqnorm = cm_trace_adjoint(puls->fstart,puls->fstart);
  dnorm = sqrt(acqnorm.re);
  acqnorm = cm_trace_adjoint(puls->fdetect, puls->fdetect);
  dnorm *= sqrt(acqnorm.re);

  for (k=2;k<=np;k++) {

      if (puls->wr != 0.0 && !(VARIOUS_NOCHECKPROPTIME & various)) {  

         double tr,tdiff,t,tpropstart;

         tr=M_PI*2.0e6/puls->wr;
         tpropstart=puls->STO_tpropstart_usec[num];
         t=puls->t_usec;

         tdiff=t-tpropstart;

         if (fabs(tdiff-floor(tdiff/tr+0.5)*tr)  > 1e-5) {
           fprintf(stderr,"error: a propagator was calculated at time %.10gusec\nrelative to the"
           " start of the rotor period but reused at time %.10gusec.\n",
            tpropstart-floor(tpropstart/tr)*tr,
            t-floor(t/tr)*tr );
            exit(1);
         }
      }

      puls->t_usec += puls->STO_tproplength_usec[num];
      puls->Uisunit=0;  

      if ( puls->U->col == 1)
        simtrans_diagprop3(puls->sigma,puls->U);
      else
        simtrans(puls->sigma,puls->U);

      ptr++; (puls->curr_nsig)++;
      if (puls->acq_adjoint == 0) {
        *ptr = CRmul(cm_trace(puls->fdetect,puls->sigma),1.0/dnorm);
      } else {
        *ptr = CRmul(cm_trace_adjoint(puls->fdetect,puls->sigma),1.0/dnorm);
      }
  
      if (phase != 0.0) {
          re=ptr->re;
          ptr->re=cosph*re+sinph*ptr->im;
          ptr->im=-sinph*re+cosph*ptr->im;
      }
  }
  _reset_prop();
}

void _delay(double duration)
{
/*  static int first=0; */
  int i,n;
  double dt,t;

  if (duration > puls->dt) {
    duration = puls->dt;
  }
  t=puls->t_usec*1.0e-6;

  /* ZT: changing this strange condition to more proper */
  /*if (puls->H->isdiag && !puls->dor) { */
  if ( (!puls->dor) && (puls->wr > 1.0e-3) && (puls->H->isdiag) ) {/* in case of a diagonal hamiltonian */
      dt=duration;
      ham_rotate2_integrate(puls->H,t+puls->gamma_add/puls->wr*(M_PI/180.0),
                            dt,puls->wr,puls->brl);
      ham_hamilton_integrate(puls->H);
      puls->dU->col = 1;
      prop_realdiag(puls->dU,puls->H->H,dt);
  } else {
      n=(int)ceil(duration/puls->dtmax);
      if (n < 1) n=1;
      dt=duration/(double)n;
      for (i=1;i<=n;i++) {
         if (puls->dor) {
/*           if (!first) {printf("Pulse: %g %g %g %g\n", puls->wr1,puls->brl1,puls->wr2,puls->brl2);first=1;} */
           ham_rotate2_dor(puls->H,puls->wr1*t*(180.0/M_PI),puls->brl1,
	                           puls->wr2*t*(180.0/M_PI),puls->brl2);
	 } else {
           ham_rotate2(puls->H,puls->wr*t*(180.0/M_PI)+puls->gamma_add,puls->brl);
         }
	 ham_hamilton(puls->H);
         if (i == 1) {
	   if (puls->H->isdiag) {
              puls->dU->col = 1;
	      prop_realdiag(puls->dU,puls->H->H,dt);
           } else { 
              puls->dU->col = puls->dU->row;
	      prop_real(puls->dU,puls->H->H,dt); 
	      /*prop_pade_real(puls->dU,puls->H->H,dt); */
           }     
         } else {
           if (puls->H->isdiag)
	      prop_realdiag(puls->tmpU,puls->H->H,dt);
           else
	      prop_real(puls->tmpU,puls->H->H,dt);
	      /*prop_pade_real(puls->tmpU,puls->H->H,dt);    */
              cm_multo_rev(puls->dU,puls->tmpU);
         }
         t += dt;
      }
  }
  update_prop(puls->dU);
  puls->Uisunit=0;  
  puls->t_usec += duration*1.0e6;
}

/* creates the propagators used for the pulses */
int _setrfprop()
{
  int i,j,spin,isany,N,NN;

  isany=0;
  N = puls->N;
  NN = N*N;
  
  memset(puls->sumUph,0,N*sizeof(double));
  memset(puls->sumHrf,0,NN*sizeof(double));

  for (i=1;i<=puls->nchan;i++) {
    if (puls->rfv[i] != 0.0) {
      isany=1;
      if (puls->isselectivepulse) {
        for (j=1;j<=puls->ss->nchanelem[i];j++) {
          spin=puls->ss->chan[i][j];      
          if (puls->spinused[spin]) {
            daxpy_(&N,&(puls->phv[i]),puls->Iz[spin],&INTONE,puls->sumUph,&INTONE);
            daxpy_(&NN,&(puls->rfv[i]),puls->Ix[spin],&INTONE,puls->sumHrf,&INTONE);
          }
        }
      } else {
        daxpy_(&N,&(puls->phv[i]),puls->chan_Iz[i],&INTONE,puls->sumUph,&INTONE);
        daxpy_(&NN,&(puls->rfv[i]),puls->chan_Ix[i],&INTONE,puls->sumHrf,&INTONE);
      }
    }
  }
  /* only use the selection of spins once */
  puls->isselectivepulse=0;
  return isany;
}

void do_Htot(int idealpulse)
{
  double *dm, *dv, *stop;
  int N, NN;
  const double done=1.0;
  
  N = puls->N;
  NN = N*N;
  dm = puls->Htot->data;
  memcpy(dm,puls->sumHrf,NN*sizeof(double));
  if (idealpulse) return;
  
  dv = puls->H->H->data;
  if ( puls->H->isdiag ) {
     /* add to diagonal */
     stop = dv + N;
     do {
        *dm += *dv;
        dv++;
        dm += N+1;
     } while ( dv != stop);
  } else {
     daxpy_(&NN,&done,dv,&INTONE,dm,&INTONE);
  }
}

void _pulse(double duration)
{
  int i,n;
  double dt,t;
  if (!_setrfprop()) {
    /* ZT: relaxation? */
    if (puls->is_relax) {
       _delay_relax(duration);
    } else {
       _delay(duration);
    }
    return;
  }
  if (duration > puls->dt) {
    duration = puls->dt;
  }

  puls->dU->col = puls->N;
  puls->tmpU->col = puls->N;

  n=(int)ceil(duration/puls->dtmax);
  if (n < 1) n=1;
  dt=duration/(double)n;
  t=puls->t_usec*1.0e-6;
  for (i=1;i<=n;i++) {
     if (puls->dor) {
           ham_rotate2_dor(puls->H,puls->wr1*t*(180.0/M_PI),puls->brl1,
	                           puls->wr2*t*(180.0/M_PI),puls->brl2);
     } else {
       ham_rotate2(puls->H,puls->wr*t*(180.0/M_PI)+puls->gamma_add,puls->brl);
     }
     ham_hamilton(puls->H);
     do_Htot(0);
     if (i == 1) {
       prop_real(puls->dU,puls->Htot,dt);
       /*prop_pade_real(puls->dU,puls->Htot,dt);    */
     } else {
       prop_real(puls->tmpU,puls->Htot,dt); 
       /*prop_pade_real(puls->tmpU,puls->Htot,dt);    */
       cm_multo_rev(puls->dU,puls->tmpU);
     }
     t += dt;
  }

  simtrans_zrot2(puls->dU,puls->sumUph);
  update_prop(puls->dU);
  
  puls->Uisunit=0;  
  puls->waspulse=1;  
  puls->t_usec += duration*1.0e6;
}


void _pulseid(double duration)
{
  int i,n;

  if (!_setrfprop()) {
    /* no delay if there is no rf-pulses */
    return;
  }

  if (duration > puls->dt) {
    duration = puls->dt;
  }
  puls->dU->col = puls->N;
  do_Htot(1);
  prop_real(puls->dU,puls->Htot,duration);
  /*prop_pade_real(puls->dU,puls->Htot,duration);*/
  simtrans_zrot2(puls->dU,puls->sumUph);
  update_prop(puls->dU);
 
  puls->Uisunit=0;  
  puls->waspulse=1;  
}

/* ZT: modified to allow this also in relaxation mode */
void _maxdt(double dt) {
  if ( (puls->wr != 0.0) || (puls->is_relax) )
    puls->dtmax= dt;
}


/****
 *  ZT: propagation with simple, static, average hamiltonian - 'avgham_static'
 *      This implementation is as I think it should look like.
 ****/
void _avgham_static(double duration,char* expr)
{
  mv_complx *mtx;

  if (duration > puls->dt) {
    duration = puls->dt;
  }

  mtx = ss_readoper(puls->ss,expr);
  if (!cm_ishermit(mtx)) {
     fprintf(stderr,"avgham_static error: expression doesn't give hermitian operator\n");
     exit(1);
  }
  
  prop_pade_complx(puls->dU,mtx,duration*2.0*M_PI);  /* duration is scaled to take care of Hz->rad.s-1 conversion in Ham. expression */     

  update_prop(puls->dU);
  complx_matrix_free(mtx);
  
  puls->Uisunit=0;  
  puls->waspulse=1;  
  puls->t_usec += duration*1.0e6;
}

/* ZT: about here there was macro 'check_pulse()'. It is moved to pulse.h */

int TclGetPhase(Tcl_Interp* interp,char* string,double* phase)
{

  *phase=0;
  if (!strcasecmp(string,"-X")) {
    *phase=180;
  } else if (!strcasecmp(string,"-Y")) {
    *phase=270;    
  } else if (!strcasecmp(string,"X")) {
    *phase=0;
  } else if (!strcasecmp(string,"Y")) {
    *phase=90;
  } else {
    if (Tcl_GetDouble(interp,string,phase) != TCL_OK)
      return TCL_ERROR;    
  }
  return TCL_OK;
}

int tclAcq(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  double phase;
  int np,prop;
  check_pulse();
  
  if (argc > 4)
    return TclError(interp,"Usage: acq ?phase?\n   or: acq <np> <prop> ?phase?");
  
  phase=0.0;
  if (argc == 2 || argc == 4) {
    if (TclGetPhase(interp,argv[argc-1],&phase) != TCL_OK)
       return TCL_ERROR;
  }
  if (argc > 2) {
    /* disable in relaxation mode */
    if (puls->is_relax)
       return TclError(interp,"acq: using propagators is prohibited in combination with relaxation");

    if (Tcl_GetInt(interp,argv[1],&np) != TCL_OK)
       return TCL_ERROR;
    if (np < 2)
      return TclError(interp,"acq: number of data points must be larger than 1 if a "
                              "propagator is specified");
    if (Tcl_GetInt(interp,argv[2],&prop) != TCL_OK)
       return TCL_ERROR;
    _nacq(np,prop,phase);
  } else {
    _acq(phase);
  }
  return TCL_OK;
}

/***
 * ZT: another version of acq command - with normalization
 ***/
int tclAcq2(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  double phase;
  int np,prop;
  check_pulse();
  
  if (argc > 4)
    return TclError(interp,"Usage: acq2 ?phase?\n   or: acq2 <np> <prop> ?phase?");
  
  phase=0.0;
  if (argc == 2 || argc == 4) {
    if (TclGetPhase(interp,argv[argc-1],&phase) != TCL_OK)
       return TCL_ERROR;
  }
  if (argc > 2) {
    /* disable in relaxation mode */
    if (puls->is_relax)
       return TclError(interp,"acq: using propafators in prohibited in combination with relaxation");

    if (Tcl_GetInt(interp,argv[1],&np) != TCL_OK)
       return TCL_ERROR;
    if (np < 2)
      return TclError(interp,"acq2: number of data points must be larger than 1 if a "
                              "propagator is specified");
    if (Tcl_GetInt(interp,argv[2],&prop) != TCL_OK)
       return TCL_ERROR;
    _nacq2(np,prop,phase);
  } else {
    _acq2(phase);
  }
  return TCL_OK;
}

int tclEvolve(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  check_pulse();
  /* using within pulseq optimizing propagator is prohibited */
  if (OCpar.gradmodeprop) {
     return TclError(interp,"evolve can not be used when optimizing propagator");
  }

  if (argc != 1)
    return TclError(interp,"Usage: evolve");

  /* check if called in gradient mode */
  if (OCpar.gradmode) {
     /* yes, take care of this situation */
     incr_OCmx_pos();
     store_OCprop();
     _evolve_with_prop();
     _reset_prop();
     store_OCdens();
     set_OCmx_code("P");
  } else {
     /* no, do usual things */  
     _evolve_with_prop();
     _reset_prop();
  }
  
  return TCL_OK;
}


int tclOffset(ClientData data,Tcl_Interp* interp,
      int argc, char *argv[])
{
  double offset[MAXCHAN+1],off;
  int i,nchan;

  check_pulse();
  if (argc < 2)
    return TclError(interp,"Usage: offset <offset(1)/Hz> ?<offset(2)/Hz>? ...");

  nchan=argc-1;
  if (nchan != puls->nchan)
    return TclError(interp,"\nerror: offset: arguments must match the number of channels\n");

  if (puls->nchan > MAXCHAN)
    return TclError(interp,"\nerror: offset: number of channels too large -- increase MAXCHAN\n");

  for (i=1;i<=nchan;i++) {
    if (Tcl_GetDouble(interp,argv[i],&off) != TCL_OK)
        return TCL_ERROR;

/*  if (off >= 1.0 && off <= 20.0) {
      fprintf(stderr,"warning: offset called with low value (have you specified wrong arguments)\n");
    } */
    
    offset[i]=off;
  }
  _offset(nchan,offset);
  return TCL_OK;
}

int tclMaxdt(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  double duration;

  check_pulse();

  if (argc != 2)
     return TclError(interp,"Usage: maxdt <duration/usec>");

  if (Tcl_GetDouble(interp,argv[1],&duration) != TCL_OK)
     return TCL_ERROR;

  if (duration < 1e-6)
     return TclError(interp,"maxdt: argument must be larger than 1e-6\n");
  /* ZT: modification to allow maxdt also in relaxation mode */
  /*if (puls->wr == 0.0)
     return TclError(interp,"maxdt: cannot be called if spin-rate is zero\n");
   */
  if ( (puls->wr == 0.0) && (puls->is_relax == 0) )
     return TclError(interp,"maxdt: cannot be called if spin-rate is zero and no relaxation");
     
  _maxdt(duration*1e-6);
  return TCL_OK;
}

int tclReset(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  double tim_usec;
  
  check_pulse();
  if (argc != 1 && argc != 2)
      return TclError(interp,"Usage: reset [time increment/usec]");

  tim_usec=0.0;
  if (argc == 2) {
    if (Tcl_GetDouble(interp,argv[1],&tim_usec) != TCL_OK) return TCL_ERROR;
  }
  _reset(tim_usec);
  
  /* if used within pulseq optimizing propagator set this flag */
  if (OCpar.gradmodeprop) {
     OCpar.propstatus = 0;
  }

  return TCL_OK;
}


int tclPulseid(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  double duration,rf,ph;
  int i,chan;

  check_pulse();
  if (argc < 2 || ((argc % 2) != 0))
    return TclError(interp,"Usage: pulseid <duration/usec> ?<rf(1)/Hz> <phase(1)/degrees>? ?<rf(2)> <phase(2)>? ...");

  if (Tcl_GetDouble(interp,argv[1],&duration) != TCL_OK)
    return TCL_ERROR;

  if (argc > 2) {
    if ((argc/2-1) != puls->nchan)
      return TclError(interp,"pulseid: arguments must match number of channels\n");
 
    for (chan=1,i=2;i<argc;i += 2,chan++) {
      if (Tcl_GetDouble(interp,argv[i],&rf) != TCL_OK) return TCL_ERROR;

      if (rf < 0.0)
        return TclError(interp,"pulseid: rf field must be positive (change pulse phases if you"
                               " want to rotate the other way.");
      _rf(chan,rf);
      if (TclGetPhase(interp,argv[i+1],&ph) != TCL_OK) return TCL_ERROR;

      _ph(chan,ph);            
    }
  }
  if (duration < 0.0)
      return TclError(interp,"pulseid: duration must be zero or positive");
  if (duration != 0.0) 
    _pulseid(duration*1.0e-6);
    
  /* if used within pulseq optimizing propagator set this flag */
  if (OCpar.gradmodeprop) {
     OCpar.propstatus = 1;
  }

  /* do immediate evolution in relaxation mode */
  if (puls->is_relax) {
     _evolve_with_prop();
     _reset_prop();
  }
  
  return TCL_OK;
}

int tclPulse(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  double duration,rf,ph;
  int i,chan;

  check_pulse();
  if (argc < 2 || ((argc % 2) != 0))
    return TclError(interp,"Usage: pulse <duration/usec> ?<rf(1)/Hz> <phase(1)/degrees>? ?<rf(2)> <phase(2)>? ...");

  if (Tcl_GetDouble(interp,argv[1],&duration) != TCL_OK)
     return TCL_ERROR;

  if (argc > 2) {
    if ((argc/2-1) != puls->nchan)
      return TclError(interp,"\nerror: pulse: arguments must match number of channels\n");

    for (chan=1,i=2;i<argc;i += 2,chan++) {
      if (Tcl_GetDouble(interp,argv[i],&rf) != TCL_OK) return TCL_ERROR;

      if (rf < 0.0)
        return TclError(interp,"pulse: rf field must be positive (change pulse phases if you"
                               " want to rotate the other way.");
      _rf(chan,rf);

      if (TclGetPhase(interp,argv[i+1],&ph) != TCL_OK) return TCL_ERROR;
      _ph(chan,ph);            
    }
  }
  if (duration < 0.0)
      return TclError(interp,"pulse: duration must be zero or positive");
  if (duration != 0.0) 
     _pulse(duration*1.0e-6);
     
  /* if used within pulseq optimizing propagator set this flag */
  if (OCpar.gradmodeprop) {
     OCpar.propstatus = 1;
  }

  /* do immediate evolution in relaxation mode */
  if (puls->is_relax) {
     _evolve_with_prop();
     _reset_prop();
  }

  return TCL_OK;
}

/****
 * ZT: helper function for normal pulse_shaped
 ****/
 void _pulse_shaped(int Nelem, int *mask, double steptime)
{
  int i, j;
  
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
  }
  
  /* if used within pulseq optimizing propagator set this flag */
  if (OCpar.gradmodeprop) {
     OCpar.propstatus = 1;
  }
  
}

/****
 * ZT: implementation of shaped pulse, using global variable RFshapes[]
 ****/
int tclPulseShaped(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int i, j, slot, Nelem=-1, Nch=0;
  double duration, steptime;
  int* mask;
  char buf[256], cd[128], buf2[64];;

  check_pulse();

  /* disable when relaxation is ON */
  if (puls->is_relax) {
     fprintf(stderr,"pulse_shaped error: not supported in combination with relaxation.\n");
     exit(1);
  }
  
  mask = int_vector(puls->nchan);
   /* 1-based vector */
   
  /* Debug output */
  /*
  printf("Number of arguments = %d\n", argc);
  for (i=0; i<argc; i++) {
     printf("%d .. %s\n\n",i,argv[i]);
  }
  */
  
  /* READING ARGUMENTS  */
  if (argc <= 2)
    return TclError(interp,"Usage: pulse_shaped <duration/usec> <RFshape on channel 1> ?<RFshape on channel 2>? ...");
  if (argc-2 != puls->nchan) 
    return TclError(interp,"\nerror: pulse_shaped: arguments must match number of channels\n"
                            "        use 'nothing' as place holder when no rf is applied\n");
  if (Tcl_GetDouble(interp,argv[1],&duration) != TCL_OK)
     return TCL_ERROR;
  if (duration < 0.0)
      return TclError(interp,"pulse: duration must be zero or positive");
      
  for (i=1; i<=puls->nchan; i++) {
     if (!strcmp(argv[i+1],"nothing")) {
        mask[i] = -1;
	/* Debug output */
	/*
	printf("Channel %d: no rf applied\n",i);
	*/
	continue;
     } else {
        /* read RFshape and check it */
        if (Tcl_GetInt(interp,argv[i+1],&slot) != TCL_OK) {
	   sprintf(buf,"error in pulse_shaped: argument %d must be interger <RFshape>",i+1);
           return TclError(interp, buf);
        }
	if (!RFshapes[slot]) {
	   sprintf(buf,"pulse_shaped: argument %d points to non-existing RFshape",i+1);
           return TclError(interp,buf);
	}
	
	if (Nelem == -1) {
	    /* set variable for length of RFshapes */
	    Nelem=RFshapes_len(slot);
	} else {
	    /* check consistency of RFshape lengths */
	    if ( RFshapes_len(slot) != Nelem )
	       return TclError(interp,"pulse_shaped: inconsistend number of elements in RFshapes!");
	}
	/* Debug output */
	/*
        printf("Channel %d: Number of elements in the RFshape = %d\n",i-1,Nelem);
	*/
	mask[i] = slot;
     }
   }
   
   /* Debug output - overview of rf parameters 
   for (i=1; i<=puls->nchan; i++) {
     printf("Channel No. %d\n==================\n",i);
     if (mask[i] != -1) {
       for (j=1; j<=Nelem; j++) {
         printf(" [%f %f]\n",RFshapes[mask[i]][j].ampl, RFshapes[mask[i]][j].phase);
       }
     }
   }
   */

   steptime = duration/Nelem*1.0e-6;   
   
   /* setting propagator flag is done only within _pulse_shaped function */

   if (OCpar.gradmode) {
      /* printf("pulse_shaped does gradients\n"); */
      /* determine which channels are active for gradients */
      if (!OCpar.grad_shapes) 
         return TclError(interp,"error when calculating propagators for gradients: grad_shapes not defined");
      int Nsh=LEN(OCpar.grad_shapes);      
      cd[0]='\0';    
      for (i=1; i<=puls->nchan; i++) {
         for (j=1; j<=Nsh; j++) {
	    if (OCpar.grad_shapes[j] == mask[i]) {
	       sprintf(buf2," I%dC%d",j,i);
	       strcat(cd,buf2);
	       Nch++;
	       break;
	     }
	  }
      }  
      /* check if there is any gradient to claculate */
      if (Nch==0) {
         /* no, then do usual things */
	 /* printf("pulse_shaped - no variable shapes, just creates propagator\n"); */
	 _pulse_shaped(Nelem, mask, steptime);
      } else {
         /* yes, check for type of optimization */
	 /* printf("pulse_shaped created code '%s'\n",cd); */
	 if ( OCpar.gradmodeprop == 1 ) {
	    /* do stuff for propagator optimization */
	    _pulse_shapedOCprops(cd,Nch,Nelem,mask,steptime);
	 } else {
	    /* do stuff for state to state optimization */
   	    _pulse_shapedOC(cd,Nch,Nelem,mask,steptime);
	 }
      }

   } else {
      /* do just actual pulsing */
      _pulse_shaped(Nelem, mask, steptime);
   }
   
  free_int_vector(mask);
    
  return TCL_OK;
}




int tclDelay(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  double duration;

  check_pulse();
  if (argc != 2)
      return TclError(interp,"Usage: delay <duration/usec>");

  if (Tcl_GetDouble(interp,argv[1],&duration) != TCL_OK)
      return TCL_ERROR;
  if (duration < 0.0)
      return TclError(interp,"delay: duration must be zero or positive");
  if (duration != 0.0) {
    /* ZT: relaxation? */
    if (puls->is_relax) {
       _delay_relax(duration*1.0e-6);
    } else {
       _delay(duration*1.0e-6);
    }  
    /* if used within pulseq optimizing propagator set this flag */
    if (OCpar.gradmodeprop) {
       OCpar.propstatus = 1;
    }
  }
  
  return TCL_OK;
}

int tclSelect(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int spin,i;
  
  check_pulse();

  if (argc < 2)
    return TclError(interp,"Usage: select <n> <n> ...");

  puls->isselectivepulse=1;
  
  for (i=1;i<=puls->nspins;i++) {
    puls->spinused[i]=0;
  }
  for (i=1;i<argc;i++) {

    if (Tcl_GetInt(interp,argv[i],&spin) != TCL_OK)
      return TCL_ERROR;
    if (spin < 1 || spin > puls->nspins) 
       return TclError(interp,"\nerror: select: spin out of range (%d not between 1 and %d)\n",spin,puls->nspins);    
    puls->spinused[spin]=1;
  }
  return TCL_OK;
}


/****
 * ZT: added more options so the 'store' command can be used to store propagator from a matrix
 ****/
int tclStore(ClientData data,Tcl_Interp* interp,
      int argc, char *argv[])
{
  int num;
  
  check_pulse();

  /* disable when relaxation is ON */
  if (puls->is_relax) {
     fprintf(stderr,"store error: cannot be used in combination with relaxation.\n");
     exit(1);
  }
  
  if ( (argc<2) || (argc>4) ) {
     return TclError(interp,"Usage: store <number>                 \n"
                            "       store <number> <matrix number> \n"
			    "       store <number> list {...}      \n"
			    "       store <number> elements {...}  \n"
			    "       store <number> notelements {...}\n");
  }
  
  if (Tcl_GetInt(interp,argv[1],&num) != TCL_OK)
     return TclError(interp,"store first argument must be integer");
  
  if (argc==2) {
     /* usual, original way */ 
     _store(num);
     return TCL_OK;
  } else if (argc==3) {
     int mxidx;
     
     if (Tcl_GetInt(interp,argv[2],&mxidx) != TCL_OK)
        return TclError(interp,"store second argument must be integer, usage: store <num> <matrix num>");
     if (mxidx<1 || mxidx>MAXMATRIX) 
        return TclError(interp,"store: matrix argument out of range <1,%d>",MAXMATRIX);
     if (!puls->matrix[mxidx])
        return TclError(interp,"store: matrix number %d is undefined",mxidx);
     if ( (puls->matrix[mxidx]->row != puls->N) || (puls->matrix[mxidx]->col != puls->N) ) 
        return TclError(interp,"store: mismatch in matrix dimensions");
	
     if (puls->STO[num] != NULL) complx_matrix_free(puls->STO[num]);

     puls->STO[num] = cmv_dup(puls->matrix[mxidx]);
     
  } else {
     /* case with argc==4 */
     mv_complx *m;
     if ( !strcmp(argv[2],"list") || !strcmp(argv[2],"elements") || !strcmp(argv[2],"notelements") ) {
        m = get_matrix_2(interp,argv[2],argv[3]);
	if (!m) return TclError(interp,"store could not set propagator form %s",argv[2]);
     } else {
        return TclError(interp,"store called with wrong argument %s",argv[2]);
     }
     if ( (m->row != puls->N) || (m->col != puls->N) ) 
        return TclError(interp,"store: mismatch in matrix dimensions");
     /* take care if STO[num] was already alocated */
     if (puls->STO[num] != NULL)
        complx_matrix_free(puls->STO[num]);
     /* assign poiter of already allocated matrix */
     puls->STO[num]=m;
  }
  /* this is done only for cases argc==3 or 4. As a consequence, time checking is not done when reusing this propagator */
  puls->STO_tproplength_usec[num]=0.0;
  puls->STO_tpropstart_usec[num]=-1.0; /*note negative time here */
  puls->STO_waspulse[num]=puls->waspulse;
     
  return TCL_OK;
}


/****
 * ZT: the same as store, but stores adjoint propagator/matrix 
 *     added more options so the 'store' command can be used to store propagator from a matrix
 ****/
int tclStoreAdj(ClientData data,Tcl_Interp* interp,
      int argc, char *argv[])
{
  int num;
  
  check_pulse();

  /* disable when relaxation is ON */
  if (puls->is_relax) {
     fprintf(stderr,"store_adjoint error: cannot be used in combination with relaxation.\n");
     exit(1);
  }
  
  if ( (argc<2) || (argc>4) ) {
     return TclError(interp,"Usage: store_adjoint <number>                 \n"
                            "       store_adjoint <number> <matrix number> \n"
			    "       store_adjoint <number> list {...}      \n"
			    "       store_adjoint <number> elements {...}  \n"
			    "       store_adjoint <number> notelements {...}\n");
  }
  
  if (Tcl_GetInt(interp,argv[1],&num) != TCL_OK)
     return TclError(interp,"store_adjoint first argument must be integer");
  
  if (argc==2) {
     /* usual, original way */ 
     if (puls->cannotbestored) {
       fprintf(stderr,"error: the command 'filter' cannot be stored\n");
       exit(-1);
     }

     if (puls->STO[num] != NULL) complx_matrix_free(puls->STO[num]);
     if (puls->U->col == 1) {
	puls->STO[num] = complx_vector_alloc(puls->N);
	cv_conj(puls->STO[num],puls->U);
     } else {
	puls->STO[num] = complx_matrix_alloc(puls->N,puls->N);
	cm_adjoint(puls->STO[num],puls->U);
     }
     puls->STO_tproplength_usec[num] = puls->t_usec - puls->tpropstart_usec;   
     puls->STO_tpropstart_usec[num] = puls->tpropstart_usec;

     /* save whether a pulse was applied */
     puls->STO_waspulse[num]=puls->waspulse;
     puls->STO_brl[num]=puls->brl;
     return TCL_OK;
  } else if (argc==3) {
     int mxidx;
     
     if (Tcl_GetInt(interp,argv[2],&mxidx) != TCL_OK)
        return TclError(interp,"store second argument must be integer, usage: store_adjoint <num> <matrix num>");
     if (mxidx<1 || mxidx>MAXMATRIX) 
        return TclError(interp,"store_adjoint: matrix argument out of range <1,%d>",MAXMATRIX);
     if (!puls->matrix[mxidx])
        return TclError(interp,"store_adjoint: matrix number %d is undefined",mxidx);
     if ( (puls->matrix[mxidx]->row != puls->N) || (puls->matrix[mxidx]->col != puls->N) ) 
        return TclError(interp,"store_adjoint: mismatch in matrix dimensions");
	
     if (puls->STO[num] != NULL) complx_matrix_free(puls->STO[num]);

     puls->STO[num] = complx_matrix_alloc(puls->N,puls->N);
     cm_adjoint(puls->STO[num],puls->matrix[mxidx]);
     
  } else {
     /* case with argc==4 */
     mv_complx *m;
     if ( !strcmp(argv[2],"list") || !strcmp(argv[2],"elements") || !strcmp(argv[2],"notelements") ) {
        m = get_matrix_2(interp,argv[2],argv[3]);
	if (!m) return TclError(interp,"store_adjoint could not set propagator form %s",argv[2]);
     } else {
        return TclError(interp,"store_adjoint called with wrong argument %s",argv[2]);
     }
     if ( (m->row != puls->N) || (m->col != puls->N) ) 
        return TclError(interp,"store_adjoint: mismatch in matrix dimensions");
     /* check if STO[num] was already alocated */
     if (puls->STO[num] != NULL) complx_matrix_free(puls->STO[num]);
     /* put adjoint to store and free m matrix */
     puls->STO[num] = complx_matrix_alloc(m->col,m->row);
     cm_adjoint(puls->STO[num],m);
     complx_matrix_free(m);
  }
  /* this is done only for cases argc==3 or 4. As a consequence, time checking is not done when reusing this propagator */
  puls->STO_tproplength_usec[num]=0.0;
  puls->STO_tpropstart_usec[num]=-1.0; /*note negative time here */
  puls->STO_waspulse[num]=puls->waspulse;
     
  return TCL_OK;
}




int tclProp(ClientData data,Tcl_Interp* interp,
      int argc, char *argv[])
{
  int num,ntimes;

  check_pulse();
  if (argc != 3 && argc != 2)
      return TclError(interp,"Usage: prop <number> [<times>]");
  if (Tcl_GetInt(interp,argv[1],&num) != TCL_OK)
      return TCL_ERROR;
  if (argc == 2) {
    ntimes=1;
  } else {
    if (Tcl_GetInt(interp,argv[2],&ntimes) != TCL_OK)
        return TCL_ERROR;
  }
  if (ntimes > 0) _prop(num,ntimes);
  
  /* if used within pulseq optimizing propagator set this flag */
  if (OCpar.gradmodeprop) {
     OCpar.propstatus = 1;
  }
  
  return TCL_OK;
}

 
/****
 * ZT: modification to properly deal OC situations
 ****/
int tclFilter(ClientData data,Tcl_Interp* interp,
      int argc, char *argv[])
{
  int num;

  check_pulse();
  /* using within pulseq optimizing propagator is prohibited */
  if (OCpar.gradmodeprop) {
     return TclError(interp,"filter can not be used when optimizing propagator!\n");     
  }
  
  if (argc != 2)
      return TclError(interp,"Usage: filter <number>");
  if (Tcl_GetInt(interp,argv[1],&num) != TCL_OK)
     return TCL_ERROR;
     
  /* check if called in gradient mode */
  if (OCpar.gradmode) {
     /* yes, take care of this situation */
     _filterOC(num);
  } else {
     /* no, do usual things */ 
     _filter(num);
  }
  return TCL_OK;
}



int tclMatrixoper(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int i, j[3];
  mv_complx *m;
  complx val1, val2;
  
  check_pulse();
  
  if (argc != 4 && argc != 7) 
    return TclError(interp,
    "Usage:\n"
    "matrixoper <whichmatrix> <ir> <ic> [(+|-|*|/) <jr> <jc>]\n");

  if (!(m = get_matrix_1(interp, argv[1]))) 
    return TclError(interp, "Specify valid matrix\n");
  
  for (i=1; i<=2; i++) {
    if (Tcl_GetInt(interp, argv[i+1], &j[i]) != TCL_OK)
      return TclError(interp, "Specify valid number\n");
    if (j[i] < 1 || j[i] > puls->N) 
      return TclError(interp, "Invalid range %d\n", j[i]);
  }
  val1 = m->data[j[1]-1+(j[2]-1)*(m->row)];

  if (argc == 7) {
    for (i=1; i<=2; i++) {
      if (Tcl_GetInt(interp, argv[i+4], &j[i]) != TCL_OK)
        return TclError(interp, "Specify valid number\n");
      if (j[i] < 1 || j[i] > puls->N) 
        return TclError(interp, "Invalid range %d\n", j[i]);
    }
    val2 = m->data[j[1]-1+(j[2]-1)*(m->row)];
    switch (*argv[4]) {
    case '+': val1 = Cadd(val1, val2); break;
    case '-': val1 = Csub(val1, val2); break;
    case '*': val1 = Cmul(val1, val2); break;
    case '/': val1 = Cdiv(val1, val2); break;
    default:
      return TclError(interp, "Specify valid operator (+|-|*|/)");
    }
  }
  
  TclAppendResult(interp, "%.8e %.8e", val1.re, val1.im);
  complx_matrix_free(m);
  return TCL_OK;
}

int tclMatrix(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int num;
  mv_complx * m;
  int arg;
  
  check_pulse();
  if (argc < 3 || argc > 5)
      return TclError(interp,
      "Usage:\n"
      "  matrix set <to> <from>\n"
      "  <result> matrix get <from>\n"
      "\n"
      "  <to> can be:\n"
      "       <number>\n"
      "       start\n"
      "       detect\n"
      "\n"
      "  <from> can be:\n"
      "       <number>\n"
      "       start\n"
      "       detect\n"
      "       density\n"
      "       propagator\n"
      "       hamiltonian\n"
      "       avgham\n"
      "       operator {I2x+I3y}\n"
      "       totalcoherence {dm1 dm2 ..}\n"
      "       coherence {{dm1 .. dmN} {dm1 .. dmN}} ..., where N is number of nuclei\n"
      "       list {{c11 c12 ..} {c21 c22 ..} ..}, where c is <re> or {<re> <im>}\n"
      "       elements {{i j} {i j} ..}\n"
      "       notelements {{i j} {i j} ..}\n"
      );
  
  arg=1;
  if (!strcmp(argv[arg],"get")) {
    if (arg + 2 == argc) {
      if ( !(m = get_matrix_1(interp,argv[arg+1])))
        return TCL_ERROR;
        
    } else if (arg + 3 == argc) {
      if ( !(m = get_matrix_2(interp,argv[arg+1],argv[arg+2])))
        return TCL_ERROR;
    } else
      return TclError(interp,"error: 'matrix get' must have one or two additional arguments\n");

    if (TclAppendMatrix(interp,m) == TCL_ERROR)
      return TCL_ERROR;
    complx_matrix_free(m);
    
  } else if (!strcmp(argv[arg],"set")) {
    if (arg + 3 == argc) {
      if ( !(m = get_matrix_1(interp,argv[arg+2])))
        return TCL_ERROR;
    } else if (arg + 4 == argc) {
      if ( !(m = get_matrix_2(interp,argv[arg+2],argv[arg+3])))
        return TCL_ERROR;
    } else
      return TclError(interp,"error: 'matrix set' must have two or three additional arguments\n");
    arg++;
    if (Tcl_GetInt(interp,argv[arg],&num) == TCL_OK) {    
      if (num < 1 || num > MAXMATRIX)
        return TclError(interp,"matrix: called with number larger than %d\n",MAXMATRIX);
      if (puls->matrix[num]) {
        complx_matrix_free(puls->matrix[num]);
      }
      puls->matrix[num]=m;
    } else if (!strcmp(argv[arg],"start")) {
      if (m->row != puls->N || m->col != puls->N) {
        fprintf(stderr,"error: 'matrix set start': matrix %d must be a %dx%d matrix\n",num,puls->N,puls->N);
        exit(1);
      }
      cmv_copy(puls->fstart,m); 
      complx_matrix_free(m);
    } else if (!strcmp(argv[arg],"detect")) {
      if (m->row != puls->N || m->col != puls->N) {
        fprintf(stderr,"error: 'matrix set detect': matrix %d must be a %dx%d matrix\n",num,puls->N,puls->N);
        exit(1);
      }
      cmv_copy(puls->fdetect,m); 
      complx_matrix_free(m);
    } else {
      return TclError(interp,"error: first argument to 'matrix set' must be <number>, 'start' or 'detect', "
                           "but not '%s'\n",argv[arg]);
    
    }
  } else
    return TclError(interp,"error: first argument to 'matrix' must be 'get' or 'set' "
                           "but not '%s'\n",argv[arg]);


  return TCL_OK;
}

int tclTurnon(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int i;
  
  check_pulse();

  if (argc < 2)
    return TclError(interp,"Usage: turnon <int> <int> ...");

  for (i=1;i<argc;i++) {
    ham_turnon(puls->H,argv[i]);
  }
  return TCL_OK;
}

int tclTurnoff(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int i;
  
  check_pulse();

  if (argc < 2)
    return TclError(interp,"Usage: turnoff <int> <int> ...");

  for (i=1;i<argc;i++) {
    ham_turnoff(puls->H,argv[i]);
  }
  return TCL_OK;
}

int tclRotorangle(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  double brl,gamma;
  check_pulse();

  if (argc != 2 && argc != 3)
    return TclError(interp,"Usage: rotorangle (<angle>|reset|ma) ?gamma_add?");

  if (!strcmp(argv[1], "reset")) {
    puls->brl=puls->brlorig;
  } else if (!strcmp(argv[1], "ma")) {
    puls->brl=MAGIC_ANGLE;
  } else {
    if (Tcl_GetDouble(interp,argv[1],&brl) != TCL_OK) return TCL_ERROR;
    puls->brl=brl;
  }
  if (argc == 3) {
    if (Tcl_GetDouble(interp,argv[2],&gamma) != TCL_OK) return TCL_ERROR;
    puls->gamma_add=gamma;
  }
  return TCL_OK;
}

int tclGetInteractions(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  Hamilton* h;
  char buf[256],allbuf[1024];
  int i;
  
  if (argc != 1)
    return TclError(interp,"Usage: {list} getinteractions");

  Tcl_ResetResult(interp);

  h=puls->H;

  *allbuf=0;
  for (i=1;i<=h->n;i++) {
     sprintf(buf,"%s %d",h->n_names[i],h->n_used[i]);
     Tcl_AppendElement(interp,buf);
     strcat(allbuf," ");
     strcat(allbuf,buf);
  }
  
  /* the DC offset is omitted */

  for (i=2;i<=h->nstatic;i++) {
     sprintf(buf,"%s %d",h->nstatic_names[i],h->nstatic_used[i]);
     if (!strstr(allbuf,buf)) {
       Tcl_AppendElement(interp,buf);     
       strcat(allbuf," ");
       strcat(allbuf,buf);
     }
  }

  for (i=1;i<=h->nq;i++) {
     sprintf(buf,"%s %d",h->nq_names[i],h->nq_used[i]);
     if (!strstr(allbuf,buf)) {
       Tcl_AppendElement(interp,buf);     
       strcat(allbuf," ");
       strcat(allbuf,buf);
     }
  }

  for (i=1;i<=h->nqc;i++) {
     sprintf(buf,"%s %d",h->nqc_names[i],h->nqc_used[i]);
     if (!strstr(allbuf,buf)) {
       Tcl_AppendElement(interp,buf);     
       strcat(allbuf," ");
       strcat(allbuf,buf);
     }
  }

  for (i=1;i<=h->nqd;i++) {
     sprintf(buf,"%s %d",h->nqd_names[i],h->nqd_used[i]);
     if (!strstr(allbuf,buf)) {
       Tcl_AppendElement(interp,buf);     
       strcat(allbuf," ");
       strcat(allbuf,buf);
     }
  }

  return TCL_OK;
}   
int tclCurrenttime(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  char buf[128];
  check_pulse();
  
  sprintf(buf,"%g",puls->t_usec);
  Tcl_AppendElement(interp,buf);
  return TCL_OK;
}


/****
 * ZT: Zdenek Tosner - implementation of avgham_static command
 ****/
int tclAvgHam_static(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  double duration;
  char* expr;

  check_pulse();
  if (argc != 3)
    return TclError(interp,"Usage: avgham_static <duration/usec> <expr>");

  if (Tcl_GetDouble(interp,argv[1],&duration) != TCL_OK)
     return TCL_ERROR;

  expr=argv[2];
  
  if (duration < 0.0)
      return TclError(interp,"avgham_static: duration must be zero or positive");
  if (duration != 0.0) 
     _avgham_static(duration*1.0e-6,expr);
     
  /* if used within pulseq optimizing propagator set this flag */
  if (OCpar.gradmodeprop) {
     OCpar.propstatus = 1;
  }
  
  return TCL_OK;
}

/******
 * TV: Routine to return Euler Angles
 ******/

int tclEulerAngles(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  Hamilton* h;
  char buf[256];

  check_pulse();
  if (argc != 1)
    return TclError(interp,"Usage: {alpha beta gamma} eulerangles");
  
  Tcl_ResetResult(interp);

  h=puls->H;
  sprintf(buf, "%g", h->alpha);
  Tcl_AppendElement(interp, buf);
  sprintf(buf, "%g", h->beta);
  Tcl_AppendElement(interp, buf);
  sprintf(buf, "%g", h->gamma);
  Tcl_AppendElement(interp, buf);
  return TCL_OK;
}

/****
 * ZT: implementation of Tcl zgrad_pulse_shape command 
 ****/
int tclZgrPulse(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  double duration, zcoor, offs;
  int slot, Nelem, i, j, used, Nchan, N;
  double *ovals;
  double *mx, *Hst1;

  check_pulse();
  if (argc != 3)
      return TclError(interp,"Usage: zgrad_pulse_shape <duration/usec> <zgrad shape>");

  if (Tcl_GetDouble(interp,argv[1],&duration) != TCL_OK)
      return TclError(interp,"zgrad_pulse_shape: cannot convert argument 1 to duration");
  if (duration < 0.0)
      return TclError(interp,"zgrad_pulse_shape: duration must be zero or positive");
  if (Tcl_GetInt(interp,argv[2],&slot) != TCL_OK)
      return TclError(interp,"zgrad_pulse_shape: cannot convert argument 2 to zgrad_shape");
  if (!ZgradShapes[slot])
      return TclError(interp,"zgrad_pulse_shape: zgrad_shape was not activated/does not exist");

  N = puls->N;
  Hst1 = puls->H->Hstatic[1];
  if (duration != 0.0) {
    Nelem = ZgradShapes_len(slot);
    duration *= 1.0e-6/(double)Nelem;
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
    for (i=1; i<=Nelem; i++) {
       get_chan_offset_ratios(puls->ss, 2.0*M_PI*ZgradShapes[slot][i], ovals);
       for (j=1;j<=Nchan; j++) {
          offs = zcoor*ovals[j];
          daxpy_(&N,&offs,puls->chan_Iz[j],&INTONE,Hst1,&INTONE);
       }
       puls->H->nstatic_used[1] = 1;
       /* ZT: relaxation? */
       if (puls->is_relax) {
          _delay_relax(duration);
       } else {
          _delay(duration);
       }
    }
    ham_set_offset(puls->H,mx,used);
    free((char*)mx);
    free_double_vector(ovals);
    /* if used within pulseq optimizing propagator set this flag */
    if (OCpar.gradmodeprop) {
       OCpar.propstatus = 1;
    }
  }  
  return TCL_OK;
}

/****
 * ZT: helper function for normal pulse_and_zgrad_shaped
 ****/
 void _pulse_and_zgrad_shaped(int Nelem, int *mask, int zgrslot, double steptime)
{
  int i, j, used, Nchan, N;
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
  }
  /* clean up after z grad offsets */
  ham_set_offset(puls->H,mx,used);
  free((char*)mx);
  free_double_vector(ovals);  
  /* if used within pulseq optimizing propagator set this flag */
  if (OCpar.gradmodeprop) {
     OCpar.propstatus = 1;
  }
  
}

/****
 * ZT: implementation of shaped pulse with z gradient on
 *      - using global variables RFshapes[], ZgradShapes[]
 ****/
int tclPulseZgrShaped(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int i, j, slot, zgrslot, Nelem, Nch=0;
  double duration, steptime;
  int* mask;
  char buf[256], cd[128], buf2[4];;

  check_pulse();

  /* disable when relaxation is ON */
  if (puls->is_relax) {
     fprintf(stderr,"pulse_and_zgrad_shaped error: not supported in combination with relaxation.\n");
     exit(1);
  }
  
  mask = int_vector(puls->nchan);

  /* Debug output */
  /*
  printf("Number of arguments = %d\n", argc);
  for (i=0; i<argc; i++) {
     printf("%d .. %s\n\n",i,argv[i]);
  }
  */
  
  /* READING ARGUMENTS  */
  if (argc <= 3)
    return TclError(interp,"Usage: pulse_and_zgrad_shaped <duration/usec> <zgrad shape> <RFshape on channel 1> ?<RFshape on channel 2>? ...");
  if (argc-3 != puls->nchan) 
    return TclError(interp,"\nerror: pulse_and_zgrad_shaped: arguments must match number of channels\n"
                            "        use 'nothing' as place holder when no rf is applied\n");
  if (Tcl_GetDouble(interp,argv[1],&duration) != TCL_OK)
     return TCL_ERROR;
  if (duration < 0.0)
      return TclError(interp,"pulse: duration must be zero or positive");
  if (Tcl_GetInt(interp,argv[2],&zgrslot) != TCL_OK)
     return TCL_ERROR;
  if (!ZgradShapes[zgrslot]) 
     return TclError(interp,"pulse_and_zgrad_shaped: second argument points to non-existing zgrad shape");

  Nelem = ZgradShapes_len(zgrslot);
      
  for (i=1; i<=puls->nchan; i++) {
     if (!strcmp(argv[i+2],"nothing")) {
        mask[i] = -1;
	/* Debug output */
	/*
	printf("Channel %d: no rf applied\n",i);
	*/
	continue;
     } else {
        /* read RFshape and check it */
        if (Tcl_GetInt(interp,argv[i+2],&slot) != TCL_OK) {
	   sprintf(buf,"error in pulse_and_zgrad_shaped: argument %d must be interger <RFshape>",i+2);
           return TclError(interp, buf);
        }
	if (!RFshapes[slot]) {
	   sprintf(buf,"pulse_and_zgrad_shaped: argument %d points to non-existing RFshape",i+2);
           return TclError(interp,buf);
	}
        /* check consistency of RFshape lengths */
        if ( RFshapes_len(slot) != Nelem )
          return TclError(interp,"pulse_and_zgrad_shaped: inconsistend number of elements in RFshapes and zgrad shape!");
	/* Debug output */
	/*
        printf("Channel %d: Number of elements in the RFshape = %d\n",i-1,Nelem);
	*/
	mask[i] = slot;
     }
   }
   
   /* Debug output - overview of rf parameters 
   for (i=1; i<=puls->nchan; i++) {
     printf("Channel No. %d\n==================\n",i);
     if (mask[i] != -1) {
       for (j=1; j<=Nelem; j++) {
         printf(" [%f %f]\n",RFshapes[mask[i]][j].ampl, RFshapes[mask[i]][j].phase);
       }
     }
   }
   */

   steptime = duration/Nelem*1.0e-6;   
   
   /* setting propagator flag is done only within _pulse_shaped function */

   if (OCpar.gradmode) {
      /* printf("pulse_and_zgrad_shaped does gradients\n"); */
      /* determine which channels are active for gradients */
      if (!OCpar.grad_shapes) 
         return TclError(interp,"error when calculating propagators for gradients: grad_shapes not defined");
      int Nsh=LEN(OCpar.grad_shapes);      
      cd[0]='\0';    
      for (i=1; i<=puls->nchan; i++) {
         for (j=1; j<=Nsh; j++) {
	    if (OCpar.grad_shapes[j] == mask[i]) {
	       sprintf(buf2," I%dC%d",j,i);
	       strcat(cd,buf2);
	       Nch++;
	       break;
	     }
	  }
      }  
      /* check if there is any gradient to claculate */
      if (Nch==0) {
         /* no, then do usual things */
	 /* printf("pulse_shaped - no variable shapes, just creates propagator\n"); */
	 _pulse_and_zgrad_shaped(Nelem, mask, zgrslot, steptime);
      } else {
         /* yes, check for type of optimization */
	 /* printf("pulse_shaped created code '%s'\n",cd); */
	 if ( OCpar.gradmodeprop == 1 ) {
	    /* do stuff for propagator optimization */
	    _pulse_and_zgrad_shapedOCprops(cd,Nch,Nelem,mask,zgrslot,steptime);
	 } else {
	    /* do stuff for state to state optimization */
   	    _pulse_and_zgrad_shapedOC(cd,Nch,Nelem,mask,zgrslot,steptime);
	 }
      }

   } else {
      /* do just actual pulsing */
      _pulse_and_zgrad_shaped(Nelem, mask, zgrslot, steptime);
   }
   
  free_int_vector(mask);
    
  return TCL_OK;
}


void tclcmd_pulse(Tcl_Interp* interp) {

  Tcl_CreateCommand(interp,"filter",tclFilter,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"pulse",tclPulse,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"reset",tclReset,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"pulseid",tclPulseid,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"delay",tclDelay,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"maxdt",tclMaxdt,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"acq",tclAcq,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"evolve",tclEvolve,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"prop",tclProp,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

  Tcl_CreateCommand(interp,"store",tclStore,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"offset",tclOffset,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"select",tclSelect,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

  Tcl_CreateCommand(interp,"turnon",tclTurnon,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"turnoff",tclTurnoff,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

  Tcl_CreateCommand(interp,"matrix",tclMatrix,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"matrixoper",tclMatrixoper,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"getinteractions",tclGetInteractions,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

  Tcl_CreateCommand(interp,"currenttime",tclCurrenttime,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"rotorangle",tclRotorangle,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  
/************ ZT: new commands: ***************/
  Tcl_CreateCommand(interp,"acq2",tclAcq2,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"avgham_static",tclAvgHam_static,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"pulse_shaped",tclPulseShaped,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);  
  Tcl_CreateCommand(interp,"store_adjoint",tclStoreAdj,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

/************ TV: new command: ****************/
  Tcl_CreateCommand(interp,"eulerangles",tclEulerAngles,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);  

/**** ZT: commands for zgrad pulses ****/
Tcl_CreateCommand(interp,"zgrad_pulse_shaped",tclZgrPulse,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateCommand(interp,"pulse_and_zgrad_shaped",tclPulseZgrShaped,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

}




