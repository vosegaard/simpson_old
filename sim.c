/*
    Simulation setup and calculation
    Copyright (C) 1999 Mads Bak, Jimmy T. Rasmussen

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
    
    Makes available a setup routine, and a routine for calculation of
    a single crystallite.
    Uses the 'readsys' routine to create the spin-system. Uses the
    'pulse_propagate' routine that performs evaluation of the
    pulses. Uses the functions in 'fidcalc' to perform the 
    evolution of the fid in case of a smart gamma averaged method.
    
    Called by simpson.c that makes the crystallite powder averaging.
*/

#include <stdlib.h>
#include <errno.h>
#include <tcl.h>
#include <string.h>
#include "matrix_new.h"
#include "tclutil.h"
#include "sim.h"
#include "fidcalc.h"
#include "relax.h"
#include "defs.h"
#include "new_direct.h"

/* these definitions moved to sim.h :
#define    M_GAMMAREP        2000
#define    M_GCOMPUTE        2001
#define    M_DIRECT          2002
#define    M_DIRECT_NEW      2009
*/

#define SPINRATE_SMALL 0.001

void readsys(Tcl_Interp* interp,Sim* s);

#define NOTSET -1234.56

/* finds the number of the observable nucleus given the detect operator*/
int obs_nuc(SpinSys* ss, char* det,int* was_all)
{
  char buf[256];
  int i,nuc=0,l,j,nuc2,iso2,iso=0;

  *was_all=0;
  l=strlen(det);
  for (i=0;i<l;i++) {
    if (det[i] == 'I') {
      j=0; i++;
      if (det[i] == 'n') {
        *was_all=1;
      } else {
        while (det[i] >= '0' && det[i] <= '9') { 
          buf[j++]=det[i++];
        }
        buf[j]=0;
        nuc2=strtol(buf,(char**)NULL, 10);
        if (errno) {
          fprintf(stderr,"error: unable to convert '%s' in '%s' to an integer value\n",buf,det);
          //exit(1);
        }
        iso2 = ss_isotope(ss,nuc2)->number;
        if (nuc == 0) {
          nuc=nuc2;
          iso=iso2;
        } else if (iso != iso2)
          return 0;
      }      
    }
  }
  if (*was_all && nuc == 0) return 1;
  return nuc;
}

void sim_initialize(Tcl_Interp* interp,Sim* s)
{
  char detectop[128],startop[256],pulseq[256],buf[256];
  double f;
  int i;

  s->H=(Hamilton*)malloc(sizeof(Hamilton));
  s->ss=(SpinSys*)malloc(sizeof(SpinSys));
  s->P=(Pulse*)malloc(sizeof(Pulse));
  if (!s->H || !s->ss || !s->P) {
    fprintf(stderr,"error: allocation failure in sim_initialize()\n");
    exit(1);
  } 

  f=TclGetDouble(interp,"par","proton_frequency",0,400e6);  
  if (f < 0.0) {
     fprintf(stderr,"error: proton_frequency must be positive\n");
     exit(1);
  }
  if (f != 0 && f <= 10000) {
     fprintf(stderr,"error: proton_frequency must be given in Hz\n");
     exit(1);
  }
  /* Internal in the program we use omega= -gamma B0, that is the
     proton resonance frequency is negative because gamma for at
     proton is positive. */
  s->specfreq = -f;
    
  s->wr=TclGetDouble(interp,"par","spin_rate",0,0);
  if (s->wr < 0.0) {
     fprintf(stderr,"error: spin_rate cannot be negative\n");
     exit(1);
  }
  s->sw=TclGetDouble(interp,"par","sw",0,s->wr);
  s->sw1=TclGetDouble(interp,"par","sw1",0,0);
  s->ni=TclGetInt(interp,"par","ni",0,0);

  TclGetString(interp,s->method,"par","method",0,"direct");
  
  s->dor=TclGetInt(interp,"par","dor",0,0);
  if (s->dor == 1) {
    s->brl1=TclGetDouble(interp,"par","inner_rotor_angle",1,0);
    s->brl2=TclGetDouble(interp,"par","outer_rotor_angle",1,0);
    s->wr1=TclGetDouble(interp,"par","inner_spin_rate",1,0);
    s->wr2=TclGetDouble(interp,"par","outer_spin_rate",1,0);
    s->wr=s->wr1;
  } else {
    s->brl=TclGetDouble(interp,"par","rotor_angle",0,NOTSET);
    if (s->brl == NOTSET) {    
      s->brl=(s->wr == 0.0 ? 0.0 : MAGIC_ANGLE);
      if (verbose & VERBOSE_SIMINFO) {
        printf("The rotor_angle is set to %g\n",s->brl);      
      }
    } 
    s->gamma_zero=TclGetDouble(interp,"par","gamma_zero",0,0);
  }
  /* s->ngamma=TclGetInt(interp,"par","gamma_angles",0,1); */
/* AB::beg */
/* Check, whether the flag "use_3_angle_set" is set.                   */
/* If a 3 angle set is used "ngamma" is set to 1.                      */
/* This doesn't allow the use of gamma compute and                     */
/* 3 angle sets at the same time. This might be not the ideal          */
/* choice. In principle a 3 angle set could be used with gamma compute */
/* leading to an additional (useless) averaging over gamma.            */
  s->use_3_angle_set=TclGetInt(interp,"par","use_3_angle_set",0,0);
  if( s->use_3_angle_set == 1 ) {
    s->ngamma=1;
  } else {
    s->ngamma=TclGetInt(interp,"par","gamma_angles",0,1);
  }
/* AB::end */

  s->P->acq_adjoint = TclGetInt(interp,"par","acq_adjoint",0,0);

  TclGetString(interp,buf,"par","dipole_check",0,"true");
  if (!strcmp(buf,"true")) {
    s->dipole_check=1;
  } else if (!strcmp(buf,"false")) {
    s->dipole_check=0;
  } else {
    fprintf(stderr,"error: 'dipole_check' must be 'false' or 'true'\n");
    exit(1);
  }

  s->np=TclGetInt(interp,"par","np",0,1);
  s->ntot=s->np*(s->ni > 1 ? s->ni : 1);

  s->nsampr_rep=TclGetInt(interp,"par","fixed_rep",0,0);
  /* 0 is default value for option f[na-j]=conj(f[j]) */
  s->realspec=TclGetInt(interp,"par","real_spec",0,0);
  s->blkdiag=TclGetVector(interp,"par","block_diag",0,NULL);
  TclGetString(interp,detectop,"par","detect_operator",0,"Inp");

  TclGetString(interp,s->crystfile,"par","crystal_file",0,"rep1");
  TclGetString(interp,s->rfproffile,"par","rfprof_file",0,""); 
  TclGetString(interp,startop,"par","start_operator",0,"Inz");
  TclGetString(interp,pulseq,"par","pulse_sequence",0,"pulseq");
  readsys(interp,s);
  s->matdim=s->ss->matdim;
  s->fdetect=ss_readoper(s->ss,detectop);
  if (s->fdetect == NULL) {
      fprintf(stderr,"error: unable to interpret detect operator '%s'\n",detectop);
      exit(1);
  }
  m_makeQlist(s->fdetect);
//printf("sim ini Qlist done\n");
  /* ZT: modification, start_operator can be also 'equilibrium' */
  if (!strncmp(startop,"equil",5)) {
     s->fstart = ss_eq_sigma(s->ss);
  } else {
     s->fstart=ss_readoper(s->ss,startop);
  }
  if (s->fstart == NULL) {
      fprintf(stderr,"error: unable to interpret start operator '%s'\n",startop);
      exit(1);
  }

  /* conjugate_fid: This parameter is used overrule the automatic detection of the magnetogyric
    ratio for the observable nucleus. It can be set to 'auto', 'true' or 'false.
  */
  TclGetString(interp,buf,"par","conjugate_fid",0,"auto");
  if (!strcmp(buf,"auto")) {
    int was_all;
    s->obs_nuc = obs_nuc(s->ss,detectop,&was_all);
    if ( (s->obs_nuc < 1) || (was_all && s->ss->nchan > 1) ) {
      fprintf(stderr,"error: in detect-operator '%s'. More than one type of nucleus\n"
                     "was detected. Set 'conjugate_fid' to 'true' or 'false' to turn off\n"
                     "this sanity check. 'conjugate_fid' is 'auto' per default which conjugates\n"
                     "the fid if the gyromagnetic ratio is larger than zero for the observable\n" 
                     "nucleus. That corrects for the standard axis convention of plotting spectra\n"
                     "that ignores the sign of gamma.\n"
                     ,detectop);
      exit(1);
    }
    s->conjugate_fid=  (ss_gamma(s->ss,s->obs_nuc) > 0 ? 1 : 0);
  } else if (!strcmp(buf,"true")) {
    s->conjugate_fid = 1;
  } else if (!strcmp(buf,"false")) {
    s->conjugate_fid = 0;
  } else {
    fprintf(stderr,"error: 'conjugate_fid' must be 'true', 'false' or 'auto'");
    exit(1);
  }

  if (verbose & VERBOSE_OPER) {
    printf("Acquisition data will%s be complex conjugated.\n",
       (s->conjugate_fid ? "" : " not"));
  }

  s->rhosymmetry=is_rhosymmetry(s->fstart,s->fdetect);
  if (verbose & VERBOSE_OPER) {
    cm_print(s->fstart,"Start operator");
    cm_print(s->fdetect,"Detect operator");
    printf("Rho-symmetry: %s\n",(s->rhosymmetry ? "yes" : "no"));
  }

  if (!strcmp(s->method,"gammarep")){
    s->imethod=M_GAMMAREP;
  } else if (!strcmp(s->method,"gcompute")){
    s->imethod=M_GCOMPUTE;
  } else if ( !strcmp(s->method,"direct")){
    s->imethod=M_DIRECT;
  } else if ( !strcmp(s->method,"idirect")){
    s->imethod=M_DIRECT_NEW;
  } else if ( !strcmp(s->method,"gcompute_new")){
    s->imethod=M_GCOMPUTE_NEW;
  } else if ( !strcmp(s->method,"igcompute")){
    s->imethod=M_GCOMPUTE2_NEW;
  } else {
    fprintf(stderr,"error: method '%s' not known\n",s->method);
    fprintf(stderr,"       must be one of : gammarep, gcompute, direct, idirect, igcompute\n");
    exit(1);
  }

  if ( (s->imethod == M_DIRECT) || (s->imethod == M_DIRECT_NEW) ) {
    s->nstepr=1;

    if (s->wr == 0.0) {
      if (s->ngamma != 1) {
        fprintf(stderr,"error: 'gamma_angles' must be one when calculating static spectra\n");
        exit(1);
      }
    }

  } else if ( (s->imethod == M_GCOMPUTE_NEW) || (s->imethod == M_GCOMPUTE2_NEW) ) {
    
    if (s->dor==1) {
      fprintf(stderr,"error: 'igcompute' method not compatible with calculating DOR spectra\n");
      exit(1);
    }

    if (s->wr < SPINRATE_SMALL) {
      fprintf(stderr,"error: must use the 'direct' or 'idirect' method when calculating static spectra (spin-rate=0).\n");
      exit(1);
    }
    
  } else {
    double sw0=s->sw;

    if (s->dor==1) {
      fprintf(stderr,"error: must use the 'direct' method when calculating DOR spectra\n");
      exit(1);
    }

    if (s->wr < SPINRATE_SMALL) {
      fprintf(stderr,"error: must use the 'direct' method when calculating static spectra (spin-rate=0).\n");
      exit(1);
    }

    s->nstepr= floor(s->sw/s->wr+0.5);
    s->sw = (double)(s->nstepr)*s->wr;
    if (fabs(sw0-s->sw) > 0.0001) {
      fprintf(stderr,"error: the spectralwidth must be adjusted from %g to %g to fulfill\n"
                     "       the condition: spectral_width/spin_rate == gamma_angles\n",sw0,s->sw);
      exit(1);
    }
    if (s->sw == 0.0) {
      fprintf(stderr,"error: spectral_width is set to zero\n");
      exit(1);
    }
    /* There are two possibilities for ngamma and nstepr to be multiples of each other
       and they demands different calculation procedures. */
    if (s->ngamma >= s->nstepr) {
      s->gammethod=1;
      s->ng=s->ngamma/s->nstepr;  /* ngamma reflects the real number of gamma
                                     angles, that is an multiplum of nstepr. */
      if (s->ng < 1) s->ng=1;
      if (s->ngamma != s->ng*s->nstepr) {
         s->ngamma=s->ng*s->nstepr;
      }
    } else {
      s->gammethod=2;
      s->ns=s->nstepr/s->ngamma;
      if (s->ns < 1) s->ns=1;
      if (s->nstepr != s->ns*s->ngamma) {
        /*  printf("warning: nstepr changed to %d\n",s->ns*s->ngamma); */
        s->nstepr=s->ns*s->ngamma;
      }
    }
    if (s->nsampr_rep == 0 && s->imethod == M_GAMMAREP) {
      s->nsampr_rep=rep_minimize_estimate(s->ntot,s->nstepr,s->matdim,m_getnqlist());
      if (verbose & VERBOSE_SIMINFO) printf("Replication factor (method 1) found to: %d\n",s->nsampr_rep);
      s->nsampr_rep=rep_minimize(s->matdim,s->nstepr,s->ntot);
      if (verbose & VERBOSE_SIMINFO) printf("Replication factor (method 2) found to: %d\n",s->nsampr_rep);
    }

    s->tstepr=1.0/((double)s->nstepr*s->wr);

  }
  s->wr *= M_PI*2.0;
  s->wr1 *= M_PI*2.0;
  s->wr2 *= M_PI*2.0;

  pulse_initialize(s->P,interp,s->ss,s->sw,s->wr,s->brl,s->dor,s->wr1,s->brl1,s->wr2,s->brl2);
  pulse_setpulsename(s->P,pulseq);

/*
  if (s->imethod == M_DIRECT) {
      pulse_setprop(s->P,s->fstart,s->fdetect);
  }
*/  
  pulse_setprop(s->P,s->fstart,s->fdetect);
  if (s->blkdiag != NULL) {
    if (verbose & VERBOSE_SIMINFO) printf("Preparing for block diagonal matrices\n");
    /* m_setmatrixtype_block(s->blkdiag,s->matdim); */
    fprintf(stderr,"WARNING: Block-diagonal matrices not implemented, option ignored\n");
  } else {
    /*  m_setmatrixtype_normal();  */
  }
  s->fid = complx_vector(s->ntot); /* old style */
  m_zerov(s->fid);
  
  if ( (s->imethod == M_GAMMAREP) || (s->imethod == M_GCOMPUTE) ) {
     s->un=(mv_complx**)malloc(sizeof(mv_complx*)*(s->nstepr+1));
     s->vn=(mv_complx**)malloc(sizeof(mv_complx*)*(s->nstepr+1));
     if (!s->un || !s->vn) {
       fprintf(stderr,"error: allocation failure (vector matrix 1)\n");
       exit(1);
     }
     for (i=0;i<=s->nstepr;i++) {
       s->un[i] = complx_matrix_alloc(s->matdim,s->matdim);
       s->vn[i] = complx_matrix_alloc(s->matdim,s->matdim);
     }
     s->vi = complx_matrix_alloc(s->matdim,s->matdim);
     s->vj = complx_matrix_alloc(s->matdim,s->matdim);
     s->ku = complx_matrix_alloc(s->matdim,s->matdim);
  } else {
     s->un = NULL;
     s->vn = NULL;
     s->vi = complx_matrix_alloc(s->matdim,s->matdim);
     s->vj = NULL;
     s->ku = NULL;
  }
  
  /* ZT: relaxation setup */
  char rxbuf[32];
  TclGetString(interp,rxbuf,"par","relax",0,"off");
  /* printf("RX: relax is %s\n",rxbuf); */
  if (!strncmp(rxbuf,"on",2)) {
     s->P->is_relax = 1;
     readsys_relax(interp);
     read_relax(interp,s);
  } else {
     s->P->is_relax = 0;
  }
  
  /* ZT: method for propagator calculation */
  TclGetString(interp,buf,"par","prop_method",0,"diag");
  if (!strncmp(buf,"diag",4)) {
     s->P->propmethod = 0;
  } else if (!strncmp(buf,"pade_real",9)) {
     s->P->propmethod = 1;
  } else if (!strncmp(buf,"pade_compl",10)) {
     s->P->propmethod = 2;
  } else if (!strncmp(buf,"cheby_ratio",10)) {
     s->P->propmethod = 3;
  } else if (!strncmp(buf,"cheby_poly",10)) {
     s->P->propmethod = 4;
  } else if (!strncmp(buf,"taylor",6)) {
     s->P->propmethod = 5;
  } else {
    fprintf(stderr,"error: prop_method '%s' not known\n",buf);
    fprintf(stderr,"       must be one of : diag, pade_real, pade_complex, cheby_rational\n");
    fprintf(stderr,"                        cheby_polynomial, taylor\n");
    exit(1);
  }
  if (verbose & VERBOSE_SIMINFO) printf("prop_method is %i\n",s->P->propmethod);
  
  /* ZT: new tricks with pulseq */
  if ( (s->imethod == M_DIRECT_NEW) || (s->imethod == M_GCOMPUTE_NEW) || (s->imethod == M_GCOMPUTE2_NEW) )
     new_direct_initialize(interp,s);
}

void sim_destroy(Sim* s)
{
  int i;

  /* ZT: relaxation setup */
  if (s->P->is_relax) destroy_Relax(s->ss->nspins);

  /* ZT: new tricks with pulseq */
  if ( (s->imethod == M_DIRECT_NEW) || (s->imethod == M_GCOMPUTE_NEW) || (s->imethod == M_GCOMPUTE2_NEW) )
     new_direct_destroy(s);
  
  pulse_destroy(s->P);
  ham_destroy(s->H);
  m_destroyQlist();
  free(s->H);
  free(s->P);
  free(s->ss);
  if (s->blkdiag) {
    /* m_free_dv(s->blkdiag); */
    /* ZT: hopefully never used */
  }
  complx_matrix_free(s->fdetect);
  complx_matrix_free(s->fstart);
  free_complx_vector(s->fid);
  
  if ( (s->imethod == M_GAMMAREP) || (s->imethod == M_GCOMPUTE) ) {
    for (i=0;i<=s->nstepr;i++) {
      if (s->un[i]) complx_matrix_free(s->un[i]);
      if (s->vn[i]) complx_matrix_free(s->vn[i]);
    }
    free(s->un);
    free(s->vn);
    complx_matrix_free(s->vj);
    complx_matrix_free(s->ku);
  }
  complx_matrix_free(s->vi);
  
}

/***
 * ZT: modified to include rf profile scalling factors, includes changes by AB
 ***/
int sim_calcfid(Sim* s,Omega omega,double* rfscalefact,complx* fidsum)
{
  int i,ig,nfid;
  double t;
  mv_complx *wfid, *wfidsum;
  /* old style wrappers */
  wfid = (mv_complx*)malloc(sizeof(mv_complx));
  wfid->row = LEN(s->fid); wfid->col = 1; wfid->data = &(s->fid[1]);
  wfidsum = (mv_complx*)malloc(sizeof(mv_complx));
  wfidsum->row = LEN(fidsum); wfidsum->col = 1; wfidsum->data = &(fidsum[1]);
  
  /* ZT: assign values of rf scale factors to structure Pulse */
  for (i=1;i<=s->ss->nchan;i++) {
     s->P->rfscalefactors[i] = rfscalefact[i];
  } 
  
  nfid=0;
  ham_rotate(s->H,omega.alpha, omega.beta, s->gamma_zero+omega.gamma);
  s->H->alpha = omega.alpha;
  s->H->beta = omega.beta;
  m_zerov(fidsum);
  if (s->imethod == M_GAMMAREP || s->imethod == M_GCOMPUTE) {
    int k,nstepr,ng;
    double tstepr;

    /* ZT: disable these when relaxation is invoked */
    if ( s->P->is_relax ) {
       fprintf(stderr,"error: methods gammarep and gcompute not allowed with relaxation\n");
       exit(1);
    }
    
    nstepr=s->nstepr;
    tstepr=s->tstepr;
    ng=s->ng;

    if (s->gammethod != 1) {
      fprintf(stderr,"error: gamma must be larger than nstepr for this method\n");
      exit(1);
    }
    for (ig=0;ig < ng;ig++) {
      for (k=1;k<=nstepr;k++) {
        t=(k-1)*tstepr + 2.0*M_PI*(double)ig/(double)(ng*nstepr)/s->wr;	      
        pulse_propagate(s->P,s->H,s->vi,tstepr,t,(complx*)NULL);
        cmv_copy(s->un[k],s->vi);
        if (k == 1) {
	   cmv_copy(s->vn[k],s->vi);
        } else {
	   cm_mul(s->vn[k],s->vi,s->vn[k-1]);
	}
      }
      if (s->imethod == M_GAMMAREP) {
        gammarep(s->vn,s->un,nstepr,s->fstart,s->fdetect,tstepr,s->fid,s->nsampr_rep);
      } else {
        /* ZT: I guess it should be like this and not as in old code where
	       gcompute was called after possible call of gammarep   */
	s->rhosymmetry=is_rhosymmetry(s->P->fstart,s->P->fdetect);
        gcompute(s->vn,nstepr,s->P->fstart,s->P->fdetect,tstepr,s->fid,s->rhosymmetry,s->realspec);
      }
      cmv_addto(wfidsum,wfid);
      nfid++;
    }
  } else if (s->imethod == M_DIRECT) {
    for (ig=1;ig<=s->ngamma;ig++) {

      if (s->wr == 0.0) {
        t=0;
      } else {
	      t = (ig-1)/(double)s->ngamma*2.0*M_PI/s->wr;
      }
      s->H->gamma = 360*(ig-1)/(double)s->ngamma + (s->gamma_zero+omega.gamma);
      pulse_propagate(s->P,s->H,s->vi,1.0e99,t,s->fid);
      cmv_addto(wfidsum,wfid);
      nfid++;
    }
  } else if (s->imethod == M_DIRECT_NEW) {
      new_direct_calc(s,s->gamma_zero+omega.gamma, wfidsum, wfid, &nfid);
  } else if (s->imethod == M_GCOMPUTE_NEW) {
      /* ZT: disable these when relaxation is invoked */
      if ( s->P->is_relax ) {
         fprintf(stderr,"error: methods gammarep, gcompute and igcompute not allowed with relaxation\n");
         exit(1);
      }
      new_gcompute_calc(s,s->gamma_zero, wfidsum, wfid, &nfid);
  } else if (s->imethod == M_GCOMPUTE2_NEW) {
      /* ZT: disable these when relaxation is invoked */
      if ( s->P->is_relax ) {
         fprintf(stderr,"error: methods gammarep, gcompute and igcompute not allowed with relaxation\n");
         exit(1);
      }
      new_gcompute2_calc(s,s->gamma_zero, wfidsum, wfid, &nfid);
  }
  
  cmv_muld(wfidsum,1.0/nfid);
  if (s->conjugate_fid) {
    cv_conj_in(wfidsum);
  }

  free((char*)wfid);
  free((char*)wfidsum);


  return 0;
}


