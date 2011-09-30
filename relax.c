/*
    Routines connected with relaxation
    Copyright (C) 2008 Zdenek Tosner

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
    
    Reads the spin system from the 'spinsys' section in the input
    file and creates the Hamiltonian for the spin system.
    
     Plenty os functions called from different parts of the code :-)
*/

#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <tcl.h>
#include "matrix_new.h"
#include "spinsys.h"
#include "ham.h"
#include "sim.h"
#include "pulse.h"
#include "defs.h"
#include "relax.h"

/****  Defining global variable Relax ****/
RXstruct Relax;

/****
 * Reading spinsysres again and extracting info relevant for Relax
 ****/
void readsys_relax(Tcl_Interp *interp) 
{
   int Nnuc, i, j, N, NN, n1, n2;
   char *buf, **nucnames, **elname;
   char bufstr[1024];
   Tcl_Obj **elptr, **elptr2;
   double d;

   /* get number of nuclei and their names */
   buf=Tcl_GetVar2(interp,"spinsysres","nuclei",TCL_GLOBAL_ONLY);
   if (!buf) {
      fprintf(stderr,"readsys_relax: error reading spinsysres(nuclei)");
      exit(1);
   }
   if (Tcl_SplitList(interp, buf, &Nnuc, &nucnames) != TCL_OK) {
      fprintf(stderr,"readsys_relax error: cannot decompose list spinsysres(nuclei)\n");
      exit(1);
   }
   
   printf("readsys_relax is creating interaction maps for %d nuclei\n",Nnuc);
   
   /* allocate interaction maps */
   Relax.shifts = (double**)malloc((Nnuc+1)*sizeof(double*));
   Relax.quadrupoles = (double**)malloc((Nnuc+1)*sizeof(double*));
   Relax.dipoles = (double***)malloc((Nnuc+1)*sizeof(double**));
   Relax.dipoles[1] = (double**)malloc((Nnuc*Nnuc+1)*sizeof(double*));
      for (i=2; i<=Nnuc; i++) Relax.dipoles[i]=Relax.dipoles[i-1]+Nnuc;
   for (i=1; i<=Nnuc; i++) {
      Relax.shifts[i]=NULL;
      Relax.quadrupoles[i]=NULL;
      for (j=1; j<=Nnuc; j++) {
	 Relax.dipoles[i][j] = NULL;
      }
   }
   
   /* fill these maps */
   if (Tcl_Eval(interp,"array names spinsysres") != TCL_OK) {
      fprintf(stderr,"readsys_relax: %s\n",interp->result);
      exit(1);
   }
   strcpy(bufstr,interp->result);
   if (Tcl_SplitList(interp, bufstr, &N, &elname) != TCL_OK) {
      fprintf(stderr,"readsys_relax error: cannot decompose spinsysres array names from a list\n");
      exit(1);
   }

   for (i=0; i<N; i++) {
      if (!strncmp(elname[i],"nuclei",6)) continue;
      if (!strncmp(elname[i],"channels",8)) continue;
      if (!strncmp(elname[i],"mixing",6)) continue;

      if (Tcl_ListObjGetElements(interp,Tcl_GetVar2Ex(interp,"spinsysres",elname[i],TCL_GLOBAL_ONLY), &NN, &elptr2) != TCL_OK) {
         fprintf(stderr,"readsys_relax error: cannot decompose list of parameters for spinsysres(%s)\n",elname[i]);
	 exit(1);
      }
      if (!strncmp(elname[i],"shift",5)) {
         if ( NN != 7) {
	    fprintf(stderr,"readsys_relax error in reading shift - parameter count mismatch\n");
	    exit(1);
	 }
	 if (Tcl_GetIntFromObj(interp,elptr2[0],&n1) != TCL_OK) {
	    fprintf(stderr,"readsys_relax cannot get int for shift (nucleus)\n");
	    exit(1);
	 }
	 if ( (n1<1) || (n1>Nnuc) ) {
	    fprintf(stderr,"readsys_relax error: shift declaration for nonexisting nucleus %d\n",n1);
	    exit(1);
	 }
	 if (Relax.shifts[n1] != NULL) {
	    fprintf(stderr,"readsys_relax error: duplicate shift declaration in spinsysres\n");
	    exit(1);
	 }
	 printf("shifts[%d] =",n1);
	 Relax.shifts[n1] = (double*)malloc(6*sizeof(double));
	 for (j=1; j<NN; j++) {
            if (Tcl_GetDoubleFromObj(interp,elptr2[j],&d) != TCL_OK) {
	       fprintf(stderr,"readsys_relax error in converting shift values to doubles\n");
	       exit(1);
	    }
	    Relax.shifts[n1][j-1] = d;
	    printf(" %g;",d);
	 }
	 printf("\n");
      }
      if (!strncmp(elname[i],"quadrupole",10)) {
         if ( NN != 7) {
	    fprintf(stderr,"readsys_relax error in reading quadrupole - parameter count mismatch\n");
	    exit(1);
	 }
	 if (Tcl_GetIntFromObj(interp,elptr2[0],&n1) != TCL_OK) {
	    fprintf(stderr,"readsys_relax cannot get int for quadrupole (nucleus)\n");
	    exit(1);
	 }
	 if ( (n1<1) || (n1>Nnuc) ) {
	    fprintf(stderr,"readsys_relax error: quadrupole declaration for nonexisting nucleus %d\n",n1);
	    exit(1);
	 }
	 if (Relax.quadrupoles[n1] != NULL) {
	    fprintf(stderr,"readsys_relax error: duplicate quadrupole declaration in spinsysres\n");
	    exit(1);
	 }
	 printf("quadrupoles[%d] =",n1);
	 Relax.quadrupoles[n1] = (double*)malloc(6*sizeof(double));
	 for (j=1; j<NN; j++) {
            if (Tcl_GetDoubleFromObj(interp,elptr2[j],&d) != TCL_OK) {
	       fprintf(stderr,"readsys_relax error in converting quadrupole values to doubles\n");
	       exit(1);
	    }
	    Relax.quadrupoles[n1][j-1] = d;
	    printf(" %g;",d);
	 }
	 printf("\n");
      }
      if (!strncmp(elname[i],"dipole",6)) {
         if ( NN != 6) {
	    fprintf(stderr,"readsys_relax error in reading dipole - parameter count mismatch\n");
	    exit(1);
	 }
	 if (Tcl_GetIntFromObj(interp,elptr2[0],&n1) != TCL_OK) {
	    fprintf(stderr,"readsys_relax cannot get int for dipole (1. nucleus)\n");
	    exit(1);
	 }
	 if (Tcl_GetIntFromObj(interp,elptr2[1],&n2) != TCL_OK) {
	    fprintf(stderr,"readsys_relax cannot get int for dipole (2. nucleus)\n");
	    exit(1);
	 }
	 if ( (n1<1) || (n1>Nnuc) || (n2<1) || (n2>Nnuc) || (n1==n2) ) {
	    fprintf(stderr,"readsys_relax error: wrong combination of nuclei (%d,%d) in dipole declaration\n",n1,n2);
	    exit(1);
	 }
	 if (n1>n2) {
	    int dum;
	    dum = n2;
	    n2 = n1;
	    n1 = dum;
	 }
	 if (Relax.dipoles[n1][n2] != NULL) {
	    fprintf(stderr,"readsys_relax error: duplicate dipole declaration in spinsysres\n");
	    exit(1);
	 }
	 printf("dipoles[%d][%d] =",n1,n2);
	 Relax.dipoles[n1][n2] = (double*)malloc(4*sizeof(double));
	 for (j=2; j<NN; j++) {
            if (Tcl_GetDoubleFromObj(interp,elptr2[j],&d) != TCL_OK) {
	       fprintf(stderr,"readsys_relax error in converting dipole values to doubles\n");
	       exit(1);
	    }
	    Relax.dipoles[n1][n2][j-2] = d;
	    printf(" %g;",d);
	 }
	 printf("\n");
      }
      Tcl_Free((char *) elptr2);
   }

   Tcl_Free((char *) nucnames);
   printf("readsys_relax is DONE.\n");
   
}

/****
 * Encoding spectral density of particular relaxation mechanism
 ****/
int get_Jtype(char* glm,char *jlocname)
{
   int type;
   
   if (!strncmp(glm,"non",3)) {
      type = 0;
   } else if (!strncmp(glm,"sph",3)) {
      type = 1;
   } else if (!strncmp(glm,"sym",3)) {
      type = 2;
   } else if (!strncmp(glm,"asy",3)) {
      type = 3;
   } else {
      fprintf(stderr,"get_Jtype error: unknown type of glob. motion '%s'\n",glm);
      exit(1);
   }
   if (!strncmp(jlocname,"rig",3)) {
      type += 0;
   } else if (!strncmp(jlocname,"model_free_ext",14)) {
      type += 20;
   } else if (!strncmp(jlocname,"model_free",10)) {
      type += 10;
   } else if (!strncmp(jlocname,"diffusion_on",12)) {
      type += 30;
   } else if (!strncmp(jlocname,"diffusion_in",12)) {
      type += 40;
   } else if (!strncmp(jlocname,"3_sites",7)) {
      type += 50;
   } else {
      fprintf(stderr,"get_Jtype error: unknown type of local motion '%s'\n",jlocname);
      exit(1);
   }
   if (type==0) {
      fprintf(stderr,"get_Jtype error: combination of global an local motion is not allowed\n");
      exit(1);
   }
   /* LIMIT TO LIPARI-SZABO DURING DEVELOPMENT */
   if (type != 11) {
      fprintf(stderr,"get_Jtype error: USE L-S SETTINGS DURING DEVELOPMENT! \n");
      exit(1);
   }
   return type;   
}



/****
 * Reading relax Tcl variable and initial set up of C global Relax structure
 ****/
void read_relax(Tcl_Interp *interp, Sim* s) 
{
   int N, i, j, k, n1, n2, lam, matdim, Nj, Nglm, jpars;
   char **elname, *jlocname, *glmname;
   char bufstr[256];
   Tcl_Obj *o1, *glm, *jlist;
   Tcl_Obj **elptr, **jelptr, **glmptr;
   double om1, om2, field, cseta, dum;
   SpinSys *ss;
   SecCouple *scdum;
   mv_complx *mp1, *mp2, *mm1, *mm2, *mz1, *mz2, *mx;
   
   printf("read_relax  starts now\n");
   
   ss = s->ss;
   matdim = ss->matdim;
   field = fabs((s->specfreq)/ss_gamma1H())*2.0*M_PI; /* this is in 10^-7 Tesla */
   printf("   Magnetic field is %f\n",field*1.0e-7);
   Relax.eqnorm = 1.0/(2.0*M_PI*fabs(s->specfreq));
   printf("   Relax.eqnorm is %g\n",Relax.eqnorm);
   
   /* get type of global motion into Tcl_Obj pointer */
   if ((glm=Tcl_GetVar2Ex(interp,"relax","global_motion",TCL_GLOBAL_ONLY)) == NULL) {
      fprintf(stderr,"read_relax: error reading relax(global_motion)");
      exit(1);
   }
   if (Tcl_ListObjGetElements(interp, glm, &Nglm, &glmptr) != TCL_OK) {
      fprintf(stderr,"read_relax error: cannot decompose parameters from a list for global_motion\n");
      exit(1);
   }
   glmname = Tcl_GetString(glmptr[0]);
   if (!glmname) {
      fprintf(stderr,"read_relax error: cannot convert glm to string\n");
      exit(1);
   }
   printf("   global motion is %s\n",glmname);
   
   /* get number of relaxation mechanisms */
   if ((o1=Tcl_GetVar2Ex(interp,"relax","auto",TCL_GLOBAL_ONLY)) == NULL) {
      fprintf(stderr,"read_relax: error reading relax(auto)");
      exit(1);
   }
   if (Tcl_GetIntFromObj(interp,o1,&N) != TCL_OK) {
      fprintf(stderr,"read_relax: cannot get int from relax(auto)");
      exit(1);
   }
   Relax.Nauto = N;
   printf("   Number of relaxation mechanisms is %d\n",N);
   
   if ((o1=Tcl_GetVar2Ex(interp,"relax","cross",TCL_GLOBAL_ONLY)) == NULL) {
      fprintf(stderr,"read_relax: error reading relax(cross)");
      exit(1);
   }
   if (Tcl_GetIntFromObj(interp,o1,&N) != TCL_OK) {
      fprintf(stderr,"read_relax: cannot get int from relax(cross)");
      exit(1);
   }
   Relax.Ncross = N;
   printf("   Number of cross correlations is %d\n",N);

   /* allocate rows in Q array */
   Relax.Q = (mv_complx***)malloc(Relax.Nauto*sizeof(mv_complx**));
   if (!Relax.Q) {
      fprintf(stderr,"read_relax error: cannot allocate rows in Q\n");
      exit(1);
   }
   /* allocate rows in Y array */
   Relax.Y = (mv_complx***)malloc(Relax.Nauto*sizeof(mv_complx**));
   if (!Relax.Y) {
      fprintf(stderr,"read_relax error: cannot allocate rows in Y\n");
      exit(1);
   }
   /* allocate rows in NQ vector */
   Relax.NQ = (int*)malloc(Relax.Nauto*sizeof(int));
   if (!Relax.NQ) {
      fprintf(stderr,"read_relax error: cannot allocate rows in NQ\n");
      exit(1);
   }
   /* allocate rows in omega array */
   Relax.omega = (double**)malloc(Relax.Nauto*sizeof(double*));
   if (!Relax.omega) {
      fprintf(stderr,"read_relax error: cannot allocate rows in omega\n");
      exit(1);
   }
   /* allocate rows in Jcode array */
   Relax.Jcode = (SpecDens*)malloc(Relax.Nauto*sizeof(SpecDens));
   if (!Relax.Jcode) {
      fprintf(stderr,"read_relax error: cannot allocate rows in Jcode\n");
      exit(1);
   }
   /* allocate rows in secular array */
   Relax.secular = (SecCouple**)malloc(Relax.Nauto*sizeof(SecCouple*));
   if (!Relax.secular) {
      fprintf(stderr,"read_relax error: cannot allocate rows in secular\n");
      exit(1);
   }

   /* get elements in relax array */
   if (Tcl_Eval(interp,"array names relax") != TCL_OK) {
      fprintf(stderr,"read_relax: %s\n",interp->result);
      exit(1);
   }
   strcpy(bufstr,interp->result);

   if (Tcl_SplitList(interp, bufstr, &N, &elname) != TCL_OK) {
      fprintf(stderr,"relax_read error: cannot decompose array names from a list\n");
      exit(1);
   }

   /* go through the elements and create all necessary info */
   lam = -1;
   for (i=0; i<N; i++) {
      printf("   inspecting element %s\n",elname[i]);
      
      if (!strncmp(elname[i],"auto",4)) continue;
      if (!strncmp(elname[i],"cross",5)) continue;
      if (!strncmp(elname[i],"global_motion",13)) continue;
      
      /* dipole dipole relaxation */
      if (!strncmp(elname[i],"dipole",6)) {
	 if (sscanf(elname[i],"dipole_%d_%d",&n1,&n2) != 2) {
	    fprintf(stderr,"read_relax error: cannot get nuclei numbers from %s\n",elname[i]);
	    exit(1);
	 }
         printf("read_relax found dipole between nuclei %d and %d\n", n1, n2);
	 if (n1>n2) {
	    int dum;
	    dum = n2;
	    n2 = n1;
	    n1 = dum;
	 }
	 if ( Relax.dipoles[n1][n2] == NULL ) {
	    fprintf(stderr,"read_relax error: dipole %d %d not defined in spinsys\n",n1, n2);
	    exit(1);
	 }
	 lam++;
	 printf("   it has index lam = %d\n",lam);
	 /* allocate 9 matrix representations of Q, fill them */
	 Relax.Q[lam] = (mv_complx**)malloc(9*sizeof(mv_complx*));
	 mp1 = Ip(ss,n1);
	 mp2 = Ip(ss,n2);
	 mm1 = Im(ss,n1);
	 mm2 = Im(ss,n2);
	 mz1 = Iz(ss,n1);
	 mz2 = Iz(ss,n2);
	 mx = complx_matrix_alloc(matdim,matdim);
	 cm_mul(mx,mm1,mm2);
	 cmv_muld(mx,0.5);
	 Relax.Q[lam][0] = cmv_dup(mx);
	 cm_mul(mx,mm1,mz2);
	 cmv_muld(mx,0.5);
	 Relax.Q[lam][1] = cmv_dup(mx);
	 cm_mul(mx,mz1,mm2);
	 cmv_muld(mx,0.5);
	 Relax.Q[lam][2] = cmv_dup(mx);
	 cm_mul(mx,mz1,mz2);
	 cmv_muld(mx,2.0/sqrt(6.0));
	 Relax.Q[lam][3] = cmv_dup(mx);
	 cm_mul(mx,mp1,mm2);
	 cmv_muld(mx,-0.5/sqrt(6.0));
	 Relax.Q[lam][4] = cmv_dup(mx);
	 cm_mul(mx,mm1,mp2);
	 cmv_muld(mx,-0.5/sqrt(6.0));
	 Relax.Q[lam][5] = cmv_dup(mx);
	 cm_mul(mx,mp1,mz2);
	 cmv_muld(mx,-0.5);
	 Relax.Q[lam][6] = cmv_dup(mx);
	 cm_mul(mx,mz1,mp2);
	 cmv_muld(mx,-0.5);
	 Relax.Q[lam][7] = cmv_dup(mx);
	 cm_mul(mx,mp1,mp2);
	 cmv_muld(mx,0.5);
	 Relax.Q[lam][8] = cmv_dup(mx);
	 complx_matrix_free(mx);
	 complx_matrix_free(mp1);
	 complx_matrix_free(mp2);
	 complx_matrix_free(mm1);
	 complx_matrix_free(mm2);
	 complx_matrix_free(mz1);
	 complx_matrix_free(mz2);
	 printf("   Relax.Q[%d][0-8] filled.\n",lam);
	 /* allocate space for 9 Y matrices */
	 Relax.Y[lam] = (mv_complx**)malloc(9*sizeof(mv_complx*));
	 for (j=0; j<9; j++) {
            Relax.Y[lam][j] = NULL;
         }
	 /* set number of allocated Q matrices */
	 Relax.NQ[lam] = 9;
         /* allocate and fill frequencies */
	 Relax.omega[lam] = (double*)malloc(9*sizeof(double));
	 om1 = ss_gamma(s->ss,n1)*field; /* this is now in rad.s-1 */
	 om2 = ss_gamma(s->ss,n2)*field;
	 Relax.omega[lam][0] = -om1-om2;
	 Relax.omega[lam][1] = -om1;
	 Relax.omega[lam][2] = -om2;
	 Relax.omega[lam][3] = 0.0;
	 Relax.omega[lam][4] = om1-om2;
	 Relax.omega[lam][5] = -om1+om2;
	 Relax.omega[lam][6] = om1;
	 Relax.omega[lam][7] = om2;
	 Relax.omega[lam][8] = om1+om2;
	 printf("   Relax.omega[%d][0-8] filled.\n",lam);
	 /* allocate and fill secular couples */
	 if (ss_issame(ss,n1,n2)) {
	    /* homonuclear */
	    scdum = (SecCouple*)malloc(20*sizeof(SecCouple));
	    *(int*)(scdum) = 19;
	    for (k=1; k<=9; k++) {
	       scdum[k].q = k-1;
	       scdum[k].y = k-1;
	    }
            scdum[10].q = 1; scdum[10].y = 2;
            scdum[11].q = 2; scdum[11].y = 1;
            scdum[12].q = 3; scdum[12].y = 4;
            scdum[13].q = 4; scdum[13].y = 3;
            scdum[14].q = 3; scdum[14].y = 5;
            scdum[15].q = 5; scdum[15].y = 3;
            scdum[16].q = 4; scdum[16].y = 5;
            scdum[17].q = 5; scdum[17].y = 4;
            scdum[18].q = 6; scdum[18].y = 7;
            scdum[19].q = 7; scdum[19].y = 6;
	    printf("   homonuclear scdum[1-19] filled");
	 } else {
	    /* heteronuclear */
            scdum = (SecCouple*)malloc(10*sizeof(SecCouple));
	    *(int*)scdum = 9;
	    for (k=1; k<=9; k++) {
	       scdum[k].q = k-1;
	       scdum[k].y = k-1;
	    }
	    printf("   heteronuclear scdum[1-9] filled");
	 }
	 Relax.secular[lam] = scdum;
	 printf(" and finaly assigned to Relax.secular[%d].\n",lam);
	 printf(" Cross check: %d is %d\n",*(int*)scdum,*(int*)(Relax.secular[lam]));
	 /* decide on spectral density */
	 if ((jlist=Tcl_GetVar2Ex(interp,"relax",elname[i],TCL_GLOBAL_ONLY)) == NULL) {
            fprintf(stderr,"read_relax: error reading relax(%s)",elname[i]);
            exit(1);
         }
	 if (Tcl_ListObjGetElements(interp, jlist, &Nj, &jelptr) != TCL_OK) {
            fprintf(stderr,"read_relax error: cannot decompose parameters from a list for %s\n",elname[i]);
            exit(1);
         }
         jlocname = Tcl_GetString(jelptr[0]);
         if (!jlocname) {
            fprintf(stderr,"relax_read error: cannot convert local motion to string for %s\n", elname[i]);
            exit(1);
         }
 	 Relax.Jcode[lam].type = get_Jtype(glmname,jlocname);
         printf("   Relax.Jcode[%d] asigned to %d",lam,Relax.Jcode[lam].type);
	 jpars = Nj+Nglm-2+1;
	 printf(", pars will have %d elems:\n",jpars);
	 Relax.Jcode[lam].par = double_vector(jpars);
	 /* set interaction constant */
	 dum = (Relax.dipoles[n1][n2][0])*2.0*M_PI;
	 Relax.Jcode[lam].par[1] = 0.5*dum*dum*6.0;
         printf("   [%g;",Relax.Jcode[lam].par[1]);
	 /* set other parameters */
	 j = 2;
	 for (k=1; k<Nglm; k++) {
	    if (Tcl_GetDoubleFromObj(interp,glmptr[k],&(Relax.Jcode[lam].par[j])) != TCL_OK) {
	       fprintf(stderr,"read_relax error: (glm) cannot get J parameter as double (%s)\n",elname[i]);
	       exit(1);
	    }
	    printf(" %g;",Relax.Jcode[lam].par[j]);
	    j++;
	 }
	 for (k=1; k<Nj; k++) {
	    if (Tcl_GetDoubleFromObj(interp,jelptr[k],&(Relax.Jcode[lam].par[j])) != TCL_OK) {
	       fprintf(stderr,"read_relax error: (loc) cannot get J parameter as double (%s)\n",elname[i]);
	       exit(1);
	    }
	    printf(" %g;",Relax.Jcode[lam].par[j]);
	    j++;
	 }
	 printf("\n");
	 Tcl_Free((char *) jelptr);
	 continue;
      }
      
      /* chemical shift anizotropy relaxation */
      if (!strncmp(elname[i],"shift",5)) {
	 if (sscanf(elname[i],"shift_%d",&n1) != 1) {
	    fprintf(stderr,"read_relax error: cannot get nucleus number from %s\n",elname[i]);
	    exit(1);
	 }
         printf("read_relax found shift for nucleus %d\n",n1);
	 if ( Relax.shifts[n1] == NULL ) {
	    fprintf(stderr,"read_relax error: shift %d not defined in spinsys\n",n1);
	    exit(1);
	 }
	 lam++;
	 printf("   it has index lam = %d\n",lam);
	 /* allocate 3 matrix representations of Q, fill them */
	 Relax.Q[lam] = (mv_complx**)malloc(3*sizeof(mv_complx*));
	 Relax.Q[lam][0] = Im(ss,n1);
	 cmv_muld(Relax.Q[lam][0],0.5);
	 Relax.Q[lam][1] = Iz(ss,n1);
	 cmv_muld(Relax.Q[lam][1],2.0/sqrt(6.0));
	 Relax.Q[lam][2] = Ip(ss,n1);
	 cmv_muld(Relax.Q[lam][2],-0.5);
	 printf("   Relax.Q[%d][0-2] filled.\n",lam);
	 /* allocate space for 3 Y matrices */
	 Relax.Y[lam] = (mv_complx**)malloc(3*sizeof(mv_complx*));
	 for (j=0; j<3; j++) {
            Relax.Y[lam][j] = NULL;
         }
	 /* set number of allocated Q matrices */
	 Relax.NQ[lam] = 3;
         /* allocate and fill frequencies */
	 Relax.omega[lam] = (double*)malloc(3*sizeof(double));
	 om1 = ss_gamma(s->ss,n1)*field; /* this is in rad.s-1 */
	 Relax.omega[lam][0] = -om1;
	 Relax.omega[lam][1] = 0.0;
	 Relax.omega[lam][2] = om1;
	 printf("   Relax.omega[%d][0-2] filled.\n",lam);
	 /* allocate and fill secular couples */
	 Relax.secular[lam] = (SecCouple*)malloc(4*sizeof(SecCouple));
	 *(int*)(Relax.secular[lam]) = 3;
	 Relax.secular[lam][1].q = 0;
	 Relax.secular[lam][1].y = 0;
	 Relax.secular[lam][2].q = 1;
	 Relax.secular[lam][2].y = 1;
	 Relax.secular[lam][3].q = 2;
	 Relax.secular[lam][3].y = 2;
	 printf("   Relax.secular[%d][1-3] filled\n",lam);
	 /* decide on spectral density */
	 if ((jlist=Tcl_GetVar2Ex(interp,"relax",elname[i],TCL_GLOBAL_ONLY)) == NULL) {
            fprintf(stderr,"read_relax: error reading relax(%s)",elname[i]);
            exit(1);
         }
	 if (Tcl_ListObjGetElements(interp, jlist, &Nj, &jelptr) != TCL_OK) {
            fprintf(stderr,"read_relax error: cannot decompose parameters from a list for %s\n",elname[i]);
            exit(1);
         }
         jlocname = Tcl_GetString(jelptr[0]);
         if (!jlocname) {
            fprintf(stderr,"get_Jtype error: cannot convert local motion to string for %s\n", elname[i]);
            exit(1);
         }
 	 Relax.Jcode[lam].type = get_Jtype(glmname,jlocname);
         printf("   Relax.Jcode[%d] asigned to %d",lam,Relax.Jcode[lam].type);
	 jpars = Nj+Nglm-2+1;
	 printf(", pars will have %d elems:\n",jpars);
	 Relax.Jcode[lam].par = double_vector(jpars);
	 /* set interaction constant */
	 dum = (Relax.shifts[n1][1])*2.0*M_PI;
	 cseta = Relax.shifts[n1][2];
	 Relax.Jcode[lam].par[1] = 0.5*dum*dum*(1.0+cseta*cseta/3.0);
         printf("   [%g;",Relax.Jcode[lam].par[1]);
	 /* set other parameters */
	 j = 2;
	 for (k=1; k<Nglm; k++) {
	    if (Tcl_GetDoubleFromObj(interp,glmptr[k],&(Relax.Jcode[lam].par[j])) != TCL_OK) {
	       fprintf(stderr,"read_relax error: (glm) cannot get J parameter as double (%s)\n",elname[i]);
	       exit(1);
	    }
	    printf(" %g;",Relax.Jcode[lam].par[j]);
	    j++;
	 }
	 for (k=1; k<Nj; k++) {
	    if (Tcl_GetDoubleFromObj(interp,jelptr[k],&(Relax.Jcode[lam].par[j])) != TCL_OK) {
	       fprintf(stderr,"read_relax error: (loc) cannot get J parameter as double (%s)\n",elname[i]);
	       exit(1);
	    }
	    printf(" %g;",Relax.Jcode[lam].par[j]);
	    j++;
	 }
	 printf("\n");
	 Tcl_Free((char *) jelptr);
	 continue;
      }
      
      /* quadrupolar relaxation */
      if (!strncmp(elname[i],"quadrupole",10)) {
	 if (sscanf(elname[i],"quadrupole_%d",&n1) != 1) {
	    fprintf(stderr,"read_relax error: cannot get nucleus number from %s\n",elname[i]);
	    exit(1);
	 }
         printf("read_relax found quadrupole for nucleus %d\n", n1);
	 if ( Relax.quadrupoles[n1] == NULL ) {
	    fprintf(stderr,"read_relax error: quadrupole %d not defined in spinsys\n",n1);
	    exit(1);
	 }
	 /* NOT IMPLEMENTED FURTHER !!! */
	 fprintf(stderr,"read_relax error: quadrupole not implemented yet\n");
	 exit(1);
	 lam++;
      }
      
      /* random field relaxation */
      if (!strncmp(elname[i],"random_field",12)) {
	 if (sscanf(elname[i],"random_field_%d",&n1) != 1) {
	    fprintf(stderr,"read_relax error: cannot get nucleus number from %s\n",elname[i]);
	    exit(1);
	 }
         printf("read_relax found random_field for nucleus %d\n", n1);
	 /* NOT IMPLEMENTED FURTHER !!! */
	 fprintf(stderr,"read_relax error: random field not implemented yet\n");
	 exit(1);
	 lam++;
      }
      
      /* cross correlations not implemented yet */
      if (!strncmp(elname[i],"cross_correlation",6)) {
	 /* NOT IMPLEMENTED FURTHER !!! */
         printf("read_relax found cross_correlation\n");
	 fprintf(stderr,"read_relax: cross correlations not implemented yet\n");
	 exit(1);
      }

   }
   printf("read_relax is DONE.\n");
}


void destroy_Relax(int Nnuc)
{
   int i, j;
   int N, NN, qmx, ymx;
   
   printf("Destroy Relax\n=============\n");
   printf("Nauto = %d, Ncross = %d\n",Relax.Nauto, Relax.Ncross);
   if (!(Relax.secular)) {
     printf("  problem, Relax.secular does not exist\n");
     exit(1);
   } else {
     printf("  OK 1\n");
   }
      
   for (i=0; i<Relax.Nauto; i++) {
      if (!(Relax.secular[i])) {
         printf("  problem, Relax.secular[%d] does not exist\n",i);
	 exit(1);
      } else {
         printf("  OK 2\n");
      }
      N = *(int*)(Relax.secular[i]);
      printf("  OK 3, N = %d\n",N);
      qmx = 0;
      ymx = 0;
      for (j=1; j<=N; j++) {
         if (Relax.secular[i][j].q > qmx) qmx = Relax.secular[i][j].q;
	 if (Relax.secular[i][j].y > ymx) ymx = Relax.secular[i][j].y;
	 /* printf("  OKej %d\n",j); */
      }

      printf("Mech. %d : %d secular terms, Qmax = %d, Ymax = %d\n",i,N,qmx,ymx);
      printf("           free Q[%d]:",i);
      
      for (j=0; j<=qmx; j++) {
         if (Relax.Q[i][j]) {
            complx_matrix_free(Relax.Q[i][j]);
	    printf("(%d,%d)",i,j);
	 }
      }
      free((char *)(Relax.Q[i]));
      printf(" done.\n");
      
      if (Relax.Y[i]) {
         printf("           free Y[%d]:",i);
         for (j=0; j<=ymx; j++) {
            if (Relax.Y[i][j]) {
               complx_matrix_free(Relax.Y[i][j]);
	       printf(" (%d,%d)",i,j);
	    }
         }
         free((char *)(Relax.Y[i]));
	 printf(" done.\n");
      } else {
         printf("           Y[%d] not present\n",i);
      }
      
      printf("           Jtype = %d, par = ",Relax.Jcode[i].type);
      NN = LEN(Relax.Jcode[i].par);
      for (j=1; j<=NN; j++) {
         printf(" %g;",Relax.Jcode[i].par[j]);
      }
      printf("\n           freeing Jcode[%d].par",i);
      free_double_vector(Relax.Jcode[i].par);
      printf(" done.\n");
      
      printf("           secular =");
      for (j=1; j<=N; j++) {
         printf(" (%d,%d)",Relax.secular[i][j].q, Relax.secular[i][j].y);
      }
      printf("\n           freeing secular[%d]",i);
      free((char *)(Relax.secular[i]));
      printf(" done.\n");

      printf("           omega =");
      for (j=0; j<=qmx; j++) {
         printf(" %g;",Relax.omega[i][j]);
      }
      printf("\n           freeing omega[%d]",i);
      free((char *)(Relax.omega[i]));
      printf(" done.\n");
      
   }
   printf("Freeing Jcode");
   free((char *)(Relax.Jcode));
   printf(" done.\n");
   
   printf("Freeing secular");
   free((char *)(Relax.secular));
   printf(" done.\n");
   
   printf("Freeing omega");
   free((char *)(Relax.omega));
   printf(" done.\n");
   
   printf("Freeing Q");
   free((char *)(Relax.Q));
   printf(" done.\n");
   
   if (Relax.Y)
   printf("Freeing omega[%d]",i);
   free((char *)(Relax.Y));
   printf(" done.\n");

   printf("Interaction maps:\n");
   for (i=1; i<=Nnuc; i++) {
      if (Relax.shifts[i]) {
         printf("   freeing shifts[%d]",i);
         free((char *)(Relax.shifts[i]));
	 printf(" done.\n");
      }
      if (Relax.quadrupoles[i]) {
         printf("   freeing quadrupoles[%d]",i);
         free((char *)(Relax.quadrupoles[i]));
	 printf(" done.\n");
      }
      for (j=1; j<=Nnuc; j++) {
         if (Relax.dipoles[i][j]) {
            printf("   freeing dipoles[%d][%d]",i,j);
            free((char *)(Relax.dipoles[i][j]));
	    printf(" done.\n");
         }
      }
   }
   printf("Freeing dipoles[1]");
   free((char *)(Relax.dipoles[1]));
   printf(" done.\n");

   printf("Freeing dipoles");
   free((char *)(Relax.dipoles));
   printf(" done.\n");

   printf("Freeing quadrupoles");
   free((char *)(Relax.quadrupoles));
   printf(" done.\n");

   printf("Freeing shifts");
   free((char *)(Relax.shifts));
   printf(" done.\n");
   
}

/****
 * Spectral density functions, based on key Relax.Jcode
 ****/
double spectral_density(int lam, double om)
{
   double res, kk;
   int N, type;
   
   if (!(Relax.Jcode)) {
      fprintf(stderr,"spectral_density error: Relax.Jcode does not exist.\n");
      exit(1);
   } 
   if (!(Relax.Jcode[lam].par)) {
      fprintf(stderr,"spectral_density error: Relax.Jcode[%d].par does not exist.\n",lam);
      exit(1);
   }
   
   type = Relax.Jcode[lam].type;
   N = LEN(Relax.Jcode[lam].par);
   
   switch (type)
     {
      case 11:
       { /* Lipari and Szabo */
	 double tauM, taue, SS;
	 if (N != 4) {
	    fprintf(stderr,"spectral_density error: not enough parameters in Relax.Jcode[%d].par for type %d (%d, not 4)\n",lam,type,N);
	    exit(1);
	 }
	 kk = Relax.Jcode[lam].par[1];
	 tauM = 1.0/(6.0*Relax.Jcode[lam].par[2]);
	 taue = 1.0/tauM + 1.0/(Relax.Jcode[lam].par[4]);
	 taue = 1.0/taue;
	 SS = Relax.Jcode[lam].par[3];
	 res = (SS*tauM)/(1.0+om*om*tauM*tauM)+((1.0-SS)*taue)/(1+om*om*taue*taue);
	 res = 2.0/5.0*res*kk;
         break;
	}
      default:
         {
	  fprintf(stderr,"spectral_density error: unknown type %d\n",type);
	  exit(1);
	 }
     }
     
   return res;
}

/****
 * helper function for maxtrix multiplication complex * real
 * NOT TESTED, NOT USED!!!
 ****/
void m_mulCoRe(complx **res, complx **a, double **b)
{

}

/*****
 * Make puls visible also here
 *****/
 extern Pulse* puls;


/****
 * Here we go with relaxation evolution...
 ****/
void _delay_relax(double duration)
{
  int i,j,k,l,n,q,y, Ndim, Nsec, isHdiag;
  double dt, omq, om, omsc, J;
  mv_complx *Q, *Y, *mx, *mx1, *mx2;
  mv_double *Ham;

  if (puls->wr != 0) {
    fprintf(stderr,"error: relaxation not supported for sample rotation\n");
    exit(1);
  }
  
  isHdiag = puls->H->isdiag;
  printf("DELAY: The hamiltonian is%s diagonal\n",(isHdiag ? "" : " not"));
  n=(int)ceil(duration/puls->dtmax);
  if (n < 1) n=1;
  dt=duration/(double)n;

  /* get current Hamiltonian */
  Ham = ham_hamilton(puls->H);
  /* IMPORTANT: Ham is a copy of pointer puls->H->H and shouldn't be freed ! */
  Ndim = puls->N;
  mx = complx_matrix_alloc(Ndim,Ndim);
  mx1 = complx_matrix_alloc(Ndim,Ndim);
  mx2 = complx_matrix_alloc(Ndim,Ndim);
  
  if (isHdiag) {
        /* in case of a diagonal hamiltonian things are simpler */
	for (i=0; i<Relax.Nauto; i++) {
	   Nsec = Relax.NQ[i];
	   /* printf("interaction %d, number of elements %d\n",i,Nsec); */
	   for (j=0; j<Nsec; j++) {
	      /* printf("   element %d\n",j); */
	      if (!Relax.Y[i][j]) {
	         Relax.Y[i][j] = complx_matrix_alloc(Ndim,Ndim);
		 /* printf("          -> allocating now\n"); */
	      }
	      Y = Relax.Y[i][j];
	      Q = Relax.Q[i][j];
	      /* get frequency */
	      omq = Relax.omega[i][j];
	      /* prepare Y matrix */
	      for (k=0; k<Ndim; k++) {
                 /* printf("     Ham[%d][%d] re = %g; im = %g\n",k,k,Ham[k][k].re,Ham[k][k].im);
		 */
	         for (l=0; l<Ndim; l++) {
		    /* printf("     Ham[%d][%d] re = %g; im = %g\n",l,l,Ham[l][l].re,Ham[l][l].im)
		    */;
		    om = omq + Ham->data[k] - Ham->data[l];
		    /* printf("     om = %g",om); */
		    J = spectral_density(i,om);
		    /* printf(";   J = %g\n",J); */
		    Y->data[k+l*Ndim].re = (Q->data[l+k*Ndim].re)*J;
		    Y->data[k+l*Ndim].im = -(Q->data[l+k*Ndim].im)*J;
		 }
	      }
	      /* m_print(Relax.Q[i][j],"Que");
	      m_print(Relax.Y[i][j],"Yps"); */
	   }
        }
        /* tohle je asi hotovy... */
  } else {
         /* diagonalize Hamiltonian, retrieve eigenvalues and transformation matrix */
         mv_double *r, *eig;
	 mv_complx *rc;
	 complx dum;
	 int ii, jj;
	 double pomr;

         r = double_matrix_alloc(Ndim,Ndim);
         eig = double_vector_alloc(Ndim);
	 rc = complx_matrix_alloc(Ndim,Ndim);
         dm_symdiag(Ham,eig,r);
	 dmv_complx(r,rc);
	 /* m_print(Ham,"Hamiltonian");
	 m_print(rc,"eigenvectors");
	 printf("eigenvalues: "); 
	 for (i=1; i<=Ndim; i++) printf(" %g",eig[i]);
	 printf("\n"); */
 
         for (i=0; i<Relax.Nauto; i++) {
	    Nsec = Relax.NQ[i];
	    /* printf("interaction %d, number of elements %d\n",i,Nsec); */
	    for (j=0; j<Nsec; j++) {
	       /* printf("   element %d\n",j); */
	       if (!Relax.Y[i][j]) {
	          Relax.Y[i][j] = complx_matrix_alloc(Ndim,Ndim);
		  /* printf("          -> allocating now\n"); */
	       }
	       Y = Relax.Y[i][j];
	       Q = Relax.Q[i][j];
	       omq = Relax.omega[i][j];
	       /* transform Q to new basis */
               cmv_copy(mx,Q);
	       simtrans_adj(mx,rc);
	       /* m_print(Relax.Q[i][j],"Que");
	       m_print(mx,"Que transformed"); */
	       /* prepare matrix Y */
	       cmv_zero(Y);
	       /* printf("      J = "); */
	       for (k=0; k<Ndim; k++) {
	          for (l=0; l<Ndim; l++) {
		     dum = mx->data[k+Ndim*l];
		     if ( (fabs(dum.re)<1e-8) && (fabs(dum.im)<1e-8) )continue; 
		     om = omq + eig->data[k] - eig->data[l];
		     J = spectral_density(i,om);
		     /* printf(" %g",J); */
		     for (ii=0; ii<Ndim; ii++) {
		        for (jj=0; jj<Ndim; jj++) {
			   pomr = r->data[ii+k*Ndim]*r->data[jj+l*Ndim];
			   /* if (pomr == 0.0) continue; */
		           mx2->data[jj+ii*Ndim].re = pomr*dum.re*J;
			   mx2->data[jj+ii*Ndim].im = -pomr*dum.im*J;
			 }
		     }
		     /* m_print(mx2,"dilci matice Y"); */
		     cmv_addto(Y,mx2);
		  }
	       }
	       /* printf("\n");
	       m_print(Relax.Y[i][j],"Yps"); */
	   }
	}
	double_matrix_free(r);
	double_vector_free(eig);
	complx_matrix_free(rc);
  }
  
  for (i=1; i<=n; i++) {
     cmv_copy(mx,puls->sigma);
    /* printf("   evolving increment %d\n",i);
     m_print(mx,"sig"); */
     /* evolve with Hamiltonian */
     if (isHdiag) {
     /* printf("kkk 1\n");*/
        puls->dU->col = 1;
	prop_realdiag(puls->dU,Ham,dt);
    /* printf("kkk 2\n"); */
	simtrans_diagprop3(puls->sigma,puls->dU);
    /* printf("kkk 3\n"); */
     } else {
        puls->dU->col = Ndim;
        /* prop_real(puls->dU,Ham,dt); */
	prop_pade_real(puls->dU,Ham,dt);
        simtrans(puls->sigma,puls->dU);
     }
    /* printf("kkk 4, Relax.Nauto = %d\n",Relax.Nauto); */
     /* evolve with relaxation term */
     for (j=0; j<Relax.Nauto; j++) {
        /* printf("       interaction %d , secul.",j); */
        Nsec = *(int*)(Relax.secular[j]);
        for (k=1; k<=Nsec; k++) {
	   /* printf(" %d",k); */
	   q = Relax.secular[j][k].q;
	   y = Relax.secular[j][k].y;
	   /* printf(" => (%d,%d) \n",q,y); */
	   Q = Relax.Q[j][q];
	   Y = Relax.Y[j][y];
	   cm_commutator(mx1,Y,mx);
	   /* m_print(mx1,"[Y,rho]"); */
	   cm_commutator(mx2,Q,mx1);
	   /* m_print(mx2,"[Q,[Y,rho]]"); */
	   cm_commutator(mx1,Q,Y);
	   /* m_print(mx1,"[Q,Y]"); */
	   cmv_multod(mx2,mx1,-(Relax.eqnorm)*(Relax.omega[j][y]));
	   /* m_print(mx2,"relax prispevek"); */
	   cmv_multod(puls->sigma,mx2, -dt);
	   /* tohle je blbe, zanedbavam omq vuci om v oprave na rovnovaznej stav!!! */
	}
	/* printf(" done\n"); */
     }

  }
  
  complx_matrix_free(mx);
  complx_matrix_free(mx1);
  complx_matrix_free(mx2);
  puls->t_usec += duration*1.0e6;
}

