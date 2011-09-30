/*
    Tcl/C code performing the simulation
    Copyright (C) 1999 Mads Bak. Modifications by Thomas Vosegaard, 2001.

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
    
    Makes available the 'internalsimpson' command used by the
    'fsimpson' Tcl function located in 'simpson.tcl'
    Uses the commands given in sim.c that performs the simulation.
*/

#include <stdio.h>
#include <string.h>
#include <tcl.h>
#include <time.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include "config.h"
#include "cm_new.h"
#include "tclutil.h"
#include "iodata.h"
#include "sim.h"
#include "defs.h"
#include "cryst.h"
#include "rfshapes.h"
#include "OCroutines.h"
#include "B0inhom.h"

int getbits(int* e,char *bits)
{
  int i;
  char* p;

  *e=0;
  for (p=bits,i=0;*p != 0 && i < sizeof(int)*8;p++,i++) {
    if (*p == '1')
      *e |= 1 << i; 
    else if (*p != '0')
      return 0;
  }
  return 1;
}

void putbits(char *bits,int e)
{
  int i;
  char* p;

  for (p=bits,i=0; i < sizeof(int)*8;p++,i++) {
    if (e & (1 << i))
      *p = '1';
    else
      *p = '0';
  }
  *p=0;
}


#ifndef STDIN_FILENO
#define STDIN_FILENO  0
#define STDOUT_FILENO 1
#endif

int tclDummy(ClientData data, Tcl_Interp* interp, int argc, char *argv[])
{
  return TCL_OK;
}

/***
 * ZT: modified to include averaging over rf inhomogeneities, includes also changes made by AB
 ***/
#ifndef MPI //include header files for multiprocess setup
#include <sys/wait.h>
#include <sys/shm.h>
#include <sys/ipc.h>
#include <sys/stat.h>
#include <assert.h>
#include <errno.h>
#include <limits.h>
#endif

void initialize_slaves(int *np, int *rank, int *segment_id, void **shared_memory, Sim **ss, pid_t **slave_pids){
#ifdef _WIN32
  *np = 1;
  *rank = 0;
  return;
#elif defined(MPI)
  MPI_Comm_rank(MPI_COMM_WORLD, rank);
  MPI_Comm_size(MPI_COMM_WORLD, np);
  if (verbose & VERBOSE_SIMINFO) printf("Welcome to MPI Simpson, I am slave %i of %i\n", *rank, *np);
  return;
#else
  int Ntot = *np; /* Total number of jobs to be executed. Make sure that np<=Ntot below */
  /* RA: fork vars */
  pid_t pid, slave_pid = -1;
  int slaves, slave_rank = 1;
  /* RA: shared memory segment vars */
  int shared_segment_size;
  int i, *p;
  char *env_np; /* to check if SIMPSON_NUM_CORES is set */
  Sim *sim = *ss;

  int *pid_array;
  /* Look for number of processors to use. If not set by env. var., use all: */
  *np = sysconf (_SC_NPROCESSORS_CONF);
  if ((env_np = getenv("SIMPSON_NUM_CORES")) != NULL){
    char *endptr;
    errno = 0;
    *np = strtol(env_np, &endptr, 10);
    if ((errno == ERANGE && ((long)np == LONG_MAX || (long)np == LONG_MIN)) || (errno != 0 && np == 0)) {
      perror("strtol");
      *np = sysconf (_SC_NPROCESSORS_CONF);
      errno = 0;
    }
    if (endptr == env_np) {
      fprintf(stderr, "No digits were found, using all available cores\n");
      *np = sysconf (_SC_NPROCESSORS_CONF);
    }
  }
  /* If there are more jobs than available processes, just resort to Ntot processes: */
  if (*np > Ntot){
    *np = Ntot;
  }
  slaves = (*np)-1;
  pid_array = (int *) malloc(slaves * 2 * sizeof(pid_t));

  if (verbose & VERBOSE_SIMINFO) printf("Welcome to Simpson, number of processes: %i\n\n", *np);

  shared_segment_size = sim->ntot*sizeof(complx)*slaves;

  /* Allocate a shared memory segment */
  *segment_id = shmget(IPC_PRIVATE, shared_segment_size, IPC_CREAT | S_IRUSR | S_IWUSR);
  /* Attach the shared memory segment */
  *shared_memory = (void *)shmat(*segment_id, 0, 0);
  /* Determine the segment's size */
  //shmctl(segment_id, IPC_STAT, &shmbuffer);
  //segment_size = shmbuffer.shm_segsz;

  p = (int *)*shared_memory;
  /* Fork the processes and write rank of each slave process to shared memory */
  for (i=0; i<slaves; i++){
    slave_pid = fork();
    if (slave_pid == 0){
      /* slaves just break, only master spawns new processes */
      break;
    }
    p[i*2] = slave_pid;
    pid_array[i] = slave_pid;
    p[i*2+1] = slave_rank++;
  }
  /* Master, set rank to 0 (slaves will overwrite below) */
  *rank = 0;
  /* Save the pid array for the master to wait for those processes later on */
  *slave_pids = pid_array;

  /* Slave processes, read rank from shared memory */
  if (slave_pid == 0){ 
    pid = getpid();
    while (*rank == 0){
      for (i=0; i<slaves; i++){
	if (pid == p[i*2]){
	  *rank = p[i*2+1];
	  break;
	}
      }
    }
    //    printf("slave, pid is %i rank is %i\n", pid, *rank);
  }
  return;
#endif
}

void wait_for_slaves(int np, pid_t *slave_pids, void *shared_memory, Sim *ss, mv_complx *wfid, mv_complx **fidsum, int segment_id){
#ifdef _WIN32
  return;
#elif defined(MPI)
  int i;
  for (i=1; i<np; i++) {
    MPI_Recv(wfid->data,2*(wfid->row),MPI_DOUBLE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    /* update global vector with results from each local vector */
    cmv_addto(*fidsum,wfid);
  }
  return;
#else
  int i;
  int slave_status;
  for (i=0; i<np-1; i++){
    waitpid(slave_pids[i], &slave_status, 0);
    assert(WEXITSTATUS(slave_status) == 0);
    void *slave_segment = shared_memory + (i*ss->ntot*sizeof(complx));
    memcpy(wfid->data, slave_segment, ss->ntot*sizeof(complx));
    cmv_addto(*fidsum,wfid);
  }
  free(slave_pids);
  /* Detach shared memory segment */
  shmdt(shared_memory);
  /* deallocate the shared memory segment */
  shmctl(segment_id, IPC_RMID, 0);
  return;
#endif
}

void update_slaves(int np, mv_complx *fidsum){
#ifdef MPI
  /* send updated global vector to all workers */
  int i;
  for (i=1; i<np; i++) {
    MPI_Send(fidsum->data,2*(fidsum->row),MPI_DOUBLE,i,0,MPI_COMM_WORLD);
  }
  return;
#else
  /* In multiprocess setup and on WIN32, slaves are not updated...\n */
  return;
#endif 
}

void slave_finalize(int rank, void *shared_memory, Sim *ss, mv_complx *fidsum){
#ifdef _WIN32
  return;
#elif defined(MPI)
  MPI_Send(fidsum->data, 2*(fidsum->row), MPI_DOUBLE,0,0,MPI_COMM_WORLD);
  MPI_Recv(fidsum->data, 2*(fidsum->row), MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  //printf("process %i received resiult vector, data is %f\n", rank, fidsum->data[ss->ntot-1].im);
  //MPI_Finalize();
  return;
#else
  /* copy fidsum to shared memory */
  void *my_segment = shared_memory + ((rank-1)*ss->ntot*sizeof(complx));
  memcpy(my_segment, fidsum->data, ss->ntot*sizeof(complx));
  /* detach the shared memory segment */
  shmdt(shared_memory);
  return;
#endif
}

FD* simpson(Tcl_Interp* interp){
  mv_double * read_crystfile(Sim* s);
  void read_rfproffile(const char* name,int nchan, mv_double** rfsc,mv_double** rfw,double* sumweight);
  FD* fd;
  int ncr,nrf,i;
  double totalintensity, rfsumweight, totalweight;
  mv_double *crdata, *rfsc, *rfw;
  double *rfscrow, *rfscalefactors, rfweight, weight;
  complx *fid;
  mv_complx *fidsum, *wfid;
  Sim *ss;
  Omega omega;
  double *zvals=NULL, *zoffsetvals=NULL;
  int znumber;

  /* RA+ZT: vars for distributing jobs */
  int np, rank = 0;
  int Ntot,Np,rem, Np_start,ncr_start,nrf_start,nz_start,icr,irf,iz;

  /* RA: vars for multiprocess setup */
  void *shared_memory = NULL;
  int segment_id;
  pid_t *slave_pids = NULL;

  ss=(Sim*)malloc(sizeof(Sim));
  if (!ss) {
    fprintf(stderr,"error: unable to allocate Sim structure\n");
    exit(1);
  }
  sim_initialize(interp,ss);
  fid = complx_vector(ss->ntot); /* old style */
  wfid = (mv_complx*)malloc(sizeof(mv_complx)); /* wrapper to new style */
  wfid->row = ss->ntot; wfid->col = 1; wfid->data = &fid[1];
  fidsum = complx_vector_alloc(ss->ntot); /* new style */

  cmv_zero(wfid);
  cmv_zero(fidsum);

  /* Read crystal file */
  crdata = read_crystfile(ss);
  /* ZT: read also rf profile(s) */
  read_rfproffile(ss->rfproffile,ss->ss->nchan, &rfsc, &rfw, &rfsumweight);
  /* initializations for averaging over z coordinate */
  prepare_zaveraging(interp,&zvals,&zoffsetvals);

  ncr = crdata->row;
  nrf = rfsc->row;
  znumber = LEN(zvals);
  Ntot = ncr*nrf*znumber; /* Total number of jobs */
  totalintensity = 0.0;
  np = Ntot;
  initialize_slaves(&np, &rank, &segment_id, &shared_memory, &ss, &slave_pids);
  /* RA+ZT: divide the total work load into smaller task units */
  Np = Ntot/np; /* Number of jobs per processing element */
  rem = Ntot - Np*np;
  Np_start = rank*Np; /* */
  int local_rem = Ntot-rem+rank;

  /* another mess: converting 0-based vector to 1-based */
  rfscalefactors = double_vector(rfsc->col); /* old style */

  for (i=1;i<=ncr;i++) {
    totalintensity += double_get_elem(crdata,i,4);
  }
  i=0;
  while (1){
    ncr_start = Np_start/nrf/znumber;
    nrf_start = (Np_start-ncr_start*nrf*znumber)/znumber;
    nz_start = Np_start - nrf_start*znumber - ncr_start*nrf*znumber;
    icr = ncr_start+1;
    irf = nrf_start+1;
    iz = nz_start+1;
    omega.alpha = double_get_elem(crdata,icr,1);
    omega.beta = double_get_elem(crdata,icr,2);
    omega.gamma = double_get_elem(crdata,icr,3);
    weight = double_get_elem(crdata,icr,4);
    rfscrow = dm_row(rfsc,irf);
    /* another mess: converting 0-based vector to 1-based */
    memcpy(&rfscalefactors[1],rfscrow,(rfsc->col)*sizeof(double));
    free((char*)rfscrow);
    rfweight = rfw->data[irf-1];
    ss->P->zcoor = zvals[iz];
    set_inhom_offsets(ss,zoffsetvals[iz]);
    sim_calcfid(ss,omega,rfscalefactors,fid);
    totalweight = weight*rfweight/rfsumweight;
    cmv_multod(fidsum, wfid, totalweight);

    if (verbose & VERBOSE_PROGRESS) {
      printf("process %i: [cryst %d/%d] [rfsc %d/%d] [z-coor %d/%d]\n", rank,icr,ncr,irf,nrf,iz,znumber);
      fflush(stdout);
    }
    if (++i < Np){
      Np_start++;
    }
    else if ((rem > 0) && (local_rem<Ntot)){
      //printf("process %i, remainder to do is %i, Ntot is %i\n", rank, local_rem, Ntot);
      Np_start = local_rem;
      local_rem = Ntot;
    }
    else{
      break;
    }
  }
  /* Master: wait for slaves to exit and copy their
     results from the shared memory segment, or for
     MPI slaves to send their results */
  if (rank == 0){
    wait_for_slaves(np, slave_pids, shared_memory, ss, wfid, &fidsum, segment_id);
  }

  double_matrix_free(rfsc);
  double_vector_free(rfw);
  double_matrix_free(crdata);
  free_double_vector(zvals); /* old style */
  free_double_vector(zoffsetvals); /* old style */
  free_double_vector(rfscalefactors); /* old style */

  /* Slaves: copy the results to the shared memory, or
   * in MPI, send results, then cleanup and exit */
  if (rank != 0){
    slave_finalize(rank, shared_memory, ss, fidsum);
    /* MPI slaves do not exit because simpson might be called again,
       (MPI does not allow any MPI routine to called after they finalize)
    */
#ifndef MPI
    /* multiprocess slaves clean up before exit... */
    sim_destroy(ss);
    free_complx_vector(fid);
    free((char*)wfid);
    complx_vector_free(fidsum);
    free(ss);
    exit(0);
#endif
  }

  /* Master: global vector is now complete */
  if (rank == 0){
    totalintensity *= (double)(znumber);
    cmv_muld(fidsum,1.0/totalintensity);
    update_slaves(np, fidsum);
  }

  /* this is really a mess, now I need to convert it to old style again... */
  memcpy(&fid[1],fidsum->data,(fidsum->row)*sizeof(complx));
  fd=FD_data2fd(NULL,(double2*)fid,ss->np,ss->ni,ss->sw,ss->sw1);

  if (ss->sw1 != 0) fd->sw1=ss->sw1;
  if (ss->ni != 0) fd->ni=ss->ni;

  sim_destroy(ss);
  free_complx_vector(fid);
  free((char*)wfid);
  complx_vector_free(fidsum);
  /* free_wsp(); free_mv_static(); - this makes problems when fsimpson called repeatedly */
  free(ss);
  fd->rank = rank;
  return fd;
}

int tclCrystallites(ClientData data,Tcl_Interp* interp,
      int argc, char *argv[])
{
  int n=0,i;
  CRYSTALLITE* c;
  char name[64];
    
  if (argc != 2) {
    interp->result = "Usage: crystallites <crystal_file>";
    return TCL_ERROR;
  }
  
  sprintf(name, "%s_cryst", argv[1]);
  
  while (strlen(cryst_names[n])) {
    if (!strcmp(cryst_names[n],name)) {
      c = cryst_pointers[n];
      for (i=0;i<cryst_numbers[n];i++) {
        TclAppendResult(interp,"%g %g %g",c[i].alpha,c[i].beta,c[i].weight);
      }
      return TCL_OK;
    }
    n++;
  }

  return TclError(interp,"crystallites: unable to find crystal file '%s'", argv[1]);
}

int tclInternalSimpson(ClientData data,Tcl_Interp* interp,
      int argc, char *argv[]){
  char buf[256],buf1[16];
  int fnew(FD* f);
  int fidN;
  extern FD **fd;
  FD* f;

  if (argc != 1) {
    interp->result = "Usage: <desc> internalsimpson";
    return TCL_ERROR;
  }

  TclGetString(interp,buf,"par","verbose",0,"0");
  if (!getbits(&verbose,buf))
    return TclError(interp,"error: 'verbose' parameter contains illegal characters '%s' (must be a string of '0' and '1')",buf);

  TclGetString(interp,buf,"par","various",0,"0");
  if (!getbits(&various,buf))
    return TclError(interp,"error: 'various' parameter contains illegal characters '%s' (must be a string of '0' and '1')",buf);

  f=simpson(interp);
  
  fidN = fnew(f);
  /*  fprintf(stderr,"fidN = %d\n", fidN); */
  sprintf(buf1,"%d", fidN);
  sprintf(buf,"%ld", (long)fd[fidN]);
  Tcl_Eval(interp,"namespace eval FD {variable f}");
  Tcl_SetVar2(interp,"FD::f",buf1,buf,TCL_GLOBAL_ONLY|TCL_LEAVE_ERR_MSG);
  Tcl_SetVar2(interp,"FD_Internal",buf1,buf,TCL_GLOBAL_ONLY);
  sprintf(interp->result,"%d",fidN);

  return TCL_OK;
}


int tcltestsavestate(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  char *state;
  /* this is needed for duplicating RFshapes to Tcl  --- BEGIN --- */
  int a, isany=0;
  char buf[64];
   
  for (a=0; a<MAXRFSHAPES; a++) {
     if (RFshapes[a]) {
        printf("RFshape slot %d is active\n",a);
        Tcl_ResetResult(interp);
	sprintf(buf,"shape2list %d",a);
	if (Tcl_Eval(interp,buf) != TCL_OK) {
           fprintf(stderr,"error: unable to execute shape2list: '%s'\n",interp->result);
           exit(-1);
        }
	sprintf(buf,"%d",a);
	Tcl_SetVar2(interp,"__RFSHAPES",buf,Tcl_GetStringResult(interp),TCL_GLOBAL_ONLY);
	isany++;
     }
  }
  if (isany) {
     sprintf(buf,"%d",isany);
     Tcl_SetVar2(interp,"__RFSHAPES","info",buf,TCL_GLOBAL_ONLY);
  }
  /* this is needed for duplicating RFshapes to Tcl  --- END --- */
  
  /* this is to test move_OCpar_to_tcl */
  move_OCpar_to_tcl(interp);
  
  Tcl_ResetResult(interp);
  if (Tcl_Eval(interp,"savestate") != TCL_OK) {
     fprintf(stderr,"error: unable to get interpreter state '%s'\n",interp->result);
     exit(-1);
  }
  state = strdup(Tcl_GetStringResult(interp));
  
  printf("%s\n",state);

  return TCL_OK;
}

int tcltestvariable(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int c;
  char *src;
  
  if ((src=Tcl_GetVar2(interp,"__RFSHAPES","info",TCL_GLOBAL_ONLY)) == NULL) {
     c = 0;
  } else {
     if (Tcl_GetInt(interp,src,&c) != TCL_OK) TclError(interp,"cannot get int from __RFSHAPES(info)");
  }
  
  
  if (c) {
     printf("Variable __RFSHAPES exists! Copying data to RFshapes[]\n");
     int n, Nelem, i, n2;
     Tcl_Obj *listPtr, **listObjElem, **listObjElem2;
     char buf[8];
     double ampl, phase;
     
     RFshapes_reset();
     for (n=0; n<MAXRFSHAPES; n++) {
     sprintf(buf,"%d",n);
        if ( (listPtr=Tcl_GetVar2Ex(interp,"__RFSHAPES",buf,TCL_GLOBAL_ONLY)) == NULL) {
	   printf("   element %d is not present\n",n);
	   continue;
	}
	if (Tcl_ListObjGetElements(interp, listPtr, &Nelem, &listObjElem) != TCL_OK) 
	  return TclError(interp,"cannot decompose list in element %d\n",n);
	printf("   element %d decomposed into RFshape list\n",n);
	RFshapes[n] = RFshapes_alloc(Nelem);
	for (i=0; i<Nelem; i++) {
	   if (Tcl_ListObjGetElements(interp, listObjElem[i], &n2, &listObjElem2) != TCL_OK)
	     return TclError(interp,"cannot decompose ampl and phase in (%d, %d)",n,i);
	   if (Tcl_GetDoubleFromObj(interp,listObjElem2[0],&ampl) != TCL_OK)
	     return TclError(interp,"error in conversion to double");
	   if (Tcl_GetDoubleFromObj(interp,listObjElem2[1],&phase) != TCL_OK)
	     return TclError(interp,"error in conversion to double");
	RFshapes[n][i+1].ampl = ampl;
	RFshapes[n][i+1].phase = phase;
	printf("      (%f, %f)\n",ampl,phase);
	}
     }
     
  } else {
     printf("Variable __RFSHAPES doesn't exist!!! Removing all data in RFshapes[]\n");
     RFshapes_reset();
  }
   
 return TCL_OK;  
}

/* RA: After main.tcl has run the main section fron the input file,
       it will call this function which will finalize MPI
 */
int tclInternalFinalize(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
#ifdef MPI
  MPI_Finalize();
#endif
  return 0;
}

void tclcmd_simpson(Tcl_Interp* interp)
{
  int i;
  Tcl_CreateCommand(interp,"internalfinalize",tclInternalFinalize,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"internalsimpson",tclInternalSimpson,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
#ifndef DISABLE_NETWORK
  Tcl_CreateCommand(interp,"crystallites",tclCrystallites,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
#endif

/* ZT: testing savestate in noncluster mode:  */
Tcl_CreateCommand(interp,"testsavestate",tcltestsavestate,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateCommand(interp,"testvariable",tcltestvariable,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
} 
