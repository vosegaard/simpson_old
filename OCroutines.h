#ifndef _OCroutines_H
#define _OCroutines_H

#ifndef __COMPLX_H
#include "complx.h"
#endif 

#define MAXOCPROPS 5119 


/* ZT: structure to hold and pass OC parameters */
typedef struct _OCoptPars {

  int isinit,gradmode,gradmodeprop,propstatus,verb;

  double eps, tol, mnbkstep, cut, stepmin;
  int ndim, nIterations, nreport, ncut, max_brack_eval, max_brent_eval;
  char VarSaveProc[64];
  
  mv_complx * prop[MAXOCPROPS+1]; 
  mv_complx * dens[MAXOCPROPS+1]; 

  int mx_pos; 
  char *mx_code[MAXOCPROPS+1];
  
  int *var_shapes, *grad_shapes;
  double *var_shapes_min, *var_shapes_max, *var_shapes_rmsmax;
  
} OCoptPars;
  
/* global variable holding all OC parameters */
extern OCoptPars OCpar;


void store_OCprop(void);
void store_OCdens(void);
void set_OCmx_code(char *c);
void incr_OCmx_pos(void);
void _filterOC(int num);
void _pulse_shapedOC(char *code, int Nch, int Nelem, int *mask, double steptime);
void _pulse_shapedOCprops(char *code, int Nch, int Nelem, int *mask, double steptime);
void _pulse_and_zgrad_shapedOC(char *code, int Nch, int Nelem, int *mask, int zgrs, double steptime);
void _pulse_and_zgrad_shapedOCprops(char *code, int Nch, int Nelem, int *mask, int zgrs, double steptime);

void move_OCpar_to_tcl(Tcl_Interp*);
void create_OCpar_from_tcl(Tcl_Interp* interp);

void test_print_codes(void);

#endif
