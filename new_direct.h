#ifndef _new_direct_H
#define _new_direct_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tcl.h>
#include "matrix_new.h"
#include "cm_new.h"
#include "sim.h" 
#include "pulse.h"
#include "tclutil.h"
#include "defs.h"
#include "rfshapes.h"
#include "iodata.h"
#include "OCroutines.h"
#include "B0inhom.h"
#include "relax.h"

#define NEW_DIRECT_MXPOS_INI 900

void new_direct_initialize(Tcl_Interp *interp,Sim *s);
void new_direct_destroy(Sim *s);
void new_direct_calc(Sim *s,double gamma_add, mv_complx *wfidsum, mv_complx *wfid, int *nfid);
void new_gcompute_calc(Sim *s,double gamma_add, mv_complx *wfidsum, mv_complx *wfid, int *nfid);
void new_gcompute2_calc(Sim *s,double gamma_add, mv_complx *wfidsum, mv_complx *wfid, int *nfid);

#endif
