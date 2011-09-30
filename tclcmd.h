#ifndef _TCLCMD_H
#define TCLCMD_H

#include "config.h"

#ifdef __cplusplus
extern "C" {
#endif


void tclcmd_csa(Tcl_Interp* interp);
void tclcmd_csastat(Tcl_Interp* interp);
void tclcmd_fit(Tcl_Interp* interp);
void tclcmd_ftools(Tcl_Interp* interp);
void tclcmd_pulse(Tcl_Interp* interp);
void tclcmd_restrict(Tcl_Interp* interp);
void tclcmd_simpson(Tcl_Interp* interp);
void tclcmd_spinsys(Tcl_Interp* interp);
void tclcmd_tclutil(Tcl_Interp* interp);
void tclcmd_minuit(Tcl_Interp* interp);
void tclcmd_rfshape(Tcl_Interp* interp);
void tclcmd_OCroutines(Tcl_Interp* interp);
void tclcmd_inhom(Tcl_Interp* interp);
void tclcmd_zte(Tcl_Interp* interp);
void tclcmd_new_direct(Tcl_Interp* interp);

#ifdef ENABLE_MINUIT
#define DECLARE_TCL_COMMANDS \
tclcmd_csa(interp);\
tclcmd_csastat(interp);\
tclcmd_fit(interp);\
tclcmd_ftools(interp);\
tclcmd_pulse(interp);\
tclcmd_restrict(interp);\
tclcmd_simpson(interp);\
tclcmd_spinsys(interp);\
tclcmd_minuit(interp);\
tclcmd_tclutil(interp);\
tclcmd_rfshape(interp);\
tclcmd_OCroutines(interp);\
tclcmd_inhom(interp);\
tclcmd_zte(interp);\
tclcmd_new_direct(interp);
#else
#define DECLARE_TCL_COMMANDS \
tclcmd_csa(interp);\
tclcmd_csastat(interp);\
tclcmd_fit(interp);\
tclcmd_ftools(interp);\
tclcmd_pulse(interp);\
tclcmd_restrict(interp);\
tclcmd_simpson(interp);\
tclcmd_spinsys(interp);\
tclcmd_tclutil(interp);\
tclcmd_rfshape(interp);\
tclcmd_OCroutines(interp);\
tclcmd_inhom(interp);\
tclcmd_zte(interp);\
tclcmd_new_direct(interp);
#endif /* Minuit */

#ifdef __cplusplus
}
#endif

#endif
