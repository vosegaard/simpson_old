/*
    Wigner rotation routines
    Copyright (C) 1999 Mads Bak
    
    2009: Z.T. changes for new structures of matrices

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
    
    Creating of Wigner roatation matrices -- either the full matrix
    or the zero'th column.
*/

#include <math.h>
#include "complx.h"
#include <stdio.h>
#include "defs.h"
#include "matrix_new.h"
//#include "defs_blas_lapack.h"
#include "cm_new.h"
#ifdef MKL_INTEL
#include "mkl.h"
#elif defined(__APPLE__)
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#endif

void wigner2(complx *d2,double alpha,double beta,double gamma)
{
  double cosb,sinb,cos2b,sin2b,cplus,cminus,gamma2,alpha2;
  double cplus2,cminus2,sinof2b,cplussinb,cminussinb,SQRT3BY8sinof2b,SQRT3BY8sin2b;
  complx em2am2g,em2amg,em2a,em2apg,em2ap2g,emam2g,emamg,ema;
  complx emapg,emap2g,em2g,emg,epg,ep2g,epam2g,epamg,epa;    
  complx epapg,epap2g,ep2am2g,ep2amg,ep2a,ep2apg,ep2ap2g;
  
  alpha *= DEG2RAD;
  beta *= DEG2RAD;
  gamma *= DEG2RAD;
  cosb=cos(beta);
  sinb=sin(beta);
  cos2b=cosb*cosb;
  sin2b=sinb*sinb;
  cplus=(1.0+cosb)*0.5;
  cminus=(1.0-cosb)*0.5;
  cplus2=cplus*cplus;
  cminus2=cminus*cminus;
  sinof2b=sin(2.0*beta);

  alpha2=-2.0*alpha;
  gamma2=2.0*gamma;
  
  cplussinb=cplus*sinb;
  cminussinb=cminus*sinb;
  SQRT3BY8sinof2b=SQRT3BY8*sinof2b;
  SQRT3BY8sin2b=SQRT3BY8*sin2b;
  
  em2am2g=Cexpi(alpha2-gamma2);
  em2amg =Cexpi(alpha2-gamma);
  em2a   =Cexpi(alpha2);
  em2apg =Cexpi(alpha2+gamma);
  em2ap2g=Cexpi(alpha2+gamma2);

  emam2g =Cexpi(-alpha-gamma2);
  emamg  =Cexpi(-alpha-gamma);
  ema    =Cexpi(-alpha);
  emapg  =Cexpi(-alpha+gamma);
  emap2g =Cexpi(-alpha+gamma2); 

  em2g   =Cexpi(-gamma2);
  emg    =Cexpi(-gamma);
  epg    =Conj(emg);
  ep2g   =Conj(em2g);

  epam2g =Conj(emap2g);
  epamg  =Conj(emapg);
  epa    =Conj(ema);
  epapg  =Conj(emamg);
  epap2g =Conj(emam2g);

  ep2am2g=Conj(em2ap2g);
  ep2amg =Conj(em2apg);
  ep2a   =Conj(em2a);
  ep2apg =Conj(em2amg);
  ep2ap2g=Conj(em2am2g);

  /* first column D_i-2 */
  d2[ 0] = RCmul(cplus2       ,ep2ap2g);
  d2[ 1] = RCmul(-cplussinb   ,epap2g);
  d2[ 2] = RCmul(SQRT3BY8sin2b   ,ep2g  );
  d2[ 3] = RCmul(-cminussinb  ,emap2g);
  d2[ 4] = RCmul(cminus2     ,em2ap2g );
  /* second column D_i-1 */
  d2[ 5] = RCmul(cplussinb   ,ep2apg );
  d2[ 6] = RCmul((cos2b-cminus),epapg );
  d2[ 7] = RCmul(-SQRT3BY8sinof2b,epg   );
  d2[ 8] = RCmul((cplus-cos2b) ,emapg );
  d2[ 9] = RCmul(-cminussinb,em2apg  );
  /* third column D_i0 */
  d2[10] = RCmul(SQRT3BY8sin2b  ,ep2a   );
  d2[11] = RCmul(SQRT3BY8sinof2b ,epa   );
  d2[12] = Complx(0.5*(3.0*cos2b-1.0),0.0);
  d2[13] = RCmul(-SQRT3BY8sinof2b,ema   );
  d2[14] = RCmul(SQRT3BY8sin2b ,em2a    );
  /* fourth column D_i+1 */
  d2[15] = RCmul(cminussinb  ,ep2amg );
  d2[16] = RCmul((cplus-cos2b) ,epamg );
  d2[17] = RCmul(SQRT3BY8sinof2b ,emg   );
  d2[18] = RCmul((cos2b-cminus),emamg );
  d2[19] = RCmul(-cplussinb ,em2amg  );
  /* fifth column D_i+2 */
  d2[20] = RCmul(cminus2      ,ep2am2g);
  d2[21] = RCmul(cminussinb   ,epam2g);
  d2[22] = RCmul(SQRT3BY8sin2b   ,em2g  );
  d2[23] = RCmul(cplussinb    ,emam2g);
  d2[24] = RCmul(cplus2      ,em2am2g );

  return;
}


void wigner20(complx* d20,double alpha,double beta)
{
  double cosb,sinb,sin2b,sinof2b;
  complx em2a,ema,epa,ep2a;

  alpha *= DEG2RAD;
  beta *= DEG2RAD;

  cosb   =cos(beta);
  sinb   =sin(beta);
  sin2b  =sinb*sinb;
  sinof2b=sin(2.0*beta);

  em2a   =Cexpi(-2.0*alpha);
  ema    =Cexpi(-alpha);
  epa    =Conj(ema);
  ep2a   =Conj(em2a);

  d20[4] =  RCmul(SQRT3BY8*sin2b,em2a);
  d20[3] =  RCmul(-SQRT3BY8*sinof2b,ema);
  d20[2] =  Complx(0.5*(3.0*cosb*cosb-1.0),0.0);
  d20[1] =  RCmul(SQRT3BY8*sinof2b,epa);
  d20[0] =  RCmul(SQRT3BY8*sin2b,ep2a);
}

void wig2rot(complx* res, complx* vec, complx *d2)
{
   const int N=5;
   
   zgemv_("T",&N,&N,&CPLX1,d2,&N,vec,&INTONE,&CPLX0,res,&INTONE);
}

double wig20rot(complx *vec, complx *d20)
{
   complx res;
   const int len=5;
   cblas_zdotu_sub(len,vec,INTONE,d20,INTONE,&res);
   if ( fabs(res.im) > 1e-8 ) {
      fprintf(stderr,"wig20rot error: result is not pure real\n");
      exit(1);
   }
   return res.re; 
}
