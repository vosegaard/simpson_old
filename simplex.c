/*
    Simplex minimization, extract from fortran code by Tom Rowan

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
  

    Extract of the simplex function from the subplex package
    Converted from fortran to C with f2c

    AUTHOR of the original simplex/subplex fortran code

         Tom Rowan 
         Computer Science and Mathematics Division 
         Oak Ridge National Laboratory 
         P.O. Box 2008, Bldg. 6012 
         Oak Ridge, TN 37831-6367 
         Phone: (423) 574-3131 
         Fax : (423) 574-0680 
         Email: na.rowan@na-net.ornl.gov 
         URL: http://www.epm.ornl.gov/~rowan/ 

    REFERENCE 

         Thomas H. Rowan, Functional Stability Analysis of Numerical Algorithms,
         Ph.D. dissertation, Department of Computer Sciences,
         University of Texas at Austin, 1990. 
*/

#include <stdio.h>
#include <stdlib.h>

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

struct {
    double alpha, beta, gamma, delta, psi, omega;
    int nsmin, nsmax, irepl, ifxsw;
    double bonus, fstop;
    int nfstop, nfxe;
    double fxstat[4], ftest;
    int minf, initx, newx;
} usubc_;

#define usubc_1 usubc_

struct {
    double fbonus, sfstop, sfbest;
    int new__;
} isubc_;

#define isubc_1 isubc_
#define TRUE_ 1
#define FALSE_ 0

/* Table of constant values */

static int c_true = TRUE_;
static int c__1 = 1;
static int c__0 = 0;
static int c_false = FALSE_;

static double c_b3 = 0.;
static double c_b7 = 1.;

/* Subroutine */ int fstats_(fx, ifxwt, reset)
double *fx;
int *ifxwt;
int *reset;
{
    /* System generated locals */
    double d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static double fscale;
    static int nsv;
    static double f1sv;



/*                                         Coded by Tom Rowan */
/*                            Department of Computer Sciences */
/*                              University of Texas at Austin */

/* fstats modifies the common /usubc/ variables nfxe,fxstat. */

/* input */

/*   fx     - most recent evaluation of f at best x */

/*   ifxwt  - int weight for fx */

/*   reset  - int switch */
/*            = .true.  : initialize nfxe,fxstat */
/*            = .false. : update nfxe,fxstat */

/* common */



/* local variables */



/* subroutines and functions */

/*   fortran */

/* ----------------------------------------------------------- */

    if (*reset) {
	usubc_1.nfxe = *ifxwt;
	usubc_1.fxstat[0] = *fx;
	usubc_1.fxstat[1] = *fx;
	usubc_1.fxstat[2] = *fx;
	usubc_1.fxstat[3] = 0.;
    } else {
	nsv = usubc_1.nfxe;
	f1sv = usubc_1.fxstat[0];
	usubc_1.nfxe += *ifxwt;
	usubc_1.fxstat[0] += *ifxwt * (*fx - usubc_1.fxstat[0]) / 
		usubc_1.nfxe;
	usubc_1.fxstat[1] = max(usubc_1.fxstat[1],*fx);
	usubc_1.fxstat[2] = min(usubc_1.fxstat[2],*fx);
/* Computing MAX */
	d__1 = abs(usubc_1.fxstat[1]), d__2 = abs(usubc_1.fxstat[2]), d__1 = 
		max(d__1,d__2);
	fscale = max(d__1,1.);
/* Computing 2nd power */
	d__1 = usubc_1.fxstat[3] / fscale;
/* Computing 2nd power */
	d__2 = (usubc_1.fxstat[0] - f1sv) / fscale;
/* Computing 2nd power */
	d__3 = (*fx - usubc_1.fxstat[0]) / fscale;
	usubc_1.fxstat[3] = fscale * sqrt(((nsv - 1) * (d__1 * d__1) + nsv * (
		d__2 * d__2) + *ifxwt * (d__3 * d__3)) / (usubc_1.nfxe - 1));
    }
    return 0;
} /* fstats_ */


double dist_(n, x, y)
int *n;
double *x, *y;
{
    /* System generated locals */
    int i__1;
    double ret_val, d__1;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static int i__;
    static double scale, absxmy, sum;



/*                                         Coded by Tom Rowan */
/*                            Department of Computer Sciences */
/*                              University of Texas at Austin */

/* dist calculates the distance between the points x,y. */

/* input */

/*   n      - number of components */

/*   x      - point in n-space */

/*   y      - point in n-space */

/* local variables */


/* subroutines and functions */

/*   fortran */

/* ----------------------------------------------------------- */

    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    absxmy = (d__1 = x[1] - y[1], abs(d__1));
    if (absxmy <= 1.) {
	sum = absxmy * absxmy;
	scale = 1.;
    } else {
	sum = 1.;
	scale = absxmy;
    }
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	absxmy = (d__1 = x[i__] - y[i__], abs(d__1));
	if (absxmy <= scale) {
/* Computing 2nd power */
	    d__1 = absxmy / scale;
	    sum += d__1 * d__1;
	} else {
/* Computing 2nd power */
	    d__1 = scale / absxmy;
	    sum = sum * (d__1 * d__1) + 1.;
	    scale = absxmy;
	}
/* L10: */
    }
    ret_val = scale * sqrt(sum);
    return ret_val;
} /* dist_ */


/* Subroutine */ int subopt_(n)
int *n;
{


/*                                         Coded by Tom Rowan */
/*                            Department of Computer Sciences */
/*                              University of Texas at Austin */

/* subopt sets options for subplx. */

/* input */

/*   n      - problem dimension */

/* common */




/* subroutines and functions */

/*   fortran */

/* ----------------------------------------------------------- */

/* *********************************************************** */
/* simplex method strategy parameters */
/* *********************************************************** */

/* alpha  - reflection coefficient */
/*          alpha .gt. 0 */

    usubc_1.alpha = 1.;

/* beta   - contraction coefficient */
/*          0 .lt. beta .lt. 1 */

    usubc_1.beta = .5;

/* gamma  - expansion coefficient */
/*          gamma .gt. 1 */

    usubc_1.gamma = 2.;

/* delta  - shrinkage (massive contraction) coefficient */
/*          0 .lt. delta .lt. 1 */

    usubc_1.delta = .5;

/* *********************************************************** */
/* subplex method strategy parameters */
/* *********************************************************** */

/* psi    - simplex reduction coefficient */
/*          0 .lt. psi .lt. 1 */

    usubc_1.psi = .25;

/* omega  - step reduction coefficient */
/*          0 .lt. omega .lt. 1 */

    usubc_1.omega = .1;

/* nsmin and nsmax specify a range of subspace dimensions. */
/* In addition to satisfying  1 .le. nsmin .le. nsmax .le. n, */
/* nsmin and nsmax must be chosen so that n can be expressed */
/* as a sum of positive ints where each of these ints */
/* ns(i) satisfies   nsmin .le. ns(i) .ge. nsmax. */
/* Specifically, */
/*     nsmin*ceil(n/nsmax) .le. n   must be true. */

/* nsmin  - subspace dimension minimum */

    usubc_1.nsmin = min(2,*n);

/* nsmax  - subspace dimension maximum */

    usubc_1.nsmax = min(5,*n);

/* *********************************************************** */
/* subplex method special cases */
/* *********************************************************** */
/* nelder-mead simplex method with periodic restarts */
/*   nsmin = nsmax = n */
/* *********************************************************** */
/* nelder-mead simplex method */
/*   nsmin = nsmax = n, psi = small positive */
/* *********************************************************** */

/* irepl, ifxsw, and bonus deal with measurement replication. */
/* Objective functions subject to large amounts of noise can */
/* cause an optimization method to halt at a false optimum. */
/* An expensive solution to this problem is to evaluate f */
/* several times at each point and return the average (or max */
/* or min) of these trials as the function value.  subplx */
/* performs measurement replication only at the current best */
/* point. The longer a point is retained as best, the more */
/* accurate its function value becomes. */

/* The common variable nfxe contains the number of function */
/* evaluations at the current best point. fxstat contains the */
/* mean, max, min, and standard deviation of these trials. */

/* irepl  - measurement replication switch */
/* irepl  = 0, 1, or 2 */
/*        = 0 : no measurement replication */
/*        = 1 : subplx performs measurement replication */
/*        = 2 : user performs measurement replication */
/*              (This is useful when optimizing on the mean, */
/*              max, or min of trials is insufficient. Common */
/*              variable initx is true for first function */
/*              evaluation. newx is true for first trial at */
/*              this point. The user uses subroutine fstats */
/*              within his objective function to maintain */
/*              fxstat. By monitoring newx, the user can tell */
/*              whether to return the function evaluation */
/*              (newx = .true.) or to use the new function */
/*              evaluation to refine the function evaluation */
/*              of the current best point (newx = .false.). */
/*              The common variable ftest gives the function */
/*              value that a new point must beat to be */
/*              considered the new best point.) */

    usubc_1.irepl = 0;

/* ifxsw  - measurement replication optimization switch */
/* ifxsw  = 1, 2, or 3 */
/*        = 1 : retain mean of trials as best function value */
/*        = 2 : retain max */
/*        = 3 : retain min */

    usubc_1.ifxsw = 1;

/* Since the current best point will also be the most */
/* accurately evaluated point whenever irepl .gt. 0, a bonus */
/* should be added to the function value of the best point */
/* so that the best point is not replaced by a new point */
/* that only appears better because of noise. */
/* subplx uses bonus to determine how many multiples of */
/* fxstat(4) should be added as a bonus to the function */
/* evaluation. (The bonus is adjusted automatically by */
/* subplx when ifxsw or minf is changed.) */

/* bonus  - measurement replication bonus coefficient */
/*          bonus .ge. 0 (normally, bonus = 0 or 1) */
/*        = 0 : bonus not used */
/*        = 1 : bonus used */

    usubc_1.bonus = 1.;

/* nfstop = 0 : f(x) is not tested against fstop */
/*        = 1 : if f(x) has reached fstop, subplx returns */
/*              iflag = 2 */
/*        = 2 : (only valid when irepl .gt. 0) */
/*              if f(x) has reached fstop and */
/*              nfxe .gt. nfstop, subplx returns iflag = 2 */

    usubc_1.nfstop = 0;

/* fstop  - f target value */
/*          Its usage is determined by the value of nfstop. */

/* minf   - int switch */
/*        = .true.  : subplx performs minimization */
/*        = .false. : subplx performs maximization */

    usubc_1.minf = TRUE_;
    return 0;
} /* subopt_ */


/* Subroutine */ int newpt_(ns, coef, xbase, xold, new__, xnew, small)
int *ns;
double *coef, *xbase, *xold;
int *new__;
double *xnew;
int *small;
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i__;
    static int eqold;
    static double xoldi;
    static int eqbase;



/*                                         Coded by Tom Rowan */
/*                            Department of Computer Sciences */
/*                              University of Texas at Austin */

/* newpt performs reflections, expansions, contractions, and */
/* shrinkages (massive contractions) by computing: */

/* xbase + coef * (xbase - xold) */

/* The result is stored in xnew if new .eq. .true., */
/* in xold otherwise. */

/* use :  coef .gt. 0 to reflect */
/*        coef .lt. 0 to expand, contract, or shrink */

/* input */

/*   ns     - number of components (subspace dimension) */

/*   coef   - one of four simplex method coefficients */

/*   xbase  - double precision ns-vector representing base */
/*            point */

/*   xold   - double precision ns-vector representing old */
/*            point */

/*   new    - int switch */
/*            = .true.  : store result in xnew */
/*            = .false. : store result in xold, xnew is not */
/*                        referenced */

/* output */

/*   xold   - unchanged if new .eq. .true., contains new */
/*            point otherwise */

/*   xnew   - double precision ns-vector representing new */
/*            point if  new .eq. .true., not referenced */
/*            otherwise */

/*   small  - int flag */
/*            = .true.  : coincident points */
/*            = .false. : otherwise */

/* local variables */


/* subroutines and functions */

/*   fortran */

/* ----------------------------------------------------------- */

    /* Parameter adjustments */
    --xold;
    --xbase;
    --xnew;

    /* Function Body */
    eqbase = TRUE_;
    eqold = TRUE_;
    if (*new__) {
	i__1 = *ns;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    xnew[i__] = xbase[i__] + *coef * (xbase[i__] - xold[i__]);
	    eqbase = eqbase && xnew[i__] == xbase[i__];
	    eqold = eqold && xnew[i__] == xold[i__];
/* L10: */
	}
    } else {
	i__1 = *ns;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    xoldi = xold[i__];
	    xold[i__] = xbase[i__] + *coef * (xbase[i__] - xold[i__]);
	    eqbase = eqbase && xold[i__] == xbase[i__];
	    eqold = eqold && xold[i__] == xoldi;
/* L20: */
	}
    }
    *small = eqbase || eqold;
    return 0;
} /* newpt_ */


/* Subroutine */ int start_(n, x, step, ns, ips, s, small)
int *n;
double *x, *step;
int *ns, *ips;
double *s;
int *small;
{
    /* System generated locals */
    int s_dim1, s_offset, i__1;

    /* Local variables */
    static int i__, j;
    extern /* Subroutine */ int dcopy_();



/*                                         Coded by Tom Rowan */
/*                            Department of Computer Sciences */
/*                              University of Texas at Austin */

/* start creates the initial simplex for simplx minimization. */

/* input */

/*   n      - problem dimension */

/*   x      - current best point */

/*   step   - stepsizes for corresponding components of x */

/*   ns     - subspace dimension */

/*   ips    - permutation vector */


/* output */

/*   s      - first ns+1 columns contain initial simplex */

/*   small  - int flag */
/*            = .true.  : coincident points */
/*            = .false. : otherwise */

/* local variables */


/* subroutines and functions */

/*   blas */
/*   fortran */

/* ----------------------------------------------------------- */

    /* Parameter adjustments */
    --ips;
    --step;
    --x;
    s_dim1 = *ns;
    s_offset = s_dim1 + 1;
    s -= s_offset;

    /* Function Body */
    i__1 = *ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s[i__ + s_dim1] = x[ips[i__]];
/* L10: */
    }
    i__1 = *ns + 1;
    for (j = 2; j <= i__1; ++j) {
	dcopy_(ns, &s[s_dim1 + 1], &c__1, &s[j * s_dim1 + 1], &c__1);
	s[j - 1 + j * s_dim1] = s[j - 1 + s_dim1] + step[ips[j - 1]];
/* L20: */
    }

/* check for coincident points */

    i__1 = *ns + 1;
    for (j = 2; j <= i__1; ++j) {
	if (s[j - 1 + j * s_dim1] == s[j - 1 + s_dim1]) {
	    goto L40;
	}
/* L30: */
    }
    *small = FALSE_;
    return 0;

/* coincident points */

L40:
    *small = TRUE_;
    return 0;
} /* start_ */



/* Subroutine */ int order_(npts, fs, il, is, ih)
int *npts;
double *fs;
int *il, *is, *ih;
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i__, j, il0;



/*                                         Coded by Tom Rowan */
/*                            Department of Computer Sciences */
/*                              University of Texas at Austin */

/* order determines the indices of the vertices with the */
/* lowest, second highest, and highest function values. */

/* input */

/*   npts   - number of points in simplex */

/*   fs     - double precision vector of function values of */
/*            simplex */

/*   il     - index to vertex with lowest function value */

/* output */

/*   il     - new index to vertex with lowest function value */

/*   is     - new index to vertex with second highest */
/*            function value */

/*   ih     - new index to vertex with highest function value */

/* local variables */


/* subroutines and functions */

/*   fortran */

/* ----------------------------------------------------------- */

    /* Parameter adjustments */
    --fs;

    /* Function Body */
    il0 = *il;
    j = il0 % *npts + 1;
    if (fs[j] >= fs[*il]) {
	*ih = j;
	*is = il0;
    } else {
	*ih = il0;
	*is = j;
	*il = j;
    }
    i__1 = il0 + *npts - 2;
    for (i__ = il0 + 1; i__ <= i__1; ++i__) {
	j = i__ % *npts + 1;
	if (fs[j] >= fs[*ih]) {
	    *is = *ih;
	    *ih = j;
	} else if (fs[j] > fs[*is]) {
	    *is = j;
	} else if (fs[j] < fs[*il]) {
	    *il = j;
	}
/* L10: */
    }
    return 0;
} /* order_ */


/* Subroutine */ int calcc_(ns, s, ih, inew, updatc, c__)
int *ns;
double *s;
int *ih, *inew;
int *updatc;
double *c__;
{
    /* System generated locals */
    int s_dim1, s_offset, i__1;
    double d__1;

    /* Local variables */
    static int i__, j;
    extern /* Subroutine */ int dscal_(), dcopy_(), daxpy_();



/*                                         Coded by Tom Rowan */
/*                            Department of Computer Sciences */
/*                              University of Texas at Austin */

/* calcc calculates the centroid of the simplex without the */
/* vertex with highest function value. */

/* input */

/*   ns     - subspace dimension */

/*   s      - double precision work space of dimension .ge. */
/*            ns*(ns+3) used to store simplex */

/*   ih     - index to vertex with highest function value */

/*   inew   - index to new point */

/*   updatc - int switch */
/*            = .true.  : update centroid */
/*            = .false. : calculate centroid from scratch */

/*   c      - centroid of the simplex without vertex with */
/*            highest function value */

/* output */

/*   c      - new centroid */

/* local variables */


/* subroutines and functions */

/*   blas */

/* ----------------------------------------------------------- */

    /* Parameter adjustments */
    --c__;
    s_dim1 = *ns;
    s_offset = s_dim1 + 1;
    s -= s_offset;

    /* Function Body */
    if (*updatc) {
	if (*ih == *inew) {
	    return 0;
	}
	i__1 = *ns;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    c__[i__] += (s[i__ + *inew * s_dim1] - s[i__ + *ih * s_dim1]) / *
		    ns;
/* L10: */
	}
    } else {
	dcopy_(ns, &c_b3, &c__0, &c__[1], &c__1);
	i__1 = *ns + 1;
	for (j = 1; j <= i__1; ++j) {
	    if (j != *ih) {
		daxpy_(ns, &c_b7, &s[j * s_dim1 + 1], &c__1, &c__[1], &c__1);
	    }
/* L20: */
	}
	d__1 = 1. / *ns;
	dscal_(ns, &d__1, &c__[1], &c__1);
    }
    return 0;
} /* calcc_ */


/* Subroutine */ int evalf_(f, ns, ips, xs, n, x, sfx, nfe)
double (*f) ();
int *ns, *ips;
double *xs;
int *n;
double *x, *sfx;
int *nfe;
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i__;
    static double fx;
    static int newbst;
    extern /* Subroutine */ int fstats_();



/*                                         Coded by Tom Rowan */
/*                            Department of Computer Sciences */
/*                              University of Texas at Austin */

/* evalf evaluates the function f at a point defined by x */
/* with ns of its components replaced by those in xs. */

/* input */

/*   f      - user supplied function f(n,x) to be optimized */

/*   ns     - subspace dimension */

/*   ips    - permutation vector */

/*   xs     - double precision ns-vector to be mapped to x */

/*   n      - problem dimension */

/*   x      - double precision n-vector */

/*   nfe    - number of function evaluations */

/* output */

/*   sfx    - signed value of f evaluated at x */

/*   nfe    - incremented number of function evaluations */

/* common */





/* local variables */



/* subroutines and functions */


/* ----------------------------------------------------------- */

    /* Parameter adjustments */
    --ips;
    --xs;
    --x;

    /* Function Body */
    i__1 = *ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[ips[i__]] = xs[i__];
/* L10: */
    }
    usubc_1.newx = isubc_1.new__ || usubc_1.irepl != 2;
    fx = (*f)(n, &x[1]);
    if (usubc_1.irepl == 0) {
	if (usubc_1.minf) {
	    *sfx = fx;
	} else {
	    *sfx = -fx;
	}
    } else if (isubc_1.new__) {
	if (usubc_1.minf) {
	    *sfx = fx;
	    newbst = fx < usubc_1.ftest;
	} else {
	    *sfx = -fx;
	    newbst = fx > usubc_1.ftest;
	}
	if (usubc_1.initx || newbst) {
	    if (usubc_1.irepl == 1) {
		fstats_(&fx, &c__1, &c_true);
	    }
	    usubc_1.ftest = fx;
	    isubc_1.sfbest = *sfx;
	}
    } else {
	if (usubc_1.irepl == 1) {
	    fstats_(&fx, &c__1, &c_false);
	    fx = usubc_1.fxstat[usubc_1.ifxsw - 1];
	}
	usubc_1.ftest = fx + isubc_1.fbonus * usubc_1.fxstat[3];
	if (usubc_1.minf) {
	    *sfx = usubc_1.ftest;
	    isubc_1.sfbest = fx;
	} else {
	    *sfx = -usubc_1.ftest;
	    isubc_1.sfbest = -fx;
	}
    }
    ++(*nfe);
    return 0;
} /* evalf_ */



/* Subroutine */ int simplx_(f, n, step, ns, ips, maxnfe, cmode, x, fx, nfe, 
	s, fs, iflag)
double (*f) ();
int *n;
double *step;
int *ns, *ips, *maxnfe;
int *cmode;
double *x, *fx;
int *nfe;
double *s, *fs;
int *iflag;
{
    /* System generated locals */
    int s_dim1, s_offset, i__1;
    double d__1, d__2;

    /* Local variables */
    static int inew;
    extern double dist_();
    static int npts;
    extern /* Subroutine */ int calcc_();
    static int i__, j;
    extern /* Subroutine */ int evalf_();
    static int icent;
    static int small;
    extern /* Subroutine */ int order_();
    static int itemp;
    extern /* Subroutine */ int dcopy_(), newpt_(), start_();
    static double fc, fe;
    static int ih, il;
    static double fr;
    static int is;
    static int updatc;
    static double dum, tol;



/*                                         Coded by Tom Rowan */
/*                            Department of Computer Sciences */
/*                              University of Texas at Austin */

/* simplx uses the Nelder-Mead simplex method to minimize the */
/* function f on a subspace. */

/* input */

/*   f      - function to be minimized, declared external in */
/*            calling routine */

/*   n      - problem dimension */

/*   step   - stepsizes for corresponding components of x */

/*   ns     - subspace dimension */

/*   ips    - permutation vector */

/*   maxnfe - maximum number of function evaluations */

/*   cmode  - int switch */
/*            = .true.  : continuation of previous call */
/*            = .false. : first call */

/*   x      - starting guess for minimum */

/*   fx     - value of f at x */

/*   nfe    - number of function evaluations */

/*   s      - double precision work array of dimension .ge. */
/*            ns*(ns+3) used to store simplex */

/*   fs     - double precision work array of dimension .ge. */
/*            ns+1 used to store function values of simplex */
/*            vertices */

/* output */

/*   x      - computed minimum */

/*   fx     - value of f at x */

/*   nfe    - incremented number of function evaluations */

/*   iflag  - error flag */
/*            = -1 : maxnfe exceeded */
/*            =  0 : simplex reduced by factor of psi */
/*            =  1 : limit of machine precision */
/*            =  2 : reached fstop */

/* common */





/* local variables */



/* subroutines and functions */

/*   blas */
/*   fortran */

/* ----------------------------------------------------------- */

    /* Parameter adjustments */
    --x;
    --step;
    --fs;
    s_dim1 = *ns;
    s_offset = s_dim1 + 1;
    s -= s_offset;
    --ips;

    /* Function Body */
    if (*cmode) {
	goto L50;
    }
    npts = *ns + 1;
    icent = *ns + 2;
    itemp = *ns + 3;
    updatc = FALSE_;
    start_(n, &x[1], &step[1], ns, &ips[1], &s[s_offset], &small);
    if (small) {
	*iflag = 1;
	return 0;
    }
    if (usubc_1.irepl > 0) {
	isubc_1.new__ = FALSE_;
	evalf_(f, ns, &ips[1], &s[s_dim1 + 1], n, &x[1], &fs[1], nfe);
      if (fs[1] <= 0.0) return 0; /* added by Mads Bak */
    } else {
	fs[1] = *fx;
    }
    isubc_1.new__ = TRUE_;
    i__1 = npts;
    for (j = 2; j <= i__1; ++j) {
	evalf_(f, ns, &ips[1], &s[j * s_dim1 + 1], n, &x[1], &fs[j], nfe);
      if (fs[j] <= 0.0) return 0; /* added by Mads Bak */
/* L10: */
    }
    il = 1;
    order_(&npts, &fs[1], &il, &is, &ih);
    tol = usubc_1.psi * dist_(ns, &s[ih * s_dim1 + 1], &s[il * s_dim1 + 1]);

/*     main loop */

L20:
    calcc_(ns, &s[s_offset], &ih, &inew, &updatc, &s[icent * s_dim1 + 1]);
    updatc = TRUE_;
    inew = ih;

/*       reflect */

    newpt_(ns, &usubc_1.alpha, &s[icent * s_dim1 + 1], &s[ih * s_dim1 + 1], &
	    c_true, &s[itemp * s_dim1 + 1], &small);
    if (small) {
	goto L40;
    }
    evalf_(f, ns, &ips[1], &s[itemp * s_dim1 + 1], n, &x[1], &fr, nfe);
    if (fr <= 0.0) return 0; /* added by Mads Bak */
    if (fr < fs[il]) {

/*         expand */

	d__1 = -usubc_1.gamma;
	newpt_(ns, &d__1, &s[icent * s_dim1 + 1], &s[itemp * s_dim1 + 1], &
		c_true, &s[ih * s_dim1 + 1], &small);
	if (small) {
	    goto L40;
	}
	evalf_(f, ns, &ips[1], &s[ih * s_dim1 + 1], n, &x[1], &fe, nfe);
  if (fe <= 0.0) return 0; /* added by Mads Bak */
	if (fe < fr) {
	    fs[ih] = fe;
	} else {
	    dcopy_(ns, &s[itemp * s_dim1 + 1], &c__1, &s[ih * s_dim1 + 1], &
		    c__1);
	    fs[ih] = fr;
	}
    } else if (fr < fs[is]) {

/*         accept reflected point */

	dcopy_(ns, &s[itemp * s_dim1 + 1], &c__1, &s[ih * s_dim1 + 1], &c__1);
	fs[ih] = fr;
    } else {

/*         contract */

	if (fr > fs[ih]) {
	    d__1 = -usubc_1.beta;
	    newpt_(ns, &d__1, &s[icent * s_dim1 + 1], &s[ih * s_dim1 + 1], &
		    c_true, &s[itemp * s_dim1 + 1], &small);
	} else {
	    d__1 = -usubc_1.beta;
	    newpt_(ns, &d__1, &s[icent * s_dim1 + 1], &s[itemp * s_dim1 + 1], 
		    &c_false, &dum, &small);
	}
	if (small) {
	    goto L40;
	}
	evalf_(f, ns, &ips[1], &s[itemp * s_dim1 + 1], n, &x[1], &fc, nfe);
  if (fc <= 0.0) return 0; /* added by Mads Bak */
/* Computing MIN */
	d__1 = fr, d__2 = fs[ih];
	if (fc < min(d__1,d__2)) {
	    dcopy_(ns, &s[itemp * s_dim1 + 1], &c__1, &s[ih * s_dim1 + 1], &
		    c__1);
	    fs[ih] = fc;
	} else {

/*           shrink simplex */

	    i__1 = npts;
	    for (j = 1; j <= i__1; ++j) {
		if (j != il) {
		    d__1 = -usubc_1.delta;
		    newpt_(ns, &d__1, &s[il * s_dim1 + 1], &s[j * s_dim1 + 1],
			     &c_false, &dum, &small);
		    if (small) {
			goto L40;
		    }
		    evalf_(f, ns, &ips[1], &s[j * s_dim1 + 1], n, &x[1], &fs[
			    j], nfe);
        if (fs[j] <= 0.0) return 0; /* added by Mads Bak */
		}
/* L30: */
	    }
	}
	updatc = FALSE_;
    }
    order_(&npts, &fs[1], &il, &is, &ih);

/*       check termination */

L40:
    if (usubc_1.irepl == 0) {
	*fx = fs[il];
    } else {
	*fx = isubc_1.sfbest;
    }
L50:
  /*DEBUG CODE*/
	goto L20;
    if (usubc_1.nfstop > 0 && *fx <= isubc_1.sfstop && usubc_1.nfxe >= 
	    usubc_1.nfstop) {
	*iflag = 2;
    } else if (*nfe >= *maxnfe) {
	*iflag = -1;
    } else if (dist_(ns, &s[ih * s_dim1 + 1], &s[il * s_dim1 + 1]) <= tol || 
	    small) {
	*iflag = 0;
    } else {
	goto L20;
    }

/*     end main loop, return best point */

    i__1 = *ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[ips[i__]] = s[i__ + il * s_dim1];
/* L60: */
    }

    return 0;
} /* simplx_ */



/* The rest is coded by Mads Bak, Copyright 2000 */

int simplex_(double (*f) (), int *n, double *tol, int *maxnfe,
             int *mode,  double *scale,  double *x,  double *fx,
             int *nfe,   int *iflag, double *awork,double *bwork,int* iwork)
{

    /* Local variables */
    static int i__;
    static int cmode;
    static int ns;
    static double dum;
    static double sfx;

    --x;
    --scale;

	  subopt_(n);

	  usubc_1.ftest = 0.;
	  cmode = FALSE_;
	  isubc_1.new__ = TRUE_;
	  usubc_1.initx = TRUE_;
	  evalf_(f, &c__0, &scale[1], &dum, n, &x[1], &sfx, nfe);
    if (sfx <= 0.0) return 0;
	  usubc_1.initx = FALSE_;

    ns= *n;

    for (i__ = 1; i__ <= ns; i__++) {
      iwork[i__] = i__;
    }
    simplx_(f, n, &scale[1], n, &iwork[1], maxnfe, &cmode, &x[1], 
      &sfx, nfe, awork, bwork, iflag);

  return 0;
} /* subplx_ */


double (*simplex_driver)(double* v);

#define MAXN 256

double simplex_fun(int *n,double *v)
{
   int i;
   double x[MAXN];
   for (i=1;i<= *n;i++) {
     x[i]=v[i-1];
   }  
   return simplex_driver(x);   
}

int simplx(double (*func)(double* v),int n,double *x,double *scale, double *fx)
{
  int i;
  static double awork[MAXN*(MAXN+4)],bwork[MAXN+2];
  static int iwork[MAXN+1];  
  static int ns,nfe,mode,iflag,_n,maxnfe;
  static double _fx,_tol,_scale[MAXN],_x[MAXN];
  
  if (n >= MAXN) {
    fprintf(stderr,"error: maximum number of variables reached in subplex, increase MAXN\n");
    exit(1);
  }  
  _n=n;
  ns=n;
  simplex_driver=func;
  for (i=1;i<=n;i++) {
    _scale[i-1]=scale[i];
    _x[i-1]=x[i];      
  }
  _tol=0;
  mode=0;
  nfe=1;
  maxnfe=1e6;
  iflag=0;

  simplex_(simplex_fun,&_n,&_tol,&maxnfe,&mode,_scale,
            _x,&_fx,&nfe,&iflag,awork,bwork,iwork);

  if (iflag == -2) {
    fprintf(stderr,"error: invalid input in simplex\n");
    exit(1);
  }
  for (i=1;i<=n;i++) {
    x[i]=_x[i-1];      
  }

  *fx = _fx;
  return iflag;
}







