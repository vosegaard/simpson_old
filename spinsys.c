/*
    Spinsystem routines
    Copyright (C) 1999 Mads Bak

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
    
    Routines for calculation of spin-tensors for the spin-system.
    Used by readsys.c that creates the spinsystem, and 
    sim.c to create the start and detect operators.
*/

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <tcl.h>
#include <math.h>
#include "matrix_new.h"
#include "cm_new.h"
#include "spinsys.h"
#include "tclutil.h"
#include "ham.h"

double ss_gamma1H()
{
   if (isotopes[0].number != 1) {
     fprintf(stderr,"error: isotope list does not have 1H as first entry\n");
     exit(1);
   }
   return isotopes[0].gamma;
}

void ss_showspins(SpinSys* S)
{
  int i;
  
  for (i=1;i<=S->nspins;i++) {
    printf("Spin[%d] : I = %g\n",i,S->iso[i]->spin);
  }
}

ISOTOPE* ss_findisotope(char* name)
{
  char buf[256],*src,*dst;
  int number,niso, gotnumber=0;
  
  if (!name) {
    fprintf(stderr,"ss_findisotope: argument must not be NULL\n");
    exit(1);
  }
  src=name;
  dst=buf;
    
  while (isdigit(*src) && *src != 0) {
    *dst++ = *src++;
  }

  if (src == name) {
    fprintf(stderr,"error: name of isotope '%s' is wrong, '23Na' is a correct example\n",name);
    exit(1);
  }

  *dst=0;
  number=atoi(buf);

  niso=0;
  while (isotopes[niso].number) {
    if (isotopes[niso].number == number) {
      gotnumber=1;
      sprintf(buf,"%d%s",isotopes[niso].number,isotopes[niso].name);
      if (!strcmp(name,buf)) {
/*
        fprintf(stderr,"error: name of isotope '%s' is wrong, maybe '%s' is correct ?\n",name,buf);
        exit(1);
*/
        return &isotopes[niso];   
      }
    }
    niso++;
  }
  fprintf(stderr,"error: unable to find name '%s' in the internal isotope table.\n", name);
  if (gotnumber==1) {
    niso=0;
    fprintf(stderr, "       My suggestions are ");
    while (isotopes[niso].number) {
      if (isotopes[niso].number == number)
        fprintf(stderr, "'%d%s' ", number, isotopes[niso].name);
      if (!strcmp(isotopes[niso].name, src))
        fprintf(stderr, "'%d%s' ", isotopes[niso].number, isotopes[niso].name);
      niso++;
    }
    fprintf(stderr, "\n");
  }
  exit(1);
}

double ss_qn(SpinSys* S,int spin)
{ 
  if (spin < 1 || spin > S->nspins) {
    fprintf(stderr,"error: invalid range of spin (was %d)\n",spin);
    exit(1);
  }
  return S->iso[spin]->spin;
}

double ss_gamma(SpinSys* S,int spin)
{ 
  if (spin < 1 || spin > S->nspins) {
    fprintf(stderr,"error: invalid range of spin (was %d)\n",spin);
    exit(1);
  }
  return S->iso[spin]->gamma;
}

ISOTOPE* ss_isotope(SpinSys* S,int spin)
{ 
  if (spin < 1 || spin > S->nspins) {
    fprintf(stderr,"error: invalid range of spin (was %d)\n",spin);
    exit(1);
  }
  return S->iso[spin];
}

int ss_matdim(SpinSys* S)
{
  int i;
  double N;

  if (S->nspins < 1) return 0;
  N = 1.0;
  for (i=1;i<=S->nspins;i++) 
    N *= (double)(2.0*(S->iso[i]->spin)+1.0);
  return (int)N;
}


void ss_addspin(SpinSys* S,char* name)
{
  S->nspins++;

  if (S->nspins >= MAXSPINS) {
    fprintf(stderr,"error: max number of spins reached, increase MAXSPINS\n");
    exit(1);
  }
  S->iso[S->nspins]=ss_findisotope(name);  
  S->matdim=ss_matdim(S);
}

void ss_initialize(SpinSys* S)
{
  S->nspins=0;
  S->nchan=0;
}

int ss_issame(SpinSys* S,int spin1,int spin2)
{
  if (spin1 < 1 || spin1 > S->nspins) {
    fprintf(stderr,"error: invalid range of spin (was %d)\n",spin1);
    exit(1);
  }
  if (spin2 < 1 || spin2 > S->nspins) {
    fprintf(stderr,"error: invalid range of spin (was %d)\n",spin2);
    exit(1);
  }
  return (S->iso[spin1] == S->iso[spin2]);
}

mv_complx * _Ip(double I)
{
  int i,j,S;  
  double mm,m;
  mv_complx * M;
  
  S= (int)(2*I+1);
  M=complx_matrix_alloc(S,S);
  for (i=0;i<S;i++) {
    mm= I-i;
    for (j=0;j<S;j++) {
      m= I-j;
      if (mm == m + 1)
        M->data[i+j*S]=Complx(sqrt(I*(I+1.0)-m*(m+1.0)),0.0);
      else
        M->data[i+j*S]=Complx(0.0,0.0);
    }
  }
  return M;
}

mv_complx * _Im(double I)
{
  int i,j,S;  
  double mm,m;
  mv_complx * M;
  
  S= (int)(2*I+1);
  M=complx_matrix_alloc(S,S);
  for (i=0;i<S;i++) {
    mm= I-i;
    for (j=0;j<S;j++) {
      m= I-j;
      if (mm == m - 1)
        M->data[i+j*S]=Complx(sqrt(I*(I+1.0)-m*(m-1.0)),0.0);
      else
        M->data[i+j*S]=Complx(0.0,0.0);
    }
  }
  return M;
}

mv_complx * _Iz(double I)
{
  int i,j,S;  
  double mm,m;
  mv_complx * M;
  
  S= (int)(2*I+1);
  M=complx_matrix_alloc(S,S);
  for (i=0;i<S;i++) {
    mm= I-i;
    for (j=0;j<S;j++) {
      m= I-j;
      if (mm == m)
        M->data[i+j*S]=Complx(m,0.0);
      else
        M->data[i+j*S]=Complx(0.0,0.0);
    }
  }
  return M;
}

mv_complx * _Iq(double I)
{
  int i,j,S;  
  double mm,m;
  mv_complx * M;
  
  S= (int)(2*I+1);
  M=complx_matrix_alloc(S,S);
  for (i=0;i<S;i++) {
    mm= I-i;
    for (j=0;j<S;j++) {
      m= I-j;
      M->data[i+j*S]=Complx(m,mm);
    }
  }
  return M;
}

/* this probably wont work in multi-spin cases. */
mv_complx * _Ic(double I)
{
  int S;
  mv_complx * M;
  /* ZT: adding test for I>1 and half-integer spin */
  if ( (I<1.4) || ( fabs(I-floor(I)-0.5)>0.01 ) ) {
     fprintf(stderr,"error in central transition operator, spin=%f is not allowed\n",I);
     exit(1);
  }

  S= (int)(2*I+1);
  M=complx_matrix_alloc(S,S);
  cmv_zero(M);
  M->data[S/2-1+(S/2)*S]=Complx(sqrt(I*(I+1)+0.25),0.0);
  return M;
}

mv_complx * _Ix(double I)
{
/* return (Ip(I)+Im(I))*complx(0.5,0); */
  mv_complx *ip,*im;

  ip=_Ip(I);
  im=_Im(I);
  cmv_addto(ip,im);
  cmv_muld(ip,0.5);
  complx_matrix_free(im);
  return ip;
}

mv_complx * _Iy(double I)
{
/* return (Ip(I)-Im(I))*complx(0.0,-5.0); */
  mv_complx *ip,*im;

  ip=_Ip(I);
  im=_Im(I);
  cmv_subfrom(ip,im);
  cmv_mulc(ip,Complx(0.0,-0.5));
  complx_matrix_free(im);
  return ip;
}

mv_complx * _Ie(double I)
{
  int S;
  mv_complx * M;

  S= (int)(2*I+1);
  M=complx_matrix_alloc(S,S);
  cm_unit(M);
  return M;
}

/* creates the direct product, for example creates
  for a four spin-system with spin=3 the matrix:
    1 (x) 1 (x) m (x) 1
*/
mv_complx * In(SpinSys* S,mv_complx * m,int spin)
{
  int i;
  mv_complx *curr=NULL,*tmp,*tmp2;
  
  for (i=1;i <= S->nspins;i++) {
    if (i == spin) {
      if (i == 1) {
        curr=cmv_dup(m);
      } else {
        tmp=cm_direct(curr,m);
        complx_matrix_free(curr);
        curr=tmp;
      }
    } else {
      if (i == 1) {
        curr=_Ie(ss_qn(S,i)); 
      } else {
        tmp2=_Ie(ss_qn(S,i));
        tmp=cm_direct(curr,tmp2);
        complx_matrix_free(curr);
        complx_matrix_free(tmp2);
        curr=tmp;
      }
    }
  }
  return curr;
}

/****
 *  ZT: new and faster version of In (not using direct product)
 ****/
mv_complx * In_new(SpinSys* S, mv_complx * m,int spin)
{
	  int k,l,i,j,r,c,K,L,N,matdim,nnz=0;
	  complx z;
	  mv_complx *res;

	  /*printf("spin = %d\n",spin);
	  cm_print(m,"matrix is");*/
	  N = m->row;
	  matdim = S->matdim;
	  K = 1;
	  for (i=1;i<spin;i++) {
	     K *= (int)(2.0*ss_qn(S,i)+1.0);
	  }
	  L = 1;
	  for (i=spin+1; i<=S->nspins; i++) {
	     L *= (int)(2.0*ss_qn(S,i)+1.0);
	  }
	  /*printf("sp_In report: K = %d; N = %d; L = %d; matdim = %d\n",K,N,L,matdim);*/
	  if (K*N*L != matdim) {
	     fprintf(stderr,"sp_In error: mismatch in dimensions (%d * %d * %d != %d)\n",K,N,L,matdim);
	     exit(1);
	  }
	  res = complx_matrix_alloc(matdim,matdim);
	  cmv_zero(res);

	  for (i=0; i<N; i++) {
	     for (j=0; j<N; j++) {
	        z = m->data[i+j*N];
		if ( fabs(z.re) < 1e-6 && fabs(z.im) < 1e-6) continue;
		for (k=0; k<K; k++) {
		   for (l=0; l<L; l++) {
		      r = i*L +l + k*L*N;
		      c = j*L +l + k*L*N;
		      res->data[r+c*matdim] = z;
		      nnz++;
		      /*printf("[%d, %d] = (%f, %f)\n",r+1,c+1,z.re,z.im);*/
		   }
		}
	     }
	  }
	  /*cm_print(res,"result");
	  exit(1); */

	  return res;
}


mv_complx * Icoherence(SpinSys* S, double* coh)
{
  int i,j,k,N;  
  mv_complx *m=NULL,*curr=NULL,*tmp=NULL;

  for (k=1;k <= S->nspins;k++) {

    m = _Iq(ss_qn(S,k));
    N=m->row;
    for (i=0;i<N;i++) {    
      for (j=0;j<N;j++) {
        double ok = 0;
	int pp = i+j*N;
        if ( (m->data[pp].im - m->data[pp].re) == coh[k] ) ok=1;
        m->data[pp]=Complx(ok,0.0);        
      }
    }
    if (k == 1) {
      curr=m; 
    } else {
      tmp=cm_direct(curr,m);
      complx_matrix_free(curr);
      complx_matrix_free(m);
      curr=tmp;      
    }
  }
  return curr;
}


mv_complx * Iq(SpinSys* S)
{
  int i;
  mv_complx *curr=NULL,*tmp,*tmp2;

  for (i=1;i <= S->nspins;i++) {
    if (i == 1) {
      curr=_Iq(ss_qn(S,i)); 
    } else {
      tmp2=_Iq(ss_qn(S,i));
      tmp=cm_directadd(curr,tmp2);
      complx_matrix_free(curr);
      complx_matrix_free(tmp2);
      curr=tmp;      
    }
  }
  return curr;
}



mv_complx* Iqdelta(SpinSys* S)
{
  int i,j,N;
  mv_complx * M;
  
  N=S->matdim;

  M=Iq(S);

  for (i=0;i<N;i++) {
    for (j=0;j<N;j++) {
      int pp = i+j*N;
      M->data[pp] = Complx( M->data[pp].im - M->data[pp].re, 0.0);
    }
  }  
  return M;
}

mv_complx * Ie(SpinSys* S)
{
  mv_complx *a;

  a=complx_matrix_alloc(S->matdim,S->matdim);
  cm_unit(a);
  return a;
}

/****
 * ZT: adding I_alpha operator for spin=1/2
 ****/
mv_complx * Ia(SpinSys* S, int spin) {
  mv_complx *a, *b;
  double I;

  I = ss_qn(S,spin);
  if ( fabs(I-0.5)>0.01 ) {
     fprintf(stderr,"error in alpha-state operator, it is undefined for spin=%f\n",I);
     exit(1);
  }
  a = complx_matrix_alloc(2,2);
  cmv_zero(a);
  a->data[0].re = 1.0;
  b = In_new(S,a,spin);
  complx_matrix_free(a);
  return b;
}

/****
 * ZT: adding I_beta operator for spin=1/2
 ****/
mv_complx * Ib(SpinSys* S, int spin) {
  mv_complx *a, *b;
  double I;

  I = ss_qn(S,spin);
  if ( fabs(I-0.5)>0.01 ) {
     fprintf(stderr,"error in alpha-state operator, it is undefined for spin=%f\n",I);
     exit(1);
  }
  a = complx_matrix_alloc(2,2);
  cmv_zero(a);
  a->data[3].re = 1.0;
  b = In_new(S,a,spin);
  complx_matrix_free(a);
  return b;
}

#define DECLARE_SIFUNC(TYPE) \
mv_complx * TYPE(SpinSys* S,int spin)\
{\
  mv_complx *a,*b;\
  a= _##TYPE(ss_qn(S,spin));\
  b=In_new(S,a,spin);\
  complx_matrix_free(a);\
  return b;\
}

DECLARE_SIFUNC(Ic)
DECLARE_SIFUNC(Ip)
DECLARE_SIFUNC(Im)
DECLARE_SIFUNC(Ix)
DECLARE_SIFUNC(Iy)
DECLARE_SIFUNC(Iz)


#ifdef DEBUG
#define DEBUG_PARSER(x) printf(x)
#define DEBUG_PARSER2(x,y) printf(x,y)
#else
#define DEBUG_PARSER(x)
#define DEBUG_PARSER2(x,y)
#endif

mv_complx * ss_oper(SpinSys* S,char* name)
{
  int i,spin,n,N;
  char buf[16],*pname,*pbuf;
  mv_complx *sum,*tmp;
  
  pname=name;
  if (*pname != 'I') {
     fprintf(stderr,"operator `%s` not known, must be of type "
             "Ina, where n=1,2,.., or n='n' for sum of all spins and "      
             "a=x,y,z,p,m\n",name);
     exit(1);
  }
  pname++;

  if (*pname == 'n') {
    pname++;

    n=S->nspins;
    N=S->matdim;
    sum=complx_matrix_alloc(N,N);
    cmv_zero(sum);
    if (*(pname+1) == 'c') {
      /* The syntax I2xc is not properply implemented yet. */
      fprintf(stderr,"oper: unknown operator: %s\n",name);
      exit(1);
      pname++;
    } else {
      switch (*pname) {
      case 'x': for (i=1;i<=n;i++) { tmp=Ix(S,i); cmv_addto(sum,tmp); complx_matrix_free(tmp); }; break;
      case 'y': for (i=1;i<=n;i++) { tmp=Iy(S,i); cmv_addto(sum,tmp); complx_matrix_free(tmp); }; break;
      case 'z': for (i=1;i<=n;i++) { tmp=Iz(S,i); cmv_addto(sum,tmp); complx_matrix_free(tmp); }; break;
      case 'p': for (i=1;i<=n;i++) { tmp=Ip(S,i); cmv_addto(sum,tmp); complx_matrix_free(tmp); }; break;
      case 'm': for (i=1;i<=n;i++) { tmp=Im(S,i); cmv_addto(sum,tmp); complx_matrix_free(tmp); }; break;
      case 'c': for (i=1;i<=n;i++) { tmp=Ic(S,i); cmv_addto(sum,tmp); complx_matrix_free(tmp); }; break;
      case 'a': for (i=1;i<=n;i++) { tmp=Ia(S,i); cmv_addto(sum,tmp); complx_matrix_free(tmp); }; break;
      case 'b': for (i=1;i<=n;i++) { tmp=Ib(S,i); cmv_addto(sum,tmp); complx_matrix_free(tmp); }; break;
       default:
         fprintf(stderr,"oper: unknown operator: '%s'\n",name);
         exit(1);
      }
    }
    pname++;
    if (*pname != 0) {
       fprintf(stderr,"oper: trailing garbage at operator: '%s'\n",name);
       exit(1);
    }
    return sum;

  } else {
    if (!isdigit(*pname) || *pname == '0') {
       fprintf(stderr,"oper: unknown operator: '%s'\n",name);
       exit(1);
    }
    pbuf=buf;
    *pbuf++ = *pname++;
    while (isdigit(*pname)) {
      *pbuf++ = *pname++;
      if (pname-name > 6) {
         fprintf(stderr,"oper: unknown operator: '%s'\n",name);
         exit(1);
      }
    }
    *pbuf=0;
    spin=atoi(buf);
    
    if (*(pname+1) == 'c') {
      /* The syntax I2xc is not properply implemented yet. */
      fprintf(stderr,"oper: unknown operator: %s\n",name);
      exit(1);
      pname++;
    } else {
      switch (*pname) {
      case 'x': return Ix(S,spin); break;
      case 'y': return Iy(S,spin); break;
      case 'z': return Iz(S,spin); break;
      case 'p': return Ip(S,spin); break;
      case 'm': return Im(S,spin); break;
      case 'c': return Ic(S,spin); break;
      case 'a': return Ia(S,spin); break;
      case 'b': return Ib(S,spin); break;
      default:
         fprintf(stderr,"oper: unknown operator: %s\n",name);
         exit(1);
      }
    }
    if (*(pname+1) != 0) {
       fprintf(stderr,"oper: trailing garbage at operator: '%s'\n",name);
       exit(1);
    }

  }
  return NULL;
}


/* A simple hash-list for the spin-operator parser 
   Currently only a list.
*/

#define HASH_MAXN 128

int _hash_n=0;
char* _hash_names[HASH_MAXN];
complx _hash_complx[HASH_MAXN];

void ss_hashinit()
{
  _hash_n=0;
}

void ss_hashinsert(char* ptr,complx c)
{
  printf("insert %s %g %g\n",ptr,c.re,c.im);
  _hash_names[_hash_n]=strdup(ptr);  
  _hash_complx[_hash_n]=c;
  _hash_n++;
}

int ss_hashlookup(char* ptr,complx* c)
{
  int i;
  for (i=0;i<_hash_n;i++) {
     if (!strcmp(_hash_names[i],ptr)) {
        *c = _hash_complx[i];
        return 1;
     }
  } 
  return 0;
}

void ss_hashdestroy()
{
  int i;
  for (i=0;i<_hash_n;i++) {
     free(_hash_names[i]);
  } 
}

/*
 Description:  A parser with the following grammar

 Grammar:
 ----------------------------
 program:
	expression RETURN
	
 expression:
	expression + term
	expression - term
	term
	
 term:
	term * primary
	primary
	
 primary:
	SPINOPERATOR
	NUMBER
	IMAG
    -primary
	( expression )	
 ----------------------------
*/ 
/* Defines of tokens */
#define	sop_token_SPINOPERATOR 1
#define	sop_token_NUMBER 2
#define	sop_token_IMAG 3
#define	sop_token_LP 4
#define	sop_token_RP 5

#define	sop_token_MINUS 6
#define	sop_token_MULT 7
#define	sop_token_PLUS 8
#define	sop_token_RETURN 9

char      sop_initial[256];
complx    sop_complexval;
mv_complx*  sop_operator;
int       sop_current_Token;
int       sop_matrix_dimension;
int       sop_nspins;
char*     sop_current_chrptr;
SpinSys*  sop_spinsysptr;

int sop_get_token()
{
	char *chr,*chr0;
	char buf[100],*bufp;
	
	chr = sop_current_chrptr;
	chr0 = chr;

	while (*chr==' ') {
		chr++;
	}

	switch (*chr)
	{
	case '\0':
		return sop_token_RETURN;
	case '+':
		DEBUG_PARSER("PLUS\n");
		sop_current_chrptr = ++chr;
		return sop_token_PLUS;
	case '*':
		DEBUG_PARSER("MULT\n");
		sop_current_chrptr = ++chr;
		return sop_token_MULT;
	case '-':
		DEBUG_PARSER("MINUS\n");
		sop_current_chrptr = ++chr;
		return sop_token_MINUS;
	case 'i':
		DEBUG_PARSER("IMAG\n");
		sop_current_chrptr = ++chr;
		return sop_token_IMAG;
	case 'I':
		bufp=buf;
		*bufp++ = *chr++;
		if (*chr=='n') {
			*bufp++ = *chr++;
		} else if (isdigit(*chr) && *chr != '0') {
			*bufp++ = *chr++;
			while (isdigit(*chr)) {
				*bufp++ = *chr++;
			}
		} else {
  		fprintf(stderr,"Cannot parse spinoperator at character"
                     " '%s' in operator '%s'\n",chr,sop_initial);
			exit(-1);
		}
		switch (*chr)
		{
		case 'x': *bufp ='x';break;
		case 'y': *bufp ='y';break;
		case 'z': *bufp ='z';break;
		case 'p': *bufp ='p';break;
		case 'm': *bufp ='m';break;
		case 'c': *bufp ='c';break;
		case 'a': *bufp ='a';break;
		case 'b': *bufp ='b';break;
    default:
  		fprintf(stderr,"Cannot parse spinoperator at character"
                     " '%s' in operator '%s'\n",chr,sop_initial);
			exit(-1);
		}
		bufp++;
    if (*(chr+1) == 'c') {    
      *bufp++ = *(++chr);
    }
		*bufp=0;
				
		DEBUG_PARSER2("SPINOPERATOR: %s\n",buf);
		sop_current_chrptr = ++chr;
		sop_operator = ss_oper(sop_spinsysptr,buf);
		return sop_token_SPINOPERATOR;
		case '(':
			DEBUG_PARSER2("LP : %c\n",*chr);
			sop_current_chrptr = ++chr;
			return sop_token_LP;
		case ')':
			DEBUG_PARSER2("RP : %c\n",*chr);
			sop_current_chrptr = ++chr;
			return sop_token_RP;
		case '1': case '2': case '3': case '4': case '5':
		case '6': case '7': case '8': case '9': case '0': case '.':
			bufp = buf;
			if (*chr != '.') {
				*bufp++ = *chr++;
				while (isdigit(*chr))
					*bufp++ = *chr++;
			}
			if (*chr == '.') {
				*bufp++ = *chr++;
				if (!isdigit(*chr)) {
  		    fprintf(stderr,"Cannot parse spinoperator at character"
                         " '%s' in operator '%s'\n",chr,sop_initial);
					exit(-1);
				}
				while (isdigit(*chr))
					*bufp++ = *chr++;
			}
			if (*chr == 'E' || *chr == 'e') {
				*bufp++ = *chr++;
				if (*chr == '+' || *chr == '-')
					*bufp++ = *chr++;
				if (!isdigit(*chr)) {
  		    fprintf(stderr,"Cannot parse spinoperator at character"
                         " '%s' in operator '%s'\n",chr,sop_initial);
					exit(-1);
				}
				while (isdigit(*chr))
					*bufp++ = *chr++;
			}
			if (!isdigit(*chr)) chr--;
			*bufp=0;
			DEBUG_PARSER2("NUMBER : %s\n",buf);
			bufp=buf;
			sop_complexval.re = strtod(buf,&bufp);
      sop_complexval.im = 0.0;
			sop_current_chrptr = ++chr;
			return sop_token_NUMBER;
		default:
      if (isalpha(*chr) || *chr == '_') {
			  bufp = buf;
			  while (isalnum(*chr) || *chr == '_') {
				  *bufp++ = *chr++;
			  }
        *bufp=0;
        if (!ss_hashlookup(buf,&sop_complexval)) {
  		    fprintf(stderr,"Unknown variable '%s' in operator '%s'\n",buf,sop_initial);
			    exit(-1);
        }
        sop_current_chrptr = chr;
        return sop_token_NUMBER;
      } else {
  		  fprintf(stderr,"Cannot parse spinoperator at character"
                     " '%s' in operator '%s'\n",chr,sop_initial);
			  exit(-1);
      }
			return sop_token_RETURN;
	}
}

mv_complx * sop_term(int get);

mv_complx * sop_expression(int get)
{
	mv_complx * M;
	mv_complx * left;

	left = sop_term(get);

	for (;;)
	{
		switch (sop_current_Token)
		{
		case sop_token_PLUS:
			M = sop_term(1);
			cmv_addto(left,M);
			complx_matrix_free(M);
			break;
		case sop_token_MINUS:
			M = sop_term(1);
			cmv_subfrom(left,M);
			complx_matrix_free(M);
			break;
		default:
			return left;
		}
	}
}

mv_complx * sop_prim(int get);

mv_complx * sop_term(int get)
{
	mv_complx* M;
	mv_complx* M2;
	mv_complx* left;
	
	left = sop_prim(get);

	for (;;)
	{
		switch (sop_current_Token)
		{
		case sop_token_MULT:
			M = complx_matrix_alloc(sop_matrix_dimension,sop_matrix_dimension);

			M2 = sop_prim(1);
			cm_mul(M,left,M2);
			complx_matrix_free(M2);
			complx_matrix_free(left);
			left = M;
			break;
		default:
			return left;
		}
	}
}

mv_complx * sop_prim(int get)
{
	int i;
	mv_complx * M;

	if (get) 
	{
		sop_current_Token = sop_get_token();
	}

	switch (sop_current_Token)
	{
	case sop_token_SPINOPERATOR:
		M = complx_matrix_alloc(sop_matrix_dimension,sop_matrix_dimension);
		cmv_copy(M,sop_operator);
		complx_matrix_free(sop_operator);
		sop_current_Token = sop_get_token();
		if(sop_current_Token < sop_token_RP)
		{
			fprintf(stderr,"Parse error: primary not expected in operator '%s'\n",sop_initial);
			exit(1);
		}
		return M;

	case sop_token_NUMBER:
		M = complx_matrix_alloc(sop_matrix_dimension,sop_matrix_dimension);
		cmv_zero(M);
		for (i=0;i<sop_matrix_dimension;i++)
			M->data[i+i*sop_matrix_dimension] = sop_complexval;
		sop_current_Token = sop_get_token();
		if(sop_current_Token < sop_token_RP)
		{
			fprintf(stderr,"Parse error: primary not expected in operator '%s'\n",sop_initial);
			exit(1);
		}
		return M;

	case sop_token_IMAG:
		M = complx_matrix_alloc(sop_matrix_dimension,sop_matrix_dimension);
		cmv_zero(M);
		for (i=0;i<sop_matrix_dimension;i++)
			M->data[i+i*sop_matrix_dimension] = Complx(0.0,1.0);
		sop_current_Token = sop_get_token();
		if(sop_current_Token < sop_token_RP)
		{
			fprintf(stderr,"Parse error: primary not expected in operator '%s'\n",sop_initial);
			exit(1);
		}
		return M;

	case sop_token_RP:
	  fprintf(stderr,"Parse error: \")\" not expected in operator '%s'\n",sop_initial);
	  exit(1);
    return 0;
    
	case sop_token_MINUS:
		M = sop_prim(1);
		cmv_muld(M,-1.0);
		return M;

	case sop_token_LP:
		M = sop_expression(1);
		if (sop_current_Token != sop_token_RP) 
		{
		  fprintf(stderr,"Parse error: expected \")\" in operator '%s'\n",sop_initial);
			exit(1);
		}
		sop_current_Token = sop_get_token();
		if(sop_current_Token < 6)
		{
			fprintf(stderr,"Parse error: primary not expected in operator '%s'\n",sop_initial);
			exit(1);
		}
		return M;
	default:
	  fprintf(stderr,"Parse error: primary expected in operator '%s'\n",sop_initial);
		exit(1);
	}
	return 0;
}

mv_complx* ss_readoper(SpinSys * S,char* psop)
{  
	sop_spinsysptr = S;
	sop_matrix_dimension = S->matdim;
    sop_nspins = S->nspins;

	sop_current_chrptr = psop;
  strcpy(sop_initial,sop_current_chrptr);
	sop_current_Token = sop_get_token();
	if (sop_current_Token==sop_token_RETURN)
	{
		fprintf(stderr,"Parse error: operator is empty\n");
		exit(1);
	}
	sop_operator = sop_expression(0);
	if (sop_current_Token != sop_token_RETURN)
	{
		fprintf(stderr,"Parse error: invalid syntax in '%s'\n",sop_initial);
		exit(1);
	}

#if DEBUG
	cm_print(sop_operator,"testmatrix");
#endif
	return sop_operator;
}

/* 
END OF PARSERCODE
  */
  
  
/******
 * ZT: this creates Boltzman equilibrium density matrix for a given spin system
 *     High temperature limit, normalized according to 1H polarization
 *****/
mv_complx * ss_eq_sigma(SpinSys* S)
{
   int n, i;
   mv_complx *mx, *sum;
   double b;
   
   n = S->nspins;
   sum = complx_matrix_alloc(S->matdim,S->matdim);
   cmv_zero(sum);
   for (i=1; i<=n; i++) {
      mx = Iz(S,i);
      b = ss_gamma(S,i)/ss_gamma1H();
      cmv_multod(sum,mx,b);
      complx_matrix_free(mx);
   }
   return sum;
}

int tclIsotopes(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  char buf[256];
  ISOTOPE* ptr;
  if (argc != 1) {
    interp->result = "Usage: isotopes  (returns a list of isotopes: number, nuclei, spin, and gyromag ratio)";
    return TCL_ERROR;
  }

  Tcl_ResetResult(interp);
  ptr=isotopes;
  while (ptr->number != 0) {
     sprintf(buf,"%d %s %g %g",ptr->number,ptr->name,ptr->spin,ptr->gamma);
     Tcl_AppendElement(interp,buf);
     ptr++;
  }
  return TCL_OK;
}    

#define HBAR 1054.59198

int tclDist2Dip(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  char* nuc1,*nuc2;
  double dist,dip,gamma1,gamma2;
  ISOTOPE* isop;

  if (argc != 4) {
    interp->result = "dist2dip <nuc1> <nuc2> <distance> : calculates dipolar coupling (Hz) from internuclear distance (Aangstrom)";
    return TCL_ERROR;
  }
  nuc1=argv[1];
  nuc2=argv[2];
  if (Tcl_GetDouble(interp,argv[3],&dist) != TCL_OK)
    return TCL_ERROR;   

  isop=ss_findisotope(nuc1);
  gamma1=isop->gamma;
  isop=ss_findisotope(nuc2);
  gamma2=isop->gamma;
 
  dip= -gamma1*gamma2*HBAR/(2.0*M_PI*dist*dist*dist);

  sprintf(interp->result,"%g",dip);
  return TCL_OK;
}

int tclDip2Dist(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  char* nuc1,*nuc2;
  double dist,dip,gamma1,gamma2;
  ISOTOPE* isop;

  if (argc != 4)
    return TclError(interp,"dip2dist <nuc1> <nuc2> <dipolar coupling> :"
    " calculates internuclear distance (Aangstrom) from dipolar coupling (Hz)");

  nuc1=argv[1];
  nuc2=argv[2];
  if (Tcl_GetDouble(interp,argv[3],&dip) != TCL_OK)
    return TCL_ERROR;   

  isop=ss_findisotope(nuc1);
  gamma1=isop->gamma;
  isop=ss_findisotope(nuc2);
  gamma2=isop->gamma;
  if (-dip/gamma1*gamma2 < 0.0)
     return TclError(interp,"dip2dist: the dipolar coupling for nuclei %s and %s must be %s",nuc1,nuc2,
           (gamma1*gamma2 < 0.0 ? "positive" : "negative"));

  dist=pow(fabs(gamma1*gamma2*HBAR/(M_PI*2.0*dip)),1.0/3.0);  
  sprintf(interp->result,"%g",dist);
  return TCL_OK;
}

int tclGamma(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  ISOTOPE* isop;

  if (argc != 2)
    return TclError(interp,
    "gamma <nuc> : returns the magnetogyric ratio of the nucleus in the unit 10^7 rad/(T s)");
    
  isop=ss_findisotope(argv[1]);
  sprintf(interp->result,"%g",isop->gamma);
  return TCL_OK;
}


int tclResfreq(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  double protonfield;
  ISOTOPE* isop;

  if (argc != 2 && argc != 3)
    return TclError(interp,
     "resfreq <nuc> ?<proton_frequency/Hz default= 1e6Hz>?: returns the absolute resonance\n"
     " frequency in Hz for a nucleus optionally at a given proton frequency");
    
  isop=ss_findisotope(argv[1]);
  protonfield=1.0;
  if (argc == 3) {
    if (Tcl_GetDouble(interp,argv[2],&protonfield) != TCL_OK)
      return TCL_ERROR;   
    if (protonfield < 10000.0)
      return TclError(interp,"resfreq: illegal value of proton frequency (must be positive and > 10000Hz)\n");
  }
  isop=ss_findisotope(argv[1]);
  sprintf(interp->result,"%g", fabs(isop->gamma/ss_gamma1H()*protonfield));
  return TCL_OK;
}


void tclcmd_spinsys(Tcl_Interp* interp)
{

  Tcl_CreateCommand(interp,"isotopes",tclIsotopes,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"gamma",tclGamma,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"resfreq",tclResfreq,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"dist2dip",tclDist2Dip,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"dip2dist",tclDip2Dist,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

}

              



