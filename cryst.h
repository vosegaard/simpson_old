#ifndef __CRYST_H
#define __CRYST_H

typedef struct _CRYSTALLITE {
  double alpha,beta,weight;
} CRYSTALLITE;

extern char* cryst_names[];
extern int cryst_numbers[];
extern CRYSTALLITE* cryst_pointers[];

typedef struct _TRIANGLE {
  int a, b, c;
  double weight;
} TRIANGLE;

typedef struct _Omega {
   double alpha, beta, gamma;
} Omega;

#endif
