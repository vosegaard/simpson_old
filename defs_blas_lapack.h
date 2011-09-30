#ifndef __DEFS_BLAS_LAPACK_H
#define __DEFS_BLAS_LAPACK_H



extern zdotusub_(const int *,const void *,const int *,const void *,const int *, void *);
extern zdotcsub_(const int *,const void *,const int *,const void *,const int *, void *);
extern zaxpy_(const int *, const void *, const void *, const int *,const void *, const int *);
extern zscal_(const int *, const void *, const void *, const int *);
extern zdscal_(const int *, const double *, const void *, const int *);
extern zgemm_(char *, char *, const int *, const int *, const int *, const void *, 
              const void *, const int *, const void *, const int *, const void *,
	      const void *, const int *);
extern zgemv_(char *,const int *,const int *,const void *,const void *,const int *,
              const void *,const int *,const void *,const void *,const int *);
extern zgerc_(const int *,const int *,const void *,const void *,const int *,const void *,
              const int *,const void *,const int *);
extern zhemm_(char *,char *,const int *,const int *,const void *,const void *,const int *,
              const void *,const int *,const void *,const void *, const int *);
extern dscal_(const int *, const double *, const double *, const int *); 
extern dsyev_(char *,char *,const int *,const double *,const int *, const double *,const double *,
              const int *, const int *);
extern zlarcm_(const int *,const int *,const double *,const int *, const void *,const int *,
               const void *,const int *, const double *);
extern dsymm_(char *,char *,const int *,const int *,const double*,const double *,
              const int *,const double *,const int *,const double *,const double *,const int *);
extern dgemm_(char *,char *,const int *,const int *,const int *,const double *,const double *,
              const int *,const double *,const int *,const double *,const double *,const int *);
extern zgesv_(const int *,const int *,const void *,const int *,const long int *,const void *,
              const int *, const int *);
extern zlaset_(char *,const int *,const int *,const void *,const void *,const void *,const int *);
extern daxpy_(const int *,const double *,const double *,const int *,const double *,const int *);
extern dasumsub_(const int *,const double *,const int *,const double *);
extern dznrm2sub_(const int *,const void *,const int *,const double *);
extern daxpy_(const int *, const double *, const double *, const int *,const double *, const int *);
extern dcopy_(const int *,const double *,const int *,const double *,const int *);



extern dlansy_(char *,char *,const int *,const double *,const int *,const double *);
extern zlanhe_(char *,char *,const int *,const void *,const int *,const double *);









#endif
