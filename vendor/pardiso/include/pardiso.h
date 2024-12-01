/* Copyright 2023, Panua Technologies, Switzerland.
 *
 * All rights reserved.
 *
 * This header file defines all functions that can be called in the Pardiso
 * library.
 */

/* TODO: Here should be also the documentation of the arguments for the functions. */

/* Structure for storing double precision complex numbers. */
typedef struct{
	double re; 
	double i;}
doublecomplex;

/*
#ifdef PARDISO_COMPLEX
typedef p_double p_doublecomplex;
#else
typedef p_double double;
#endif
*/

#ifdef __cplusplus
extern "C"
{
#endif

#define DllImport   __declspec( dllimport )

DllImport void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
// "void" in the following means that it can be double or doublecomplex
DllImport void pardiso     (void   *, int    *,   int *, int *,    int *, int *,
                  void   *, int    *,    int *, int *,   int *, int *,
                     int *, void    *, void   *, int *, double *);
DllImport void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
DllImport void pardiso_chkvec     (int *, int *, double *, int *);
DllImport void pardiso_printstats (int *, int *, double *, int *, int *, int *,
                           double *, int *);

DllImport void pardisoinit_z (void   *, int    *,   int *, int *, double *, int *);
// "void" in the following means that it can be double or doublecomplex
DllImport void pardiso_z     (void   *, int    *,   int *, int *,    int *, int *,
                  void   *, int    *,    int *, int *,   int *, int *,
                     int *, void    *, void   *, int *, double *);
DllImport void pardiso_chkmatrix_z  (int *, int *, void *, int *, int *, int *);
DllImport void pardiso_chkvec_z     (int *, int *, void *, int *);
DllImport void pardiso_printstats_z (int *, int *, void *, int *, int *, int *,
                           void *, int *);
DllImport void pardiso_get_schur_z(void*, int*, int*, int*, void*, int*, int*);

DllImport void pardisoinit_d (void   *, int    *,   int *, int *, double *, int *);
// "void" in the following means that it can be double or doublecomplex
DllImport void pardiso_d     (void   *, int    *,   int *, int *,    int *, int *,
                  void   *, int    *,    int *, int *,   int *, int *,
                     int *, void    *, void   *, int *, double *);
DllImport void pardiso_chkmatrix_d  (int *, int *, void *, int *, int *, int *);
DllImport void pardiso_chkvec_d     (int *, int *, void *, int *);
DllImport void pardiso_printstats_d (int *, int *, void *, int *, int *, int *,
                           void *, int *);
DllImport void pardiso_get_schur_d(void*, int*, int*, int*, void*, int*, int*);
  
#ifdef __cplusplus
}
#endif
