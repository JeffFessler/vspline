/*
* defs-vspline.h
*
* Definitions for vspline that may need to be changed
* to make the code more transportable.
* Copyright 1989, Jeff Fessler
*/
#ifndef DEFS_VSPLINE
#define DEFS_VSPLINE

#ifdef sparc
#define YES_ALLOCA
#include <alloca.h>
#endif

#include "defs-env.h"

// matrix/mexdep.c
extern char *alloc_check(cint, cint);

#define Void void

#define	ALLOC_CHECK(n,s)	alloc_check(n,s)

#ifdef YES_ALLOCA
	extern char *alloca();
#	define TMP_Alloc(n,s)	alloca((int) ((n)*(s)))
#	define TMP_Free(p)	{ }
#else
#	define TMP_Alloc(n,s)	Calloc(n,s)
#	define TMP_Free(p)	Free(p)
#endif

#define mycalloc(i, j)		calloc(i,j)
#define myfree(p)		Free(p)
#define message(str, val)	Note2("message\n%s %g", str, val)
#define message0		message
#define jumpship(msg)		Exit(msg)

// matrix/debug.c
extern void print_vector_d(cdouble *, cint);
extern sof scan_vector_d(double *, cint);
extern void print_vector_f(cfloat *, cint);
extern sof scan_vector_f(float *, cint);

// matrix/io_d.c
extern void print_matrix_d(cdouble **mat, cint N, cint M);
extern sof write_vector_d(cchar *name, cdouble *vect, cint N);
extern sof fprint_matrix_d(cchar *name, cint *ext, cdouble **mat, cint N, cint M);
extern sof fscan_vector_d(double *vect, cchar *name, cint N);
extern sof read_vector_d(double *vect, cchar *name, cint N);



// USE "ABORT" FOR MEX-linked ROUTINES, "EXIT" FOR BATCH C
#define ABORT(msg)	{ char str[512]; \
			sprintf(str, "ERROR: %s, in %s line %d\n", msg, __FILE__, __LINE__); \
			jumpship(str);}

#define MAX_LINE_LENGTH 256 // MAX LENGTH FOR LINES OF ASCII FILES


// ?
extern	double brent(float ax, float bx, float cx, double (*f)(float),
	float tol, float *xmin);
extern	void linmin(float *pnt, float *dir, int MM, float *fret, double (*func)(float *));
extern	void mnbrak(float *ax, float *bx, float *cx,
	float *fa, float *fb, float *fc, double (*func)(float));
extern	void powell(float *p, float **xi, int MM,
	float ftol, int *iter, float *fret, double (*func)(float *));

#endif // DEFS_VSPLINE
