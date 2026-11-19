// mexdep.h
// Header files for code that uses MEX dependent subroutines.
// For compiling Matlab MEX files, must use the MEX library rather than stdio.
// Quite a pain, frankly.
// Affected functions: malloc(), calloc(), free(), printf()
// Jeff Fessler
// COPYRIGHT 1989, THE BOARD OF TRUSTEES OF THE
// LELAND STANFORD JUNIOR UNIVERSITY. ALL RIGHTS RESERVED.

#include "defs-vspline.h"

#ifdef MATLAB
#include "cmex.h"
#endif


#ifdef MATLAB

#ifdef DoublePrec

/* This is the right way...
#define	MEX_TO_VECTOR(v, mat)	Bcopy(mat->pr, &v[1], mat->n*mat->m)
#define	VECTOR_TO_MEX(mat, v)	Bcopy(&v[1], mat->pr, mat->n*mat->m)
*/

#define	MEX_TO_VECTOR(v, mat)	{v = (mat->pr - 1);}
#define	VECTOR_TO_MEX(mat, v)

#else

#define	MEX_TO_VECTOR(v, mat)	vector_d_to_f(v, mat->pr-1, mat->n*mat->m)
#define	VECTOR_TO_MEX(mat, v)	vector_f_to_d(mat->pr-1, v, mat->n*mat->m)

#endif // DoublePrec

#define	MEX_TO_FLOAT(v, mat)	vector_d_to_f(v, mat->pr-1, (int) (mat->n*mat->m))
#define	FLOAT_TO_MEX(mat, v)	vector_f_to_d(mat->pr-1, v, (int) (mat->n*mat->m))

#endif // MATLAB
