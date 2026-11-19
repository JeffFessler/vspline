/***************************** smooth.h **********************************
* Declarations for vector-spline smoothing algorithm.
* Jeff Fessler
* COPYRIGHT 1989,  THE BOARD OF TRUSTEES OF THE
* LELAND STANFORD JUNIOR UNIVERSITY.  ALL RIGHTS RESERVED.
*************************************************************************/
#include "defs_vmatrix.h"

/********************** BAND MATRIX FORMAT *******************************
*	IF B IS A nxN SYMMETRIC BAND MATRIX WITH 2*J+1 NON-ZERO DIAGONALS
*	THEN I AM USING THE "TRANSPOSE" OF THE LINPACK BAND MATRIX FORMAT.
*	(IE, THE GOLUB - VAN LOAN FORMAT).
*	SUCH A MATRIX IS DENOTED BAND(J,N) IN COMMENTS.
*
*	B(i,j) = B[i][j-i+J+1],		i >= j
*
*	B LOOKS LIKE:
*              |   0   0                 b(1,1)     |
*              |   0              b(2,1) b(2,2)     |
*              |               .           :        |
*              |        .                  :        |
*              |  b(J+1,1)     ...       b(J+1,J+1) |
*              |     :                     :        |
*              |  b(N,N-J)     ...       b(N,N)     |
*
*	SO THE COLUMNS OF B[][] ARE THE DIAGONALS OF B()
*	AND THE ROWS OF B[][] ARE THE (LOWER BAND) ROWS OF B().
*********************** BAND MATRIX FORMAT ******************************/

#define OPT_ARG	float			/* Args to optimizing routines	*/
#define OPT_FUNC double			/* functions to optimize	*/

#ifdef DoublePrec
//	typedef double DATA;
#define DATA double

#define ALLOC_MATRIX			ALLOC_MATRIX_D
#define ALLOC_VECTOR			ALLOC_VECTOR_D
#define TMP_ALLOC_MATRIX		TMP_ALLOC_MATRIX_D

#define SCALE_VECTOR			SCALE_VECTOR_D
#define SCALE_MATRIX			SCALE_MATRIX_D
#define INC_SCALED_VECTOR		INC_SCALED_VECTOR_D
#define DOT_VECTOR			DOT_VECTOR_D
#define NORM_VECTOR			NORM_VECTOR_D
#define SUB_VECTOR			SUB_VECTOR_D
#define DEC_VECTOR			DEC_VECTOR_D
#define INC_VECTOR			INC_VECTOR_D

#define PRINT_VECTOR			PRINT_VECTOR_D
#define PRINT_MATRIX			PRINT_MATRIX_D
#define WRITE_VECTOR			WRITE_VECTOR_D
#define FPRINT_MATRIX			FPRINT_MATRIX_D
#define SCAN_VECTOR			SCAN_VECTOR_D
#define SCAN_MATRIX			SCAN_MATRIX_D
#define READ_VECTOR			READ_VECTOR_D
#define FSCAN_VECTOR			FSCAN_VECTOR_D

#else
	NOT IMPLEMENTED SINCE DOUBLE PRECISION SEEMS BEST!
#endif

#define ALLOC_BAND(n,j)			ALLOC_MATRIX((n),(j+1))
#define FREE_BAND			FREE_MATRIX
#define COPY_BAND(src,dst,n,j)		COPY_MATRIX(src,dst,n,(j+1))
#define ZERO_BAND(mat,n,j)		ZERO_MATRIX(mat,(n),(j+1))
#define PRINT_BAND(m,n,j)		PRINT_MATRIX(m,(n),(j+1))
#define TMP_ALLOC_VECT			TMP_ALLOC_VECTOR_D
#define TMP_ALLOC_BAND(m,n,j)		TMP_ALLOC_MATRIX(m,n,(j+1))
#define TMP_FREE_BAND			TMP_FREE_MATRIX

#define JOBUR	1		// CONTROL BITS FOR cross_val()
#define JOBCV	2
#define JOBGCV	4

typedef DATA **OMATRIX;		// ORDINARY MATRIX [1..rows][1..cols]
#define SYMBAND OMATRIX		// SYMMETRIC BAND MATRIX
#define cSYMBAND Const SYMBAND
#define BAND SYMBAND		// BAND MATRIX
#define cBAND Const BAND
typedef DATA *VECT;		// VECTOR[1..N]
#define cVECT Const DATA *
#define CMAT VECT		// MATRIX STORED COLUMN-WISE AS A VECT
#define cCMAT Const DATA *

#define ALPHA		OPT_ARG
typedef ALPHA		*ALPHA_SET;
#define ALLOC_ALPHA	ALLOC_VECTOR_F
#define PRINT_ALPHA	PRINT_VECTOR_F
#define SCAN_ALPHA	SCAN_VECTOR_F
#define FSCAN_ALPHA	FSCAN_VECTOR_F
#define SETUP_DIRSET	identity_f
typedef OPT_ARG **DIRSET;

/*************************************************************************
*	Structure for a MxM covariance matrix and its inverse.
*	(SYMMETRIC, POSITIVE DEFINITE)
*	inv(cov) = inv = Linv Dinv Linv'
*	For nonlinear problems, inv is typically the Hessian.
*	cov2 = cov * cov (used for unbiased risk score).
*************************************************************************/

typedef struct
{
	SYMBAND cov;				/* SYMBAND(M-1,M)	*/
	SYMBAND cov2;				/* SYMBAND(M-1,M)	*/
	SYMBAND inv;				/* SYMBAND(M-1,M)	*/
	SYMBAND Linv;				/* SYMBAND(M-1,M)	*/
	VECT Dinv;				/* [1..M]		*/
	DATA trace;				/* trace(cov)		*/
} COVAR;

typedef Const COVAR cCOVAR;

#include "decs-smooth.h"
