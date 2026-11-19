/****************************** band_util.c *******************************
* UTILITY SUBROUTINES FOR BAND MATRICES
* Jeff Fessler
* COPYRIGHT 1989, THE BOARD OF TRUSTEES OF THE
* LELAND STANFORD JUNIOR UNIVERSITY. ALL RIGHTS RESERVED.
*************************************************************************/
#include "smooth.h"


/****************************** band_trace() ******************************
*	TRACE OF A SYMMETRIC BAND MATRIX
*	INPUT:
*		A: SYMBAND(J,M)
*	RETURNS:
*		trace(A)
*************************************************************************/
double band_trace(cSYMBAND A, int J, int M)
{
	register double trace = 0.;
	J++;
	for (; M; M--)
		trace += A[M][J];
	return trace;
}


/****************************** trace_prod() ******************************
*	TRACE OF PRODUCT OF TWO SYMMETRIC BAND MATRICES
*	INPUT:
*		A: SYMBAND(M-1,M)
*		B: SYMBAND(M-1,M)
*	RETURNS:
*		trace( A * B )
*	FLOPS:
**************************************************************************/
double trace_prod(cSYMBAND A, cSYMBAND B, cint M)
{
	register double trace = 0.;
	register int m, i;

	for (m=M; m >= 1; m--)
	{
		for (i=M; i >= m; i--)
			trace += A[i][m-i+M] * B[i][m-i+M];
		for (   ; i >= 1; i--)
			trace += A[m][i-m+M] * B[m][i-m+M];
	}

	return trace;
}


/****************************** mult_band_vect() **************************
*	MULTIPLY A BAND MATRIX TIMES A VECTOR
*	INPUT:
*		A: SYMBAND(J,N)
*		x[1..N]
*	OUTPUT:
*		y[1..N]		y = A x
*	FLOPS:	N(2J+1)2J - ?
**************************************************************************/
void mult_band_vect(VECT y, cSYMBAND A, cVECT x, cint J, cint N)
{
	register int n, k;
	/* ASSUMES J+J <= N	*/
	if (J+J > N)
		ABORT("Band matrix too fat for this algorithm")

	for (n=J; n >= 1; n--)
	{
		y[n] = A[n][J+1] * x[n];
		for (k=J; k >= 1; k--)
			y[n] += A[n+k][-k+J+1] * x[n+k]; // A(n,n+k) A(n+k,n)
		for (k=n-1; k >= 1; k--)
			y[n] += A[n][-k+ J+1] * x[n-k]; // A(n, n-k)
	}

	for (n=J+1; n <= N; n++)
	{
		int k0 = (n+J > N) ? N-n : J;
		y[n] = A[n][J+1] * x[n];
		for (k=k0; k >= 1; k--)
			y[n] += A[n+k][-k+J+1] * x[n+k]; // A(n,n+k) A(n+k,n)
		for (k=J; k >= 1; k--)
			y[n] += A[n][-k+ J+1] * x[n-k]; // A(n, n-k)
	}
}


/****************************** mult_full_vect() **************************
*	MULTIPLY A FULL BAND MATRIX TIMES A VECTOR
*	INPUT:
*		A: SYMBAND(M-1,M)
*		z[1..M]
*	OUTPUT:
*		v[1..M]		v = A z
*	FLOPS:	M(2M-1)
**************************************************************************/
void mult_full_vect(VECT v, cSYMBAND A, cVECT z, cint M)
{
	register int m, k;

	for (m=M; m >= 1; m--)
	{
		v[m] = A[m][M] * z[m];
		for (k=m-1; k >= 1; k--)
			v[m] += A[m][k-m + M] * z[k];
		for (k=m+1; k <= M; k++)
			v[m] += A[k][m-k + M] * z[k];
	}
}


/****************************** mult_full_band() **************************
*	MULTIPLY TWO FULL M*M MATRICES STORED IN BAND FORMAT
*	INPUT:
*		A: SYMBAND(M-1, M)
*		B: SYMBAND(M-1, M)
*	OUTPUT:
*		C: SYMBAND(M-1, M)	C = A * B
*************************************************************************/
void mult_full_band(SYMBAND C, cSYMBAND A, cSYMBAND B, cint M)
{
	register int i, j, k;
	register DATA *sum;

	for (i=M; i >= 1; i--)
	{
		for (j=i; j >= 1; j--)
		{
			*(sum = &C[i][j-i + M]) = 0.;
			for (k=1; k <= j; k++)
				*sum += A[i][k-i + M] * B[j][k-j + M];
			for (k=j+1; k <= i; k++)
				*sum += A[i][k-i + M] * B[k][j-k + M];
			for (k=i+1; k <= M; k++)
				*sum += A[k][i-k + M] * B[k][j-k + M];
		}
	}
}

/*************************************************************************
*	SUBROUTINES FOR CONVERTING BETWEEN BAND AND COLUMN-WISE FORMATS.
*************************************************************************/

/****************************** sym_band_to_mat() ************************
*	CONVERT A SYMMETRIC BAND MATRIX B TO CMAT FORMAT MATRIX A
*	INPUT:
*		B: SYMBAND(J, N)
*	OUTPUT:
*		A: CMAT[1..nxN]
*************************************************************************/
void sym_band_to_mat(CMAT A, cSYMBAND B, cint J, cint N)
{
	register int i, j;

	for (i=1; i <= J; i++)
		for (j=1; j <= i; j++)
			A[(i-1)*N + j] = A[(j-1)*N + i] = B[i][j-i+J+1];

	for (i=J+1; i <= N; i++)
		for (j=i-J; j <= i; j++)
			A[(i-1)*N + j] = A[(j-1)*N + i] = B[i][j-i+J+1];
}


/****************************** mat_to_sym_band() ************************
*	CONVERT A SYMMETRIC COLUMN-WISE MATRIX A TO BAND FORMAT MATRIX B
*	INPUT:
*		A: CMAT[1..nxN]
*	OUTPUT:
*		B: SYMBAND(J, N)
*************************************************************************/
void mat_to_sym_band(SYMBAND B, cCMAT A, cint J, cint N)
{
	for (int j=J+1; j >= 1; j--)
		for (int i=J+2-j; i <= N; i++)
			B[i][j] = A[(i-1)*N + i+j-(J+1)];
}
