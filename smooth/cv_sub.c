// cv_sub.c
// SUBROUTINES FOR COMPUTING CROSS-VALIDATION LIKE SCORES.
// Jeff Fessler
// COPYRIGHT 1989, THE BOARD OF TRUSTEES OF THE
// LELAND STANFORD JUNIOR UNIVERSITY. ALL RIGHTS RESERVED.

#include "smooth.h"


/***************************** select_nth_qbiqt() ************************
*	COMPUTE THE n'TH CENTRAL M*M INFLUENCE MATRIX OF
*	Fnn = [(Q @ I_M) IB (Qt @ I_M)]nn
*	INPUT:
*		IB: SYMBAND(3M-1, M(N-2))	central inverse
*		Qt:				see vs_init.c
*		n:				WHICH ONE TO DO (1,..., N).
*	OUTPUT:
*		F: SYMBAND(M-1, M)		FULL SYMMETRIC M*M MATRIX
*	FLOPS: M*(M+1)/2* 3 * (3*3)	(1*1) for n=1, (2*2) for n=2
*	all N calls:	N(27/2M(M+1)) - 39/2M(M+1)
*************************************************************************/
void select_nth_qbiqt(
SYMBAND F,
cSYMBAND IB, cSYMBAND Qt, cint M, cint N, cint n)
{
	int i, j, m, k, n1, nN, i1, i2;
	double *sum;
	int J1 = 3*M-1 + 1;

	n1 = (n-2 <= 1) ? 1 : n-2;
	nN = (n <= N-2) ? n : N-2;

	for (i=1; i <= M; i++) // LOWER TRIANGLE OF F
	{
	for (j=1; j <= i; j++)
	{
		sum = &F[0*M+i][j-i+M]; // 0 instead of n-1
		*sum = 0.;

		for (k=n1; k <= nN; k++)
			for (m=n1; m <= nN; m++)
			{
				i1 = (k-1)*M+i;
				i2 = (m-1)*M+j;

				*sum += Qt[m][3-(n-m)] * Qt[k][3-(n-k)] *
				((i2 < i1) ? IB[i1][i2-i1+J1] : IB[i2][i1-i2+J1]);
			}
	}
	}
}


/********************* rss_term() ****************************************
*	COMPUTE A SINGLE TERM IN THE RESIDUAL SUM OF SQUARES
*	rss_n = err_n' * inv(cov_n) * err_n
*	INPUT:
*		err[1..M]		nth error term = yi^n - ys^n
*		inv: SYMBAND(M-1,M)	inverse of covariance_n
*	RETURNS:
*		rss_n
*	FLOPS:
*************************************************************************/
double rss_term(cVECT err, cSYMBAND inv, cint M)
{
	register double rss;
	VECT werr;
	TMP_ALLOC_VECT(werr, M)

	mult_full_vect(werr, inv, err, M);
	rss = DOT_VECTOR(err, werr, M);

	TMP_FREE_VECTOR(werr);
	return rss;
}


/********************* ur_term() *****************************************
*	COMPUTE A SINGLE TERM IN THE UNBAISED RISK SUMMATION.
*	UR = (1/N) [ || err ||^2 - 2 * trace(cov * (I - A)) + trace(cov) ]
*	ur_n = trace(cov_n * cov_n * Fn)
*	INPUT:
*		Fn: SYMBAND(M-1,M)	nth qbiqt term
*		cov2: SYMBAND(M-1,M)	covariance_n^2
*	RETURNS:
*		ur_n:			(part of) nth term of ur score
*	FLOPS:
*************************************************************************/
double ur_term(cSYMBAND Fn, cSYMBAND cov2, cint M)
{
	return trace_prod(cov2, Fn, M);
}


/********************* cv_term() *****************************************
*	COMPUTE A SINGLE TERM IN THE CROSS VALIDATION SUMMATION.
*	CV = (1/N) sum_1^N || sign^(-1/2) inv(I - An) err_n ||^2
*	cv = werr' * inv(cov) * werr, Fn werr = inv(cov) * err
*	INPUT:
*		err[1..M]		nth error term = yi^n - ys^n
*		FL: SYMBAND(M-1,M)
*		FD[1..M]		Fn = FL' FD FL
*		inv: SYMBAND(M-1,M)	inverse of covariance_n
*	RETURNS:
*		nth term of cv score
*	FLOPS:
*************************************************************************/
double cv_term(cVECT err, cSYMBAND FL, cVECT FD, cSYMBAND inv, cint M)
{
	double sum;
	VECT werr, invwerr;

	TMP_ALLOC_VECT(werr, M)
	TMP_ALLOC_VECT(invwerr, M)

	mult_full_vect(werr, inv, err, M);
	solve_sym_band((VECT) 0, FL, FD, werr, M-1, M, 1);
	mult_full_vect(invwerr, inv, werr, M);

	sum = DOT_VECTOR(invwerr, werr, M);

	TMP_FREE_VECTOR(werr);	TMP_FREE_VECTOR(invwerr);
	return sum;
}


/********************* gcv_term() ****************************************
*	COMPUTE A SINGLE TERM IN THE GENERALIZED CROSS VALIDATION SUMMATION.
*	GCV = (1/N) RSS / [ (1/N) trace(I - A) ]^2
*	den_gcv_n = trace(I - An) = trace(cov_n * Fn)
*	INPUT:
*		Fn: SYMBAND(M-1,M)
*		cov: SYMBAND(M-1,M)	covariance_n
*	RETURNS:
*		den_gcv_n:		nth term of denominator of gcv score
*	FLOPS:
*************************************************************************/
double gcv_term(cSYMBAND Fn, cSYMBAND cov, cint M)
{
	return trace_prod(cov, Fn, M);
}
