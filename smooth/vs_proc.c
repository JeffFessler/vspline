/****************************** vs_proc.c *********************************
*	SUBROUTINES FOR CUBIC SPLINE SMOOTHING ALGORITHM
*	THESE ALL DEPEND ON alpha - THE SMOOTHING PARAMETERS.
*	SEE vs_smooth.c and vs_init.c FOR THE DIMENSIONS OF Qt etc.
*	Jeff Fessler
*	COPYRIGHT 1989,  THE BOARD OF TRUSTEES OF THE
*	LELAND STANFORD JUNIOR UNIVERSITY.  ALL RIGHTS RESERVED.
**************************************************************************/
#include "smooth.h"


/****************************** form_band() ******************************
*	FORM THE BAND MATRIX USED IN THE REINSCH ALGORITHM
*	B = (R @ inv(alpha)) + (Qt @ I) Sigma (Q @ I)
*	INPUT:
*		QtSQ:	SYMBAND(3*M-1, M*(N-2))
*		R_A:	SYMBAND(2*M-1, M*(N-2))
*	OUTPUT:
*		sum:	SYMBAND(3*M-1, M*(N-2))		R_A + QtQS
*	FLOPS: 2MM(N-2)
*************************************************************************/
void form_band(SYMBAND sum, cSYMBAND R_A, cSYMBAND QtSQ, cint M, cint N)
{
	register int n, m;

	for (m=3*M; m > M; m--)
		for (n = M*(N-2); n >= 1; n--)
			sum[n][m] = R_A[n][m-M] + QtSQ[n][m];

	// TRICKY: THIS ALSO FILLS IN DESIRABLE 0'S IN UPPER LEFT OF sum
	for (n = M*(N-2); n >= 1; n--)
		COPY_VECTOR(QtSQ[n], sum[n], M);
}


/****************************** form_q_c() *******************************
*	COMPUTE (Q @ I) * c_hat
*	INPUT:
*		Qt:		see above
*		c[1..M(N-2)]	estimated coefficients
*	OUTPUT:
*		q_c[1..MN]	(Q @ I) * c_hat
*	FLOPS: 6M(N-2)
*************************************************************************/
void form_q_c(VECT q_c, cBAND Qt, cVECT c, cint M, cint N)
{
	ZERO_VECTOR(q_c, M*N);

	INC_SCALED_VECTOR(q_c, (double) Qt[1][3], c, M);
	q_c += M;

	INC_SCALED_VECTOR(q_c, (double) Qt[1][2], c,     M);
	INC_SCALED_VECTOR(q_c, (double) Qt[2][3], c + M, M);
	q_c += M;

	for (int n=2; n <= N-3; n++, q_c += M)
	{
		INC_SCALED_VECTOR(q_c, (double) Qt[n-1][1], c + M*(n-2), M);
		INC_SCALED_VECTOR(q_c, (double) Qt[n  ][2], c + M*(n-1), M);
		INC_SCALED_VECTOR(q_c, (double) Qt[n+1][3], c + M*n,     M);
	}

	INC_SCALED_VECTOR(q_c, (double) Qt[N-3][1], c + M*(N-4), M);
	INC_SCALED_VECTOR(q_c, (double) Qt[N-2][2], c + M*(N-3), M);
	q_c += M;

	INC_SCALED_VECTOR(q_c, (double) Qt[N-2][1], c + M*(N-3), M);
}


/****************************** form_sig_q_c() ***************************
*	COMPUTE diag(sigma) * (Q' @ I) * c_hat
*	INPUT:
*		covar[1..?]		error covariance structure(s)
*		q_c[1..MN]		(Q' @ I) * c_hat
*		flag_iid		is the noise iid?
*	OUTPUT:
*		sig_q_c[1..N][1..M]
*	FLOPS: NM(M+M-1) = N(2MM-M)
*************************************************************************/
void form_sig_q_c(OMATRIX sig_q_c, cCOVAR *covar, cCMAT q_c,
cint M, cint N,
ctruf flag_iid)
{
	for (int n=1; n <= N; n++, q_c += M)
	{
		OMATRIX cov = covar[flag_iid ? 1 : n].cov;
		mult_full_vect(sig_q_c[n], cov, q_c, M);
	}
}
