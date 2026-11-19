/****************************** vs_init.c *********************************
* SUBROUTINES FOR CUBIC SPLINE SMOOTHING ALGORITHM;
* THESE ARE ALL INDEPENT OF THE SMOOTHING PARAMETERS.
* SEE band.h FOR THE FORMAT OF THE BAND MATRICES.
* SEE vs_smooth.c FOR DIMENSIONS OF
*		R, QtQ, QtSQ
*
* EXCEPTION: SINCE Qt IS UPPER BAND AND NOT SYMMETRIC,
* THE FOLLOWING STORAGE FORMAT IS USED:
* Qt(i,i+k) = Qt[i][3-k], for k=0,1,2  i = 1, ..., N-2
* Jeff Fessler
* COPYRIGHT 1989, THE BOARD OF TRUSTEES OF THE
* LELAND STANFORD JUNIOR UNIVERSITY. ALL RIGHTS RESERVED.
**************************************************************************/
#include "smooth.h"


/****************************** form_qt_s_q_reg() ************************
*	COMPUTE QtSQ = (Q' @ I) diag(sigma) (Q @ I)
*	INPUT:
*		Qt:		see above
*		covar[1..N]	error covariance structures
*	OUTPUT:
*		QtSQ
*	FLOPS:	(N-4)(MM+1)+(N-3)(3MM+2)+(N-2)(5MM+3) = N(9MM+6)-23MM-16
*		OLD (N-2) M M 2 + (N-1) M M 5 + (N-1) M M 8 = 15MMN-17MM
*************************************************************************/
static void form_qt_s_q_reg(SYMBAND QtSQ, cBAND Qt, cCOVAR *cov, cint M, cint N)
{
	register int m, k, n;

	for (n=N-2; n >= 3; n--)
	{
		register DATA t1 = Qt[n][3] * Qt[n-2][1];

	  for (m=M; m >= 1; m--)
	  {
	    for (k=M; k >= m; k--)
		QtSQ[(n-1)*M+m][k-m+M] =	t1 * cov[n  ].cov[k][m-k+M];
	    for (; k >= 1; k--)
		QtSQ[(n-1)*M+m][k-m+M] =	t1 * cov[n  ].cov[m][k-m+M];
	  }
	}

	for (n=N-2; n >= 2; n--)
	{
		register DATA t1 = Qt[n][3] * Qt[n-1][2];
		register DATA t2 = Qt[n][2] * Qt[n-1][1];

	  for (m=M; m >= 1; m--)
	  {
	    for (k=M; k >= m; k--)
		QtSQ[(n-1)*M+m][k-m+2*M] =	t1 * cov[n  ].cov[k][m-k+M] +
						t2 * cov[n+1].cov[k][m-k+M];
	    for (; k >= 1; k--)
		QtSQ[(n-1)*M+m][k-m+2*M] =	t1 * cov[n  ].cov[m][k-m+M] +
						t2 * cov[n+1].cov[m][k-m+M];
	  }
	}

	for (n=N-2; n >= 1; n--)
	{
		register DATA t1 = Qt[n][3] * Qt[n][3];
		register DATA t2 = Qt[n][2] * Qt[n][2];
		register DATA t3 = Qt[n][1] * Qt[n][1];

	  for (m=M; m >= 1; m--)
	  {
	    for (k=m; k >= 1; k--)
		QtSQ[(n-1)*M+m][k-m+3*M] =	t1 * cov[n  ].cov[m][k-m+M] +
						t2 * cov[n+1].cov[m][k-m+M] +
						t3 * cov[n+2].cov[m][k-m+M];
	  }
	}
}


/****************************** form_qt_s_q_iid() ************************
*	COMPUTE QtSQ = (Q' @ I) diag(sigma) (Q @ I)
*	INPUT:
*		Qt:		see above
*		covar		one covariance matrix structure
*	OUTPUT:
*		QtSQ
*	FLOPS:
*************************************************************************/
static void form_qt_s_q_iid(SYMBAND QtSQ, cBAND Qt, cCOVAR *covar, cint M, cint N)
{
	register int m, k, n;
	SYMBAND cov = covar->cov;

	for (n=N-2; n >= 3; n--)
	{
		register DATA t1 = Qt[n][3] * Qt[n-2][1];

	  for (m=M; m >= 1; m--)
	  {
	    for (k=M; k >= m; k--)
		QtSQ[(n-1)*M+m][k-m+M] =	t1 * cov[k][m-k+M];
	    for (; k >= 1; k--)
		QtSQ[(n-1)*M+m][k-m+M] =	t1 * cov[m][k-m+M];
	  }
	}

	for (n=N-2; n >= 2; n--)
	{
		register DATA t1 =	Qt[n][3] * Qt[n-1][2] +
					Qt[n][2] * Qt[n-1][1];

	  for (m=M; m >= 1; m--)
	  {
	    for (k=M; k >= m; k--)
		QtSQ[(n-1)*M+m][k-m+2*M] =	t1 * cov[k][m-k+M];
	    for (; k >= 1; k--)
		QtSQ[(n-1)*M+m][k-m+2*M] =	t1 * cov[m][k-m+M];
	  }
	}

	for (n=N-2; n >= 1; n--)
	{
		register DATA t1 =	Qt[n][3] * Qt[n][3] +
					Qt[n][2] * Qt[n][2] +
					Qt[n][1] * Qt[n][1];

	  for (m=M; m >= 1; m--)
	  {
	    for (k=m; k >= 1; k--)
		QtSQ[(n-1)*M+m][k-m+3*M] =	t1 * cov[m][k-m+M];
	  }
	}
}


/****************************** spline_mat() *****************************
*	COMPUTE R, QtQ = Q'*Q, AND Qt =  Q' GIVEN KNOTS ti
*	IF ti IS NULL, THEN THE SAMPLES ARE SPACED EQUALLY.
*	INPUT:
*		ti[1..N]	sample points (increasing order) (or null)
*		N: # of points
*	OUTPUT:
*		R, QtQ, Qt
*	RETURNS:	0 if ti's not in increasing order, 1 otherwise
*	FLOPS: N-1 + 6(N-1) = 7N - 7
*		PLUS (N-1) * 9 for QtQ - NOT NEEDED FOR NON-IID
*************************************************************************/
sof spline_mat(SYMBAND R, BAND Qt, cVECT ti, cint N)
{
	register int i, i1;
	VECT hi;
	TMP_ALLOC_VECT(hi, N-1)

	// COMPUTE DIFFERENCES OF THE KNOTS
	if (ti)
	{
		for (i=N-1; i >= 1; i--)
			if ( (hi[i] = ti[i+1]-ti[i]) <= 0. )
				return Failure;
	}
	else
		for (i=N-1; i >= 1; i--)
			hi[i] = 1.;

	for (i=N-2, i1=N-1; i >= 1; i--, i1--)
	{
		R[i][1] =  hi[i] / 6.;
		R[i][2] = (hi[i] + hi[i1]) / 3.;

		Qt[i][1] = 1. / hi[i1];
		Qt[i][3] = 1. / hi[i];
		Qt[i][2] = -Qt[i][3] - Qt[i][1];
	}

	R[1][1] = 0.; // CORRECTION

#ifdef NOTUSED
	i = N-2;
	QtQ[i ][3] = Qt[i][3]*Qt[i ][3] + Qt[i][2]*Qt[i ][2]+ Qt[i][1]*Qt[i][1];
	i1 = i--;
	QtQ[i1][2] = Qt[i][2]*Qt[i1][3] + Qt[i][1]*Qt[i1][2];
	QtQ[i ][3] = Qt[i][3]*Qt[i ][3] + Qt[i][2]*Qt[i ][2]+ Qt[i][1]*Qt[i][1];
	i2 = i1;
	i1 = i--;
    for ( ; i >= 1; i--, i1--, i2--)
    {
	QtQ[i2][1] = Qt[i][1]*Qt[i2][3];
	QtQ[i1][2] = Qt[i][2]*Qt[i1][3] + Qt[i][1]*Qt[i1][2];
	QtQ[i ][3] = Qt[i][3]*Qt[i ][3] + Qt[i][2]*Qt[i ][2]+ Qt[i][1]*Qt[i][1];
    }
#endif
	TMP_FREE_VECTOR(hi);
	Ok
}


/****************************** form_qt_Y() ******************************
*	COMPUTE (Q' @ I) Y
*	INPUT:
*		Qt:		see above
*		Y[1..MN]	measured data (Y1_1 .. YM_1, .., Y1_1 .. YM_N)
*	OUTPUT:
*		QtY[1..M(N-2)]	(Qt @ I) * Y
*	FLOPS: M(5N-5)
*************************************************************************/
void form_qt_Y(VECT QtY, cBAND Qt, cVECT Y, cint M, cint N)
{
	register int m, i;

	for (m=M; m >= 1; m--)
	    for (i=N-2; i >= 1; i--)
		QtY[(i-1)*M+m] =	Y[(i-1)*M+m] * Qt[i][3] +
					Y[(i  )*M+m] * Qt[i][2] +
					Y[(i+1)*M+m] * Qt[i][1];
}


/****************************** form_qt_s_q() ****************************
*	COMPUTE (Qt @ I) Sigma (Q @ I)
*	INPUT:
*		Qt			see above
*		covar[1..?]		error covariance structure(s)
*		flag_iid:		more than one covariance?
*	OUTPUT:
*		QtSQ
*************************************************************************/
void form_qt_s_q(SYMBAND QtSQ, cBAND Qt, cCOVAR *covar, cint M, cint N,
ctruf flag_iid)
{
	if (flag_iid)
		form_qt_s_q_iid(QtSQ, Qt, &covar[1], M, N);
	else
		form_qt_s_q_reg(QtSQ, Qt, covar, M, N);
}
