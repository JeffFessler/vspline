/*********************** example1.c **************************************
* EXAMPLE FRONT END ROUTINE FOR USING VECTOR-SPLINE SUBROUTINES TO
* SMOOTH VECTOR MEASUREMENTS, COMPUTE CV SCORES, ETC.
* Jeff Fessler
* COPYRIGHT 1989, THE BOARD OF TRUSTEES OF THE
* LELAND STANFORD JUNIOR UNIVERSITY. ALL RIGHTS RESERVED.
*************************************************************************/
#include "smooth.h"


/*
*	VECTOR SPLINE SMOOTH AND COMPUTE SCORES
*	INPUT:
*		alpha[1..M]		smoothing parameters
*		ti[1..N]		sample points
*		yi[1..nxM]		measurements
*		isigma[1..MxMx(N)]	inverse of measurement error covariances
*		flag_iid		0: one covariance, 1: N covariances
*	OUTPUT:
*		ys[1..nxM]		smoothed estimates
*	RETURNS:			1 if successful, 0 if error
*/
sof vspline(CMAT ys, ALPHA_SET alpha, VECT ti, CMAT yi, CMAT isigma,
ctruf flag_iid, cint MM, cint nn)
{
	COVAR *covar;

	// convert the flattened covariance inverses to structure
	Call(covar = inv_to_covar_set, (isigma, MM, nn, flag_iid))
	Call(vs_setup_alloc, (MM, nn))
	Call(vs_setup_ti, (ti, nn))
	Call(vs_setup_yi, (yi, MM, nn))
	Call(vs_setup_cov, (covar, flag_iid))
	Call(vs_setup_alpha, (alpha))
	Call(vs_smooth, (ys, 1))
	vs_cleanup();
	covar_free(covar, flag_iid ? 1 : nn);

	Ok
}


/*
* VECTOR SPLINE SMOOTH AND COMPUTE SCORES
* INPUT:
*	alpha[1..M]		smoothing parameters
*	ti[1..N]		sample points
*	yi[1..nxM]		measurements
*	covar[1..?]		measurement error covariances
*	flag_iid		covariance varies?
* OUTPUT:
*	ys[1..nxM]		smoothed estimates
*	ur			unbiased risk score
*	cv			cross validation score
*	gcv			generalized CV score
*	penin			initial roughness penalty
*	penout			final roughness penalty
*	rough			final roughness penalty
*/
sof vect_smooth(CMAT ys, double *ur, double *cv, double *gcv,
double *penin, double *penout, double *rough, ALPHA_SET alpha,
VECT ti, CMAT yi, COVAR *covar, ctruf flag_iid, cint MM, cint nn)
{
	// SETUP: INDEPENDENT OF ALPHA
	Call(vs_setup_alloc, (MM, nn))
	Call(vs_setup_ti, (ti, nn))
	Call(vs_setup_yi, (yi, MM, nn))
	Call(vs_setup_cov, (covar, flag_iid))
	Call(cv_setup_alloc, ())
	Call(cv_setup_covar, ())

	// THE ROUTINES BELOW ALL DEPEND ON ALPHA
	Call(vs_setup_alpha, (alpha))
	Call(vs_smooth, (ys, 1))
	Call(cross_val, (ur, cv, gcv, JOBUR | JOBCV | JOBGCV))

	*rough = roughness();
	Call(penalty, (penout, alpha, ys, MM, nn))
	Call(penalty, (penin, alpha, yi, MM, nn))

	vs_cleanup();
	cv_cleanup();
	Ok
}
