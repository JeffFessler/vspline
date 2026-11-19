/*********************** example2.c **************************************
* EXAMPLE FRONT END ROUTINE FOR USING VECTOR-SPLINE SUBROUTINES TO
* SMOOTH VECTOR MEASUREMENTS, COMPUTE CV SCORES, ETC.
* THIS EXAMPLE CHOOSES THE SMOOTHING PARAMETER AUTOMATICALLY,
* Jeff Fessler
* COPYRIGHT 1989, THE BOARD OF TRUSTEES OF THE
* LELAND STANFORD JUNIOR UNIVERSITY.  ALL RIGHTS RESERVED.
*************************************************************************/
#include "smooth.h"


/*
* USE CV SCORE TO DETERMINE SMOOTHING PARAMETER
*	alpha:		initial value
*	all others as in vspline() above
* OUTPUT:
*	alpha:		optimal alpha
*	ys:		best estimates
* RETURNS:		1 if ok, 0 if not
*/
#define TOL	0.01

static sof bestcv(
CMAT ys,
ALPHA_SET alpha,
VECT ti,
CMAT yi, CMAT isigma,
ctruf flag_iid,
cint MM, cint nn)
{
	COVAR *covar;
	Call(covar = inv_to_covar_set, (isigma, MM, nn, flag_iid))

	best_score(ys, alpha, ti, yi, covar, flag_iid, MM, nn, JOBCV, 0.01, True);
	covar_free(covar, flag_iid ? 1 : nn);
	Ok
}
