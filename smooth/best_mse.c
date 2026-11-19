/*********************** best_mse.c **************************************
* VECTOR SPLINE SMOOTHING WITH SMALLEST MSE
* Jeff Fessler
* COPYRIGHT 1989, THE BOARD OF TRUSTEES OF THE
* LELAND STANFORD JUNIOR UNIVERSITY. ALL RIGHTS RESERVED.
*************************************************************************/
#include "smooth.h"

	static int MM, nn;
	static VECT ys, yr;


/*
* RETURN PARTICULAR MSE
*/
static OPT_FUNC func_mse(ALPHA_SET alpha)
{
	int n;
	VECT pys = ys+1, pyr = yr+1;
	double mse = 0.;

	if ( !vs_setup_alpha(alpha) )			ABORT("8.")
	if ( !vs_smooth(ys, 1) )			ABORT("9.")

	for (n=nn*MM; n >= 1; n--, pys++, pyr++)
		mse += (*pys - *pyr) * (*pys - *pyr);

	return (OPT_FUNC)( mse / nn);
}


/*
*	INPUT:
*		yr[1..nxM]		true values
*		alpha[1..M]		initial value
*		tol:			termination tolerance
*	OUTPUT:
*		ys:		best estimates
*	EFFECT:
*		alpha:		optimal alpha
*	RETURNS:
*				best mse
*/
double best_mse(CMAT l_ys, CMAT l_yr, ALPHA_SET alpha,
VECT ti, CMAT yi, COVAR *covar, ctruf flag_iid, cint M, cint N, cdouble tol)
{
	OPT_ARG mse;
	DIRSET dirset;
	int iter;

	// "SAVE" THE PASSED VARIABLES
	ys = l_ys;	yr = l_yr;	MM = M;		nn = N;

	// SETUP: INDEPENDENT OF ALPHA
	if (!vs_setup_alloc(MM, nn))		ABORT("1.")
	if (!vs_setup_ti(ti, nn))		ABORT("2.")
	if (!vs_setup_yi(yi, MM, nn))		ABORT("3.")
	if (!vs_setup_cov(covar, flag_iid))	ABORT("4.")

	// THE ROUTINES BELOW ALL DEPEND ON ALPHA

	dirset = SETUP_DIRSET(MM);
	powell(alpha, dirset, MM, (OPT_ARG) tol, &iter, &mse, func_mse);

	vs_cleanup();

	return (double) mse;
}
