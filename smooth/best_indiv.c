/*********************** best_indiv.c ************************************
* SMOOTH CURVES OPTIMALLY - INDIVIDUALLY
* Jeff Fessler
* COPYRIGHT 1989, THE BOARD OF TRUSTEES OF THE
* LELAND STANFORD JUNIOR UNIVERSITY. ALL RIGHTS RESERVED.
*************************************************************************/
#include "smooth.h"

	static int score_id;


/*
* RETURN PARTICULAR SCORE
*/
static OPT_FUNC score_func(ALPHA pass_alpha)
{
	double ur, cv, gcv;
	double score = 0.;
	ALPHA alpha = pass_alpha;
	ALPHA_SET pseudo_alpha = &alpha-1;		/* tricky...	*/

	if ( !vs_setup_alpha(pseudo_alpha) )		ABORT("sf 1")
	if ( !vs_smooth((VECT) 0, 0) )			ABORT("sf 2")
	if ( !cross_val(&ur, &cv, &gcv, score_id) )	ABORT("sf 3")

	if (score_id == JOBUR)		score = ur;
	else if (score_id == JOBCV)	score = cv;
	else if (score_id == JOBGCV)	score = gcv;
	else				ABORT("Invalid score id")

	return (OPT_FUNC) score;
}


/*
*	INPUT:
*		score_id:	which score to do
*		alpha:		initial value
*		tol:		termination tolerance
*	OUTPUT:
*		alpha:		optimal alpha
*		ys:		best estimates
*	RETURNS:
*				best score
*/
static double best_single(VECT ys, ALPHA *alpha, VECT yi, cint nn, cdouble tol)
{
	double score;
	OPT_ARG ax, bx, cx, fa, fb, fc;

	/* "SAVE" THE PASSED VARIABLES	*/
	if (!vs_setup_yi(yi, 1, nn) )			ABORT("bs 1")


	/* THE ROUTINES BELOW ALL DEPEND ON ALPHA			*/

	ax = *alpha-1;
	bx = *alpha;
	mnbrak(&ax, &bx, &cx, &fa, &fb, &fc, score_func);
	score = brent(ax, bx, cx, score_func, tol, alpha);

	if ( !vs_setup_alpha(alpha-1) )			ABORT("bs 2")
	if ( !vs_smooth(ys, 1) )			ABORT("bs 3")

	return score;
}


/*
*	SETUP FAKE 1x1 COVARIANCE MATRIX
*	RETURNS:
*		covar[1]
*/
static COVAR *fake_covar(void)
{
	DATA one = 1.0;
	VECT isig = &one - 1;

	return inv_to_covar_set(isig, 1, 1, True);
}


/*
*	CURRENTLY ONLY IMPLEMENTED FOR iid CASE
*	INPUT:
*		alpha[1..M]		initial values
*		ti[1..N]		sample points
*		yi[1..nxM]		NOISY MEASUREMENTS TO BE SMOOTHED
*		score_id:		WHICH SCORE
*		tol:			CONVERGENCE TOLERANCE
*	OUTPUT:
*		ys[1..nxM]		best estimates
*		trace[1..M]		trace of influence matrices
*	EFFECT:
*		alpha:			best alpha values
*	RETURNS:			sum of scores
*/
double best_indiv(VECT ys, VECT trace, ALPHA_SET alpha, VECT ti, VECT yi,
cint l_score_id, cdouble tol, cint MM, cint nn)
{
	COVAR *covar = fake_covar();
	double score_sum = 0.;

	// SETUP: INDEPENDENT OF ALPHA
	if (!vs_setup_alloc(1, nn))		ABORT("bi 1")
	if (!vs_setup_ti(ti, nn))		ABORT("bi 2")
	if (!vs_setup_yi(yi, 1, nn))		ABORT("bi 3")
	if (!vs_setup_cov(covar, 1))		ABORT("bi 4")
	if (!cv_setup_alloc())			ABORT("bi 5")
	if (!cv_setup_covar())			ABORT("bi 6")

	score_id = l_score_id;

	// LOOP THROUGH INDIVIDUAL VECTOR COMPONENTS
	for (int m=1; m <= MM; m++) {
		score_sum += best_single(ys+(m-1)*nn, &alpha[m], yi+(m-1)*nn, nn, tol);
		if (!trace_i_a(&trace[m], alpha+m-1))	ABORT("bi 7")
	}

	covar_free(covar, 1);
	vs_cleanup();
	cv_cleanup();

	return score_sum;
}
