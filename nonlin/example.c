/*********************** example() ***************************************
*	Simple example using NONLIN subroutines.
*	This subroutine solves the following problem:
*	\hat{g} = \arg \min_g \sum_{n=1}^N \| y_n - h(g(t_n)) \|^2 +
*	\sum_{m=1}^M \int \ddot{g}^2 dt.
*	This is a penalized maximum likelihood problem, and it reverts
*	to the ordinary MLE when alpha==0.  However, since the vector
*	spline smoothing subroutine called below actually uses Nexp(alpha)
*	for alpha, if the MLE is desired, then set alpha to be NULL.
*	The function "h" is passed in the structure prob.
*************************************************************************/
#include "nonlin.h"

#if 0
	void nonlin_est();

	extern void nl_alloc(), nl_setup(), nl_init(), nl_rough(), nl_cleanup();
	extern sof nl_iter();
	extern Nonlin *nl_info();
	extern double nl_score(), nl_score_old();
#endif


/*
*	Non-linear smoothing iterative estimation
*	INPUT:
*		alpha[1..M]		smoothing parameters
*		ti[1..N]		sample points (or null for uniform)
*		meas[1..N][1..L]	original noisy measurements
*		istates[1..N][1..M]	initial estimates
*		prob			problem dependent structure
*		score_id		one of JOBUR JOBCV JOBGCV
*		L			length of each measurement vector
*		M			# of unknowns
*		N			# of sample points
*		O			# of offsets
*		W			maximum signal width
*	EFFECT:
*		score			a cross validation score
*		ostates			final estimates
*/

static void nonlin_est(
double *score,
PARAM_SET ostates,
ALPHA_SET alpha,
VECT ti,
MEAS meas,
PARAM_SET istates,
Problem *prob,
int score_id, int LL, int MM, int nn, int PP, int WW,
cchar *estfile)
{
	Nonlin *init;
	truf flag_mle = alpha == 0;
	truf flag_mod = prob->modify != 0;

	nl_alloc(LL, MM, nn, PP, WW, flag_mle, flag_mod);
	nl_setup(ti, nn);
	init = nl_info();
	nl_init(init, istates, meas, prob);
	nl_rough(init, alpha);

	/*
	* Iterate to find best estimate
	*/
	if (!nl_iter(init, alpha, meas, prob, estfile))
		message("Warning: iterations not succesful", 0.);

	if (alpha) {
		double s;
		message("score new: ", s = nl_score(alpha, init, score_id, False));
		message("score lam: ", nl_score(alpha, init, score_id, True));
		message("score old: ", nl_score_old(alpha, init, score_id));
		*score = s;
	}

	/*
	*	Save results
	*/
	INFO_TO_PARAM(ostates, init, MM, nn);

	nl_cleanup();
	free_info(init, nn);
}
