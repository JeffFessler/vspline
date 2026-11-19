/*********************** lm_test.c ***************************************
* TEST LEVENBERG-MARQUARDT ALGORITHM.
* As lambda -> Infinity, does xhat_lambda -> x0 ?
*************************************************************************/
#include <math.h>
#include "nonlin.h"

	void lm_test(), lm_test_info();

	extern void nl_alloc(), nl_setup(), nl_init(), nl_cleanup();
	extern Nonlin *nl_info();

	extern sof penalty();
	extern sof nonlin_update();
	extern Nonlin *alloc_info();	extern void free_info();


/*
*	Non-linear smoothing iterative estimation
*	INPUT:
*		alpha[1..M]		smoothing parameters
*		ti[1..N]		sample points (or null for uniform)
*		meas[1..N][1..L]	original noisy measurements
*		istates[1..N][1..M]	initial estimates
*		prob
*		L			length of each measurement vector
*		M			# of unknowns
*		N			# of sample points
*		O			# of offsets
*		W			maximum signal width
*		lambda			specified
*	EFFECT:
*		ostates			final estimates
*/

void lm_test(ostates, alpha, ti, meas, istates, prob,
	LL, MM, NN, OO, WW, estfile, lambda)
PARAM_SET ostates, istates;
ALPHA_SET alpha;
VECT ti;
MEAS meas;
Problem *prob;
int LL, MM, NN, OO, WW;
char *estfile;
double lambda;
{
	Nonlin *init;
	truf flag_mle = alpha == 0;
	truf flag_mod = prob->modify != 0;

	nl_alloc(LL, MM, NN, OO, WW, flag_mle, flag_mod);
	nl_setup(ti, NN, flag_mle);
	init = nl_info();
	nl_init(init, istates, meas, prob);

	lm_test_info(ostates, alpha, meas, init, prob, estfile, LL, MM, NN, OO, WW, lambda);

	nl_cleanup(flag_mle);
	free_info(init, NN);
}


/*
*	INPUT:
*		alpha[1..M]		smoothing parameters
*		meas[1..N][1..L]
*		istates[1..N][1..M]	initial estimates
*		init			initial info structure
*		prob
*	OUTPUT:
*		ostates			final estimates
*		init			final (best) info structure
*/

void lm_test_info(ostates, alpha, meas, init, prob,
	estfile, LL, MM, NN, OO, WW, lambda)
PARAM_SET ostates;
ALPHA_SET alpha;
MEAS meas;
Nonlin *init;
Problem *prob;
char *estfile;
int LL, MM, NN, OO, WW;
double lambda;
{
	static Nonlin *best, *inew, *work;
	int decade;
	double pcnt=0., fnew, fold;
	truf change;

	/* COPY INITIAL TO BEST		*/
	best = init;
	work = nl_info();

	best->meas_norm = NORM_MEAS(meas, NN, LL);

	/* COMPUTE INITIAL ROUGHNESS PENALTY AND BADNESS CRITERION	*/
	if (alpha)
	    if (!penalty(&best->roughness, alpha, best->par[1], MM, NN))
		ABORT("ne4");
	fold = best->meas_norm + best->norm_diff + best->roughness;

	message0("", -1.); message0("", floor(0.5+log10(lambda)));
	message0("", best->norm_diff); message0("", best->roughness);
	message("", best->norm_diff+best->roughness);

	change =
	nonlin_update(best, work, meas, prob, alpha, &lambda, MM, NN, OO, WW);

	inew	= (change) ? best : work;
	fnew	= inew->meas_norm + inew->norm_diff + inew->roughness;
	pcnt	= (fnew - fold) / fold * 100.;
	decade	= (int) floor(0.5+log10(lambda));

	message0("", 0.);	message0("", (double) decade);
	message0("", inew->meas_norm + inew->norm_diff);	message0("", inew->roughness);
	message0("", fnew);
	if (change) message("Y", pcnt); else	message("N (", pcnt);

	if (change)
	{
		Nonlin tmp;
		tmp = *best;
		*best = *work;
		*work = tmp;
	}

	free_info(work, NN);
}
