/*********************** nonlin_best.c ***********************************
*	NONLINEAR SMOOTHING VECTOR MEASUREMENTS WITH AUTOMATIC SCORE CRITERIA
*************************************************************************/
#include "nonlin.h"

#if 0
	double nonlin_best(), nl_best();
	static OPT_FUNC score_func();

	extern void nl_alloc(), nl_setup(), nl_init(), nl_rough(), nl_cleanup();
	extern sof nl_iter();
	extern Nonlin *nl_info(), *info_alloc();
	extern double nl_score();
	extern void free_info(), copy_info();
#endif

	static MEAS meas;
	static int score_id;
	static Problem *prob;
	static Nonlin *init, *work;
	static int MM, nn, PP, WW;
#define flag_mle	False
#define ESTFILE		((char *) 0)


/*
*	INPUT:
*		alpha[1..M]:		initial value
*		ti[1..N]:		sample points
*		meas[1..N][1..L]:	measurements
*		xi[1..N][1..M]:		initial value
*		prob			problem dep. struct
*		score_id:		which score to do
*		tol:			termination tolerance
*	OUTPUT:
*		alpha[1..M]:		optimal alpha
*		xs[1..N][1..M]:		best estimates
*	RETURNS:			best score
*/

double nonlin_best(
PARAM_SET xs,
ALPHA_SET alpha,
VECT ti,
MEAS meas,
PARAM_SET xi,
Problem *prob,
int score_id,
double tol;
int LL, int MM, int nn, int PP, int WW)
{
	double score;
	truf flag_mod = prob->modify != 0;
	Nonlin *info;

	/*
	*	Setup for nonlinear estimation
	*/
	nl_alloc(LL, MM, nn, PP, WW, flag_mle, flag_mod);
	nl_setup(ti, nn);
	info = nl_info();
	nl_init(info, xi, meas, prob);

	/*
	*	Do score minimization
	*/
	score = nl_best(info, alpha, meas, prob, score_id, tol, MM, nn, PP, WW);

	/*
	*	Save final estimates
	*/
	INFO_TO_PARAM(xs, info, MM, nn);
	free_info(info, nn);

	return score;
}


/*
*	Must call nl_alloc() and nl_setup() before using this routine
*	INPUT:
*		info		mostly precomputed
*	OUTPUT:
*		info		best estimate ????
*		alpha		best value
*	RETURNS:		best score
*/

double nl_best(
Nonlin *info;
ALPHA_SET alpha;
MEAS l_meas;
Problem *l_prob;
int l_score_id;
double tol;
int lMM, int lnn, int lPP, int lWW)
{
	OPT_ARG score, score2;
	DIRSET dirset;
	int iter;

	/*
	*	save variables needed by score_func()
	*/
	MM = lMM; nn = lnn; PP = lPP; WW = lWW;
	score_id = l_score_id;		prob = l_prob;		meas = l_meas;
	if (!alpha)
		ABORT("Must use nonzero alpha for CV score")

	/*
	*	Do initial estimation
	*/
	work = nl_info();
	init = nl_info();
	copy_info(info, init, MM, nn, PP, WW);
	nl_rough(init, alpha);
	if ( !nl_iter(init, alpha, meas, prob, ESTFILE) )
		message("WARNING: initial estimate doesn't converge!", 0.);

	/*
	*	Do nonlinear minimization
	*/
	dirset = SETUP_DIRSET(MM);
	powell(alpha, dirset, MM, tol, &iter, &score, score_func);

	/*
	*	Now use that best alpha value
	*/
	copy_info(init, info, MM, nn, PP, WW);
	nl_rough(info, alpha);
	if ( !nl_iter(info, alpha, meas, prob, ESTFILE) )
		message("WARNING: final estimate doesn't converge!", 0.);
	score2 = nl_score(alpha, info, score_id, False);
	if (score2 != score) {
		message0("WARNING: score1: ", score);
		message("and score2 differ:", score2);
	}

	FREE_MATRIX(dirset);
	free_info(init, nn);
	free_info(work, nn);
	nl_cleanup();

	return (double) score;
}


/*
*	RETURN PARTICULAR SCORE
*/

static OPT_FUNC score_func(ALPHA_SET alpha)
{
	double score;
	/*
	*	Since nl_iter destroys work, we must make a copy`
	*/
	copy_info(init, work, MM, nn, PP, WW);
	nl_rough(work, alpha);

	if ( !nl_iter(work, alpha, meas, prob, ESTFILE) )
		message("Warning: bad convergence", 0.);

	score = nl_score(alpha, work, score_id, False);
	message0("score,alpha: ", score);      PRINT_ALPHA(alpha, MM);
	return (OPT_FUNC) score;
}
