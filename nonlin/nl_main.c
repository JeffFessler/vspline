/*********************** nl_main.c ***************************************
*	Non-parametric estimation of a smooth vector-valued function
*	observed through a non-linearity with additive white noise.
*	USES log10(), floor()
*************************************************************************/
#include <math.h>
#include "nonlin.h"
#define MAXITER		30

	void nl_alloc(), nl_setup(), nl_init(), nl_cleanup(), nl_rough();
	sof nl_iter();
	Nonlin *nl_info();
	double nl_score(), nl_score_old();

	extern sof vs_setup_alloc(), vs_setup_ti();
	extern sof cv_setup_alloc(), penalty();
	extern void vs_cleanup(), cv_cleanup();
	extern Nonlin *alloc_info();	extern void copy_info(), free_info();
	extern sof compute_sigjac();
	extern void compute_hessian(), compute_resid();
	extern double compute_score(), compute_score_old();
	extern sof nonlin_update();

	static int LL, MM, nn, PP, WW, nnI;
	static Nonlin *work;
	static truf flag_mle, flag_mod;
	static double Lambda;		/* Marquart relaxation parameter */


/*
*	SETUP DIMENSIONS AND ALLOCATE SPACE FOR INFO STRUCTURES
*	INPUT:
*		L			length of measurement vectors
*		M			# of unknowns
*		N			# of sample points
*		P			# of offsets
*		W			max signal width
*	EFFECT:
*		work, ...
*/

void nl_alloc(L, M, N, P, W, l_flag_mle, l_flag_mod)
int L, M, N, P, W;
truf l_flag_mle, l_flag_mod;
{
	LL = L; MM = M;	nnI = N; PP = P; WW = W;
	flag_mle = l_flag_mle;
	flag_mod = l_flag_mod;

	/*
	*	Allocate space for working state
	*/
	work = alloc_info(MM, nnI, PP, WW, flag_mod);
	if (!flag_mle)
	{
		if (!vs_setup_alloc(MM, nnI))		ABORT("ne1")
		if (!cv_setup_alloc())			ABORT("ne3")
	}
}


/*
*	SETUP ti's etc. for vector spline smoothing
*	INPUT:
*		ti[1..N]		sample points (or null for uniform)
*/

void nl_setup(ti, N)
VECT ti;
int N;
{
	nn = N;
	if (nn > nnI)	ABORT("smoothing too many")

	/*
	*	Setup for spline smoothing
	*/
	if (!flag_mle)
		if (!vs_setup_ti(ti, nn))		ABORT("ne2")
}


/*
*	ALLOCATE SPACE FOR INITIAL INFO
*/

Nonlin *nl_info()
{
	return alloc_info(MM, nnI, PP, WW, flag_mod);
}


/*
*	SETUP INITIAL ESTIMATES
*	INPUT:
*		istates[1..N][1..M]	initial estimates
*		meas[1..N][1..L]
*		prob			problem dependent struct
*		N:			number of measurements
*	EFFECT:
*		cv			cross validation score
*		ostates			final estimates
*/

void nl_init(info, istates, meas, prob)
Nonlin *info;
PARAM_SET istates;
MEAS meas;
Problem *prob;
{
	/*
	*	Convert initial estimates to "info" format
	*/
	PARAM_TO_INFO(info, istates, MM, nn);

	/*
	*	Compute the (total) measurement norm
	*/
	info->bad.meas_norm = NORM_MEAS(meas, nn, LL);

	/*
	*	Compute initial signals, jacobians, residuals, and hessian
	*/
	if (!compute_sigjac(info, prob, MM, nn, PP, WW))
		ABORT("The initial estimate doesn't sigjac!")
	compute_resid(info, meas, nn, PP, WW);
	compute_hessian(info, MM, nn, PP, WW);
}


/*
*	EVERYTHING ABOVE HERE IS INDEPENDENT OF ALPHA, BELOW IS DEPENDENT
*/


/*
*	COMPUTE ROUGHNESS PENALTY
*/

void nl_rough(info, alpha)
Nonlin *info;
ALPHA_SET alpha;
{
	if (alpha) {
	    if (!penalty(&info->bad.roughness, alpha, info->par[1], MM, nn))
		ABORT("rough")
	}
	else {
		message("WARNING: roughness for mle?", 0.);
	}
}


/*
*	ITERATIVE ESTIMATION, GIVEN PRECOMPUTED:
*	info-> roughness, meas_norm, norm_diff, sig_mod, resid, jac, off, xh, covar
*	BASICALLY, EVERYTHING IS PRECOMPUTED!
*	INPUT
*		best			initial estimate
*		alpha[1..M]		smoothing parameters
*		meas[1..N][1..L]	measurements
*		prob			problem dependent stuff
*		estfile			print estimates?
*	OUTPUT
*		best			final estimate
*	EVERYTHING HAS BEEN COMPUTED ON OUTPUT AS WELL (FIX: verify?)
*	RETURNS				successful descent?
*/

sof nl_iter(best, alpha, meas, prob, estfile)
Nonlin *best;
ALPHA_SET alpha;
MEAS meas;
Problem *prob;
char *estfile;
{
	double fbst, fnew;
	Badness bbst, bnew;
	int ii, decade;

	/* BADNESS CRITERION	*/
	bbst = best->bad;
	fbst = bbst.meas_norm+bbst.norm_loc_meas+bbst.norm_loc_resid+bbst.roughness;

	/* PRESET LM PARAMETER	*/
	Lambda = 0.001;

	if (prob->print) {
		message0("", -1.);
		message0("", floor(0.5+log10(Lambda)));
		message0("", bbst.meas_norm);		message0("", bbst.norm_loc_meas);
		message0("", bbst.norm_loc_resid);	message0("", bbst.roughness);
		message("\t", fbst);
	}

	/*
	*	Loop over iterations
	*/
	for (ii = 0; ii < MAXITER; ii++)
	{
		double pcnt;
		truf improve;

		decade	= (int) floor(0.5+log10(Lambda));
		improve = nonlin_update(best, work, meas, prob, alpha, &Lambda, MM, nn, PP, WW);

		if (improve && estfile)
			FPRINT_MATRIX(estfile, &ii, best->par, nn, MM);

		bnew = (improve) ? best->bad : work->bad;
		fnew = bnew.meas_norm+bnew.norm_loc_meas+bnew.norm_loc_resid+bnew.roughness;
		pcnt = BADDIFF(bnew, bbst) / fbst * 100.;

		if (improve) {
			fbst = fnew;	bbst = bnew;
		}

		if (prob->print) {
			message0("", (double) ii);
			message0("", (double) decade);
			message0("", bnew.meas_norm);		message0("", bnew.norm_loc_meas);
			message0("", bnew.norm_loc_resid);	message0("", bnew.roughness);
			message0("\t", fnew);		message((improve) ? "\tY" : "\tN  ", pcnt);
		}

		/*
		*	Done if improvement is small
		*/
		if (improve && pcnt > -0.1 && decade < 5)
			Ok
		if (decade > 8)
			Fail("decade > 8")
	}
	Fail("nl_iter")
}


/*
*	SMOOTH FINAL TO COMPUTE CV SCORE
*/

double nl_score_old(alpha, best, score_id)
ALPHA_SET alpha;
Nonlin *best;
int score_id;
{
	return compute_score_old(best, alpha, score_id, MM, nn);
}


/*
* Use final estimate to compute approximate cv-like score
*/

double nl_score(
ALPHA_SET alpha,
Nonlin *best,
int score_id,
truf flag_lam)
{
	return compute_score(best, work, alpha, score_id, (flag_lam) ? Lambda : 0., MM, nn);
}


/*
* Free everything that has been allocated
*/

void nl_cleanup(void)
{
	free_info(work, nnI);
	if (!flag_mle)
	{
		vs_cleanup();
		cv_cleanup();
	}
	LL = MM = nn = PP = WW = nnI = 0;
}
