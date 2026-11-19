/*********************** nonlin_update() *********************************
* A single Levenberg-Marquardt update of estimate
*************************************************************************/
#include "nonlin.h"

	sof nonlin_update();

	extern sof compute_ridge(), compute_sigjac();
	extern void compute_resid(), compute_hessian();
	extern void copy_info();
	extern sof vs_setup_yi(), vs_setup_cov(), vs_setup_alpha(), vs_smooth();
	extern double roughness();


/*
* COMPUTE NEW PARAMETER ESTIMATES
* During this routine, every element of "work" must be computed,
* so that if it turns out to be better, we can transfer it to
* "best" - which must be complete on return.
* INPUT:
*	best:			previous estimate
*	work:			work space
*	meas[1..N][1..L]:	measurements
*	alpha[1..M]		smoothing parameters
*	lambda:			lm parameter
*	prob			problem dependent
* EFFECT:
*	best:			set to new estimate
* RETURNS:
*	1 IF NEW ESTIMATE BETTER, 0 IF NO CHANGE (EXCEPT LAMBDA)
*	also return 0 if new estimate doesn't satisfy constraints
*/

sof nonlin_update(
Nonlin *best,
Nonlin *work,
MEAS meas,
Problem *prob,
ALPHA_SET alpha,
double *lambda,
int MM, int nn, int PP, int WW)
{
	/*
	* DO RIDGE REGRESSION TO GET NEW ESTIMATES.
	* work: (par, covar)
	*/
	if (!compute_ridge(work, best, *lambda, MM, nn, prob->borrow)) {
		work->bad.norm_loc_meas = work->bad.norm_loc_resid = work->bad.roughness = 0.;
		*lambda *= 10.;
		return Failure;
	}

	/*
	* SMOOTH ESTIMATES
	* work: (par, roughness)
	*/
	if (alpha)
	{
		if (!vs_setup_yi(work->par[1], MM, nn))		ABORT("nu1")
		if (!vs_setup_alpha(alpha))			ABORT("nu2")
		if (!vs_setup_cov(work->covar, False))		ABORT("nu3")

		if (!vs_smooth(work->par[1], 4))		ABORT("nu4")
		work->bad.roughness = roughness();
	}

	/*
	* COMPUTE NEW SIGNALS AND JACOBIANS
	* work: (sig_raw, sig_mod, jac_raw, jac_mod, off)
	*/
	if (!compute_sigjac(work, prob, MM, nn, PP, WW))
	{
		work->bad.norm_loc_meas = work->bad.norm_loc_resid = work->bad.roughness = 0.;
		*lambda *= 10.;
		return Failure;
	}

	/*
	* Compute new residuals and norm
	* work: (resid, norm_diff)
	*/
	compute_resid(work, meas, nn, PP, WW);

	/*
	* Must transfer meas_norm to work since it hasn't been computed above.
	*/
	work->bad.meas_norm = best->bad.meas_norm;

	/*
	* If new estimate is better, then save it.
	* If not, increment relaxation parameter lambda and try again.
	*/
	if ( BADDIFF(work->bad, best->bad) < 0. )
	{
		*lambda *= 0.1;
		copy_info(work, best, MM, nn, PP, WW);

		/*
		* The remaining computation need to leave best complete
		* best: (covar.inv, xh)
		*/
		compute_hessian(best, MM, nn, PP, WW);
		Ok
	}

	else
	{
		*lambda *= 10.;
		return Failure;
	}
}
