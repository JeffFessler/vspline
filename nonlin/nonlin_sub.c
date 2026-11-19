/***************************** nonlin_sub.c ******************************
*	Data type dependent subroutines for nonlin_*.c
*************************************************************************/
#include "nonlin.h"

	sof compute_sigjac();
	void compute_resid();
	void compute_hessian();
	sof compute_ridge();
	double compute_score(), compute_score_old();
	static void borrow();

	extern void hessian();
	extern sof ridge_reg1();
	extern sof chol_sym_band();
	extern void invert_sym_band(), copy_covar(), scale_covar();
	extern sof vs_setup_yi(), vs_setup_cov(), vs_setup_alpha();
	extern sof vs_smooth(), cv_setup_covar(), cross_val();


/*
*	COMPUTE SIGNALS AND JACOBIANS FOR THIS SET OF ESTIMATES
*	INPUT:
*		prob:		problem dependent stuff
*/

sof compute_sigjac(Nonlin *info, Problem *prob, int MM, int nn, int PP, int WW)
{
	SIGJAC_RETURN out;

	(*prob->sigjac)(&out, info->sig_raw, info->jac_raw, info->off, info->par, prob->sigarg, MM, nn, PP, WW);
	if (out)
	{
		message0("start1 ", (double) info->off[out].start[1]);
		message0("count1 ", (double) info->off[out].count[1]);
		PRINT_VECTOR(info->par[out], MM);
		Fail1("Bad parameters %g", (double) out)
	}

	if (prob->modify) {
		COPY_SIGNAL_SET(info->sig_raw, info->sig_mod, nn, PP, WW);
		COPY_JACOB_SET(info->jac_raw, info->jac_mod, nn, MM, PP, WW);
		(*prob->modify)(info->sig_mod, info->jac_mod, info->off, info->par, prob->modarg, PP, WW);
	}

	Ok
}


/*
*	Compute residual and norms for the current signal estimate.
*	This is tricky, since it involves both the full measurement
*	and the piecemeal signals.
*	|Y-S| = |Y|(s=0) + |Y-S|(s!0) = |Y| - |Y|(s!0) + |Y-S|(s!0)
*/

void compute_resid(info, meas, nn, PP, WW)
Nonlin *info;
MEAS meas;
int nn, PP, WW;
{
	register int jj, pp, ww;

	info->bad.norm_loc_meas = 0.;
	info->bad.norm_loc_resid = 0.;
	COPY_SIGNAL_SET(info->sig_mod, info->resid, nn, PP, WW);

	for (jj=nn; jj >= 1; jj--) {
	    for (pp=PP; pp >= 1; pp--) {

		SIGNAL_TYPE *p1 = &meas[jj][info->off[jj].start[pp]];
		SIGNAL_TYPE *p2 = &info->resid[jj][pp][1];

		for (ww=info->off[jj].count[pp]; ww >= 1; ww--, p1++, p2++) {
			*p2 = *p1 - *p2;
			info->bad.norm_loc_resid += *p2 * *p2;
			info->bad.norm_loc_meas -= *p1 * *p1;
		}
	    }
	}
}


/*
*	Compute Hessian
*	INPUT:
*/

void compute_hessian(info, MM, nn, PP, WW)
Nonlin *info;
int MM, nn, PP, WW;
{
	register int ii;

	for (ii = nn; ii >= 1; ii--)
		hessian(info->covar[ii].inv, info->xh[ii], info->jac_mod[ii],
			info->resid[ii], MM, PP, WW);
}


/*
*	Compute ridge regression
*	If a jacobian is rank deficient, "borrow" one from a
*	neighboring slice, but scale it down by JAC_SCALE
*	INPUT:
*	RETURN:		jacobians computed ok?
*/
#define JAC_SCALE	0.001

sof compute_ridge(new, old, lambda, MM, nn, flag_borrow)
Nonlin *new, *old;
double lambda;
int MM, nn;
truf flag_borrow;
{
	int jj;
	truf flag_bad_start = False;

	for (jj = 1; jj <= nn; jj++)
	{
	    if (!ridge_reg1(new->par[jj], &new->covar[jj], old->covar[jj].inv,
		old->par[jj], old->xh[jj], lambda, True, MM))
	    {
		PRINT_VECTOR(old->par[jj], MM);
		if (flag_borrow == False)
			Fail1("Bad jacobian: jj=%d", jj)
		if (jj == 1)
			flag_bad_start = True;
		else if (!flag_bad_start)
			borrow(new, old, jj-1, jj, MM, JAC_SCALE);
	    }
	    else if (flag_bad_start)
	    {
		int ii;
		for (ii = 1; ii < jj; ii++)
			borrow(new, old, jj, ii, MM, JAC_SCALE);
		flag_bad_start = False;
	    }
	}

	/*
	*	If borrowing, make sure a good one was found
	*/
	if (flag_borrow && flag_bad_start)
		ABORT("Sorry - this entire object was rank deficient")
	Ok
}


/*
*	Borrow a good jacobian to replace a rank deficient one.
*	We hope this procedure is never called!
*	INPUT
*		scale:		scale the jacobian
*				inv(cov) = J'*J
*/

static void borrow(info, old, n_good, n_bad, MM, scale)
Nonlin *info, *old;
int n_good, n_bad, MM;
double scale;
{
	COVAR *src = &info->covar[n_good];
	COVAR *dst = &info->covar[n_bad];

	message0("WARNING: borrow from ", (double) n_good);
	message(" to ", (double) n_bad);

	copy_covar(src, dst, MM);
	scale_covar(dst, 1./(scale*scale), MM);

	/* FIX: SHOULD THIS BE old OR info ?			*/
	COPY_VECTOR(old->par[n_bad], info->par[n_bad], MM);
}


/*
*	Compute approximate cross-validation like score.
*	INPUT:
*		info:	optimal nonlinear estimate for this alpha
*			and hessian (ie, inv and xh) are up-to-date
*/

double compute_score(info, work, alpha, score_id, lambda, MM, nn)
Nonlin *info, *work;
ALPHA_SET alpha;
double lambda;
int score_id, MM, nn;
{
	double ur, cv, gcv;

	/*
	*	First, generate linearized measurement
	*/
	if (!compute_ridge(work, info, lambda, MM, nn, False))
		ABORT("Oh no! Can't ridge!")

	if (!vs_setup_yi(work->par[1], MM, nn))		ABORT("")
	if (!vs_setup_alpha(alpha))			ABORT("")
	if (!vs_setup_cov(work->covar, False))		ABORT("")
	if (!vs_smooth((CMAT) 0, 0))			ABORT("")

	if (!cv_setup_covar())				ABORT("")
	if (!cross_val(&ur, &cv, &gcv, score_id))	ABORT("")

	if (score_id == JOBUR)		return ur;
	if (score_id == JOBCV)		return cv;
	if (score_id == JOBGCV)		return gcv;
	ABORT("Invalid score id");	return 0.;
}


/*
*	COMPUTE CROSS VALIDATION SCORE FOR info
*	FIX: THIS MODIFIES info
*	and - I think it is wrong.
*	INPUT:
*/

double compute_score_old(info, alpha, score_id, MM, nn)
Nonlin *info;
ALPHA_SET alpha;
int score_id, MM, nn;
{
	int ii;
	double cv, dum;

	for (ii = nn; ii >= 1; ii--) {
		COVAR *c = &info->covar[ii];
		if ( !chol_sym_band(c->Linv, c->Dinv, c->inv, MM-1, MM) ) ABORT("")
		invert_sym_band(c->cov, c->Linv, c->Dinv, MM-1, MM);
	}

	if (!vs_setup_yi(info->par[1], MM, nn))		ABORT("")
	if (!vs_setup_alpha(alpha))			ABORT("")
	if (!vs_setup_cov(info->covar, False))		ABORT("")
	if (!vs_smooth(info->par[1], 4))		ABORT("")

	if (!cv_setup_covar())				ABORT("")
	if (!cross_val(&dum, &cv, &dum, score_id))	ABORT("")

	return cv;
}
