/************************** ridge_reg1.c *********************************
*	RIDGE REGRESSION
*************************************************************************/
#include "nonlin.h"

	void hessian();
	sof ridge_reg1();

	extern sof chol_sym_band();
	extern void solve_sym_band(), invert_sym_band();
	extern void scale_covar();


/*
*	RIDGE REGRESSION
*	Xr = inv(HH' + lam diag(HH')) H (Y - h(X0)) + X0
*	INPUT:
*		inv: [1..M][1..M]	H H'
*		x0: [1..M]		X0
*		xh: [1..M]		xh = H (Y - h(X0))
*		lambda:			ridge coef. (levenberg marquardt param)
*		flag_inv:		compute inverse?
*	OUTPUT:
*		xr: [1..M]		xr = inv(hh' hh + lam) xh + x0
*		covar:			error covariance matrices
*	RETURNS:			success (not singular)?
*/

#define SCALE_COVAR	10000.

sof ridge_reg1(
VECT xr,
COVAR *covar,
SYMBAND inv,
VECT x0,
VECT xh,
double lambda,
truf flag_inv,
int MM)
{
	register int mm;

	// Copy inv to covar
	COPY_MATRIX(inv, covar->inv, MM, MM);

	// Add ridge parameter to diagonal elements
	for (mm=MM; mm >= 1; mm--)
		covar->inv[mm][MM] *= (1. + lambda);

	/*
	* Do LDL decomposition of "inverse" (hessian) matrix.
	* If it fails, try ADDING lambda to diagonal.
	*/
	if ( !chol_sym_band(covar->Linv, covar->Dinv, covar->inv, MM-1, MM) )
	{
		for (mm=MM; mm >= 1; mm--)
			covar->inv[mm][MM] += lambda;

		if ( !chol_sym_band(covar->Linv, covar->Dinv, covar->inv, MM-1, MM) )
			return Failure;
		else
			message("WARNING: Added this lambda: ", lambda);

		// I don't have too much confidence in this measurement now,
		// so scale up the covariance (uncertainty). FIX: side effects?
		scale_covar(covar, SCALE_COVAR, MM);
	}

	// Compute estimate "delta"
	solve_sym_band(xr, covar->Linv, covar->Dinv, xh, MM-1, MM, 0);

	// Add previous value to get new estimate
	INC_VECTOR(xr, x0, MM);

	if (flag_inv)
		invert_sym_band(covar->cov, covar->Linv, covar->Dinv, MM-1, MM);

	Ok
}


/*
*	COMPUTE HESSIAN
*	Note: this function relies on the fact that wherever h(X) is 0, H is
*	also zero, so H * (Y - h(X0)) is correct, even though only the
*	non-zero elements of Y-h(X0) are computed in computed_resid().
*	inv = H*H'
*	INPUT:
*		hh: [1..M][1..O][1..W]		H
*		resid: [1..O][1..W]		Y - h(X0)
*	OUTPUT:
*		inv: [1..M][1..M]		hh * hh'
*		xh: [1..M]			xh = hh' resid
*/

void hessian(inv, xh, hh, resid, MM, PP, WW)
SYMBAND inv;
VECT xh;
JACOBIAN hh;
SIGNAL resid;
int MM, PP, WW;
{
	register int ii, jj;

	for (ii=MM; ii >= 1; ii--)
	{
		/*
		*	Compute H resid
		*/
		xh[ii] = DOT_JACOB(hh[ii][1], resid[1], PP*WW);

		/*
		*	Compute H H'
		*/
		for (jj=ii; jj >= 1; jj--)
			inv[ii][jj-ii+MM] = DOT_JACOB(hh[ii][1], hh[jj][1], PP*WW);
	}
}
