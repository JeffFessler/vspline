/*********************** cv_score.c **************************************
*	SUBROUTINES FOR COMPUTING CROSS-VAILDATION LIKE SCORES FOR
*	FIXED-INTERVAL SMOOTHING OF VECTOR MEASUREMENTS WITH CUBIC SPLINES.
*	TO BE USED ONLY AFTER SMOOTHING WITH THE SUBROUTINES IN vs_smooth.c.
*	Jeff Fessler
*	COPYRIGHT 1989,  THE BOARD OF TRUSTEES OF THE
*	LELAND STANFORD JUNIOR UNIVERSITY.  ALL RIGHTS RESERVED.
*************************************************************************/
#include "smooth.h"

	/* VARIABLES DECLARED, ALLOCATED, AND PRECOMPUTED IN vs_smooth.c */
	extern int NG, MG;		/* 'G' SUFFIX FOR "GLOBAL"	*/
	extern int NI, MI;		/* 'I' SUFFIX FOR "INITIALIZED"	*/
	extern truf Flag_iid;
	extern BAND Qt, BL;
	extern VECT BD;
	extern OMATRIX sig_q_c;
#define DELTA	sig_q_c			/* DELTA[1..N][1..M] = (yi-yhat)*/

	/* POINTERS INITIALIZED IN vs_smooth.c				*/
	extern COVAR *Gcovar;

	/* GLOBAL VARIABLES FOR THIS FILE				*/
	static SYMBAND BI;		/* BI(3M-1,M(N-2))		*/
	static SYMBAND Fn, FL;		/* F, FL: SYMBAND(M-1,M)	*/
	static VECT FD;			/* FD[1..N]			*/
	static double trace_pi;		/* trace(Pi)			*/


/*
* ALLOCATE MATRICES
*/
sof cv_setup_alloc(void)
{
	Call(BI	= ALLOC_BAND,	(MI*(NI-2), 3*MI-1))
	Call(Fn	= ALLOC_MATRIX,	(MI, MI))
	Call(FL	= ALLOC_MATRIX,	(MI, MI))
	Call(FD	= ALLOC_VECTOR,	(MI))

	Ok
}


/*
*	PRECOMPUTE TERMS USED FOR CROSS VALIDATION AND OTHER SCORES
*	USES:
*		Gcovar
*/
sof cv_setup_covar(void)
{
	trace_pi = 0.;
	for (int n = Flag_iid ? 1 : NG; n >= 1; n--)
	{
		COVAR *covar = &Gcovar[n];
		covar->trace = band_trace(covar->cov, MG-1, MG);
		mult_full_band(covar->cov2, covar->cov, covar->cov, MG);
		trace_pi += covar->trace;
	}

	if (Flag_iid)
		trace_pi *= NG;

	Ok
}


/*
*	FREE ALL THE ALLOCATED MATRICES
*/
void cv_cleanup(void)
{
	FREE_BAND(BI);
	FREE_MATRIX(Fn);
	FREE_MATRIX(FL);
	FREE_VECTOR(FD);
}


/*
*	COMPUTE CROSS VALIDATION SCORES AND OTHER STATISTICS
*
*	CV = 1/N \sum_{n=1}^N | Sig_n^{-1/2} inv(I-A) (Y - Yhat)
*	UR = (1/N) [ | Y - Yhat |^2 - 2 trace(Sig (I-A)) + trace(Sig) ]
*	GCV = (1/N) RSS / [ (1/N) trace(I - A) ]^2
*
*	I - A = Pi F;  F = (Q @ I) BI (Q' @ I); BI = inv(BL BD BL')
*	I - A_n = Sig_n F_n
*	e_n = Y_n - Yhat_n
*	CV_n = 1/N \sum b_n Sig_n^{-1/2} b_n, b_n = inv(Sig_n F_n) e_n
*
*	INPUT:
*		job:			which scores to compute
*	OUTPUT:
*		ur:			(1) UNBIASED RISK SCORE
*		cv:			(2) CROSS VALIDATION SCORE
*		gcv:			(4) GENERALIZED CV SCORE
*	USES:
*		DELTA (sig_q_c) COMPUTED BY vect_smooth()
*		Gcovar			error covariance structure(s)
*/
sof cross_val(double *ur, double *cv, double *gcv, cint job)
{
	COVAR *covar = &Gcovar[1];
	double rss = 0., errn = 0.;

	// COMPUTE THE CENTRAL BANDS OF THE INVERSE MATRIX
	invert_sym_band(BI, BL, BD, 3*MG-1, MG*(NG-2));

	*cv = *ur = *gcv = 0.;

	for (int nn=1; nn <= NG; nn++)
	{
		VECT err = DELTA[nn];

		select_nth_qbiqt(Fn, BI, Qt, MG, NG, nn);
		Call(chol_sym_band, (FL, FD, Fn, MG-1, MG))

		/*
		*	Compute unbiased risk estimator
		*	Compute norm of residuals
		*/
		if (job & JOBUR) {
			errn += NORM_VECTOR(err, MG);
			*ur += ur_term(Fn, covar->cov2, MG);
		}

		/*
		*	Compute cv score
		*/
		if (job & JOBCV)
			*cv += cv_term(err, FL, FD, covar->inv, MG);

		/*
		*	Compute generalized cross validation score
		*	Compute residual sum of squares
		*/
		if (job & JOBGCV) {
			rss += rss_term(err, covar->inv, MG);
			*gcv += gcv_term(Fn, covar->cov, MG);
		}

		/*
		*	If not iid noise, increment to next covariance
		*/
		if (!Flag_iid)
			covar++;
	}

	/*
	*	Correct scores for constant factors (saves computations)
	*/
	if (job & JOBUR)
		*ur = (errn - 2. * *ur + trace_pi) / NG;
	if (job & JOBCV)
		*cv /= NG;
	if (job & JOBGCV)
		*gcv = NG * rss / (*gcv * *gcv);

	Ok
}


/*
*	COMPUTE TRACE (I - A)
*	INPUT:
*		alpha[1..M]		smoothing parameters
*	OUTPUT:
*		trace:			trace(I - A(alpha))
*	USES:
*		Gcovar
*/
sof trace_i_a(double *trace, ALPHA_SET alpha)
{
	COVAR *covar = &Gcovar[1];

	Call(vs_setup_alpha, (alpha))
	Call(vs_transfer, ())

	// COMPUTE THE CENTRAL BANDS OF THE INVERSE MATRIX
	invert_sym_band(BI, BL, BD, 3*MG-1, MG*(NG-2));

	*trace = 0.;
	for (int nn=1; nn <= NG; nn++)
	{
		select_nth_qbiqt(Fn, BI, Qt, MG, NG, nn);
		Call(chol_sym_band, (FL, FD, Fn, MG-1, MG))

		*trace += gcv_term(Fn, covar->cov, MG);

		// IF NOT IID NOISE, INCREMENT TO NEXT COVARIANCE
		if (!Flag_iid)
			covar++;
	}

	Ok
}
