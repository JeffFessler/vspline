/************************* vs_smooth.c ***********************************
*	MAIN SUBROUTINES FOR
*	FIXED INTERVAL SMOOTHING OF VECTOR MEASUREMENTS WITH CUBIC SPLINES.
*
*	REFERENCE: Non-parametric Fixed-Interval Smoothing with Vector Splines
*			Jeff Fessler,  IEEE Trans-ASSP 1989
*	PROBLEM:
*		minimize (Y - a)' inv(Sigma) (Y - a)
*			+ \sum_{m=1}^M alpha_m c' R c
*		subject to (Q' @ I) a = (R @ I) c
*		(Sigma = blockdiag(sigma_1, ..., sigma_N))
*	SOLUTION:
*		c_alpha = inv[(R@inv(alpha)) + (Q'@I) Sigma (Q@I)] (Q'@I) Y
*		ahat = Y - Sigma (Q@I) c_alpha
*	("@" IS USED TO DENOTE MATRIX TENSOR (KRONECKER) PRODUCT.)
*	TO USE:
*		Initialize:
*		vs_setup_alloc(), vs_setup_ti(), vs_setup_yi(), vs_setup_cov()
*		Then for each value of alpha:
*		vs_setup_alpha(), vs_smooth()
*		Finally:
*		vs_cleanup()
*	There are essentially two versions contained within, depending
*	on whether the noise is IID or not.  If IID, then only one covariance
*	matrix needs to be passed, ie:
*		Flag_iid:
*			1: noise is iid: Sigma[1..MxM]	covar[1]
*			0: noise varies: Sigma[1..MMxN]	covar[1..N]
*	USES:	exp()
*	Note that the actual value of alpha is: exp(alpha) * (N-1)/(tmax-tmin)
*	Jeff Fessler
*	COPYRIGHT 1989,  THE BOARD OF TRUSTEES OF THE
*	LELAND STANFORD JUNIOR UNIVERSITY.  ALL RIGHTS RESERVED.
*************************************************************************/
#include "smooth.h"

	int NI=0, MI=0;		/* 'I' SUFFIX: SIZES WHEN INITIALIZED	*/
	int NG, MG;		/* 'G' SUFFIX: MEANS "GLOBAL"		*/
	truf Flag_iid;
	truf Flag_equispace;	/* EQUALLY SPACED SAMPLE POINTS		*/
	static double ti_range;		/* ti[N] - ti[1]		*/

	/*
	*	GLOBAL VARIABLES (used in cv_score.c)
	*/
	BAND Qt;			/* SEE vs_init.c FOR FORMAT	*/
	static SYMBAND R;		/* R(1,N-2)			*/
	static SYMBAND QtSQ;		/* QtSQ(3M-1,M(N-2))		*/
	static VECT alpha_inv;		/* [1..MG]	diag(1/alpha)	*/
	static SYMBAND R_A;		/* R_A(2M-1,M(N-2))		*/
	BAND BL;			/* BL(3M-1,M(N-2))		*/
	VECT BD;			/* BD[1..M(N-2)]		*/
	static VECT c_alpha;		/* c_alpha[1..M(N-2)]		*/
	OMATRIX sig_q_c;		/* sig_q_c[1..N][1..M]		*/
	static VECT QtY;		/*, QtY[1..M(N-2)]		*/

	/*
	*	GLOBAL POINTERS
	*/
	CMAT Gyi;
	COVAR *Gcovar;


/*
*	The next few routines are indepentdent of alpha, so should
*	only be called once prior to searching for optimal alpha
*/

/*
*	ALLOCATE MATRICES, USING UPPER BOUNDS FOR THE SIZES TO
*	BE USED IN SUBSEQUENT CALLS.
*	INPUT:
*		NI:			# of sample points
*		MI:			length of sample vectors
*	ALLOCATED:
*		R, Qt, QtSQ, QtY
*		alpha_inv, R_A, BL, BD, c_alpha, sig_q_c
*	RETURN:
*		1 if OK, 0 if too big to initialize.
*/
sof vs_setup_alloc(int M, int N)
{
	NI = N;
	MI = M;

	if (NI < 3)
		Fail("need NI >= 3")

	/* fix:	if failure, free the others */
	Call(Qt		= ALLOC_BAND, (NI-2, 2))
	Call(QtSQ	= ALLOC_BAND, (MI*(NI-2), 3*MI-1))
	Call(R		= ALLOC_BAND, (NI-2, 1))
	Call(R_A	= ALLOC_BAND, (MI*(NI-2), 2*MI-1))
	Call(BL		= ALLOC_BAND, (MI*(NI-2), 3*MI-1))
	Call(BD		= ALLOC_VECTOR, (MI*(NI-2)))
	Call(QtY	= ALLOC_VECTOR, (MI*(NI-2)))
	Call(c_alpha	= ALLOC_VECTOR, (MI*(NI-2)))
	Call(alpha_inv	= ALLOC_VECTOR, (MI))
	Call(sig_q_c	= ALLOC_MATRIX, (NI, MI))

	Ok
}


/*
*	Setup given sample points
*
*	INPUT:
*		ti[1..N]		sample points (or null if equispaced)
*		NG:			# of sample points
*	RETURN:				proper ordering?
*/
sof vs_setup_ti(VECT ti, int N)
{
	NG = N;
	if (NG > NI)	Fail("NG > NI")
	if (NG < 5)	Fail("NG < 5")	/* MUST HAVE AT LEAST 3 POINTS	*/

	Flag_equispace = (ti == (VECT) 0);
	ti_range = (Flag_equispace) ? (NG-1) : (ti[NG]-ti[1]);

	Call(spline_mat, (R, Qt, ti, NG))	/* NOT ORDERED	*/
	Ok
}


/*
*	Setup given measurement data.
*	INPUT:
*		yi[1..MN]		measurements
*		N:			# of sample points
*		M:			length of sample vectors
*	RETURN:				ok dimensions?
*/
sof vs_setup_yi(VECT yi, int M, int N)
{
	MG = M;
	if (MG != MI)	Fail("MG != MI")
	if (NG != N)	Fail("NG != N")

	Gyi = yi;
	form_qt_Y(QtY, Qt, yi, MG, NG);
	Ok
}


/*
*	Setup given error covariances
*	INPUT:
*		covar[1..?]		noise covariances (see above)
*		flag_iid:		see above
*	USES:
*		Qt
*	EFFECT:
*		QtSQ, Gcovar, Flag_iid
*/
sof vs_setup_cov(COVAR *covar, ctruf flag_iid)
{
	Flag_iid = flag_iid;
	Gcovar = covar;
	ZERO_BAND(QtSQ, 3*MI-1, MI*(NI-2));
	form_qt_s_q(QtSQ, Qt, covar, MG, NG, Flag_iid);
	Ok
}


/*
*	FREE ALL THE ALLOCATED MATRICES
*/
void vs_cleanup()
{
	FREE_BAND(Qt);
	FREE_BAND(QtSQ);
	FREE_BAND(R);
	FREE_BAND(R_A);
	FREE_BAND(BL);
	FREE_VECTOR(BD);
	FREE_VECTOR(QtY);
	FREE_VECTOR(c_alpha);
	FREE_VECTOR(alpha_inv);
	FREE_MATRIX(sig_q_c);
	MI = NI = 0;
}


/*
*	Now, the following routines DO depend on alpha
*/

/*
*	Finally, setup given the smoothing parameter alpha.
*	INPUT:
*		alpha[1..M]		smoothing parameters
*/
sof vs_setup_alpha(ALPHA_SET alpha)
{
	int m;
	/* double alpha_fac = (double) NG;	OLD WAY	*/
	double alpha_fac = (double) (NG-1.) / ti_range;

	/*
	*	Compute diagonal of matrix with 1/Alpha on diagonal
	*/
	for (m=MG; m >= 1; m--)
		alpha_inv[m] = exp(-alpha[m]) / alpha_fac;

	ZERO_BAND(R_A, 2*MI-1, MI*(NI-2));
	kron_diag(R_A, R, 1, NG-2, alpha_inv, MG);
	Ok
}


/*
*	ALGORITHM FOR SMOOTHING VECTOR MEASUREMENTS
*	job OPTIONS:
*		0:	ignore result, just compute sig_q_c (for cross_val())
*		1:	result is yi
*		2:	result is (yi - yhat) on output
*		4:	result is yi on input -> yhat on output
*	INPUTS:
*		result[1..MxN]		see above
*		job
*	OUTPUT:
*		result[1..MxN]		see above
*	USES:
*		everything
*	EFFECT:
*		sig_q_c:		set to (yi - yhat)
*/
sof vs_smooth(CMAT result, int job)
{
	VECT q_c;
	TMP_ALLOC_VECT(q_c,MG*NG)

	Call(vs_transfer, ())
	solve_sym_band(c_alpha, BL, BD, QtY, 3*MG-1, MG*(NG-2), 0);

	form_q_c(q_c, Qt, c_alpha, MG, NG);
	form_sig_q_c(sig_q_c, Gcovar, q_c, MG, NG, Flag_iid);
	if (job == 1)
		SUB_VECTOR(result, Gyi, sig_q_c[1], MG*NG);
	if (job == 2)
		COPY_VECTOR(sig_q_c[1], result, MG*NG);
	if (job == 4)
		DEC_VECTOR(result, sig_q_c[1], MG*NG);

	TMP_FREE_VECTOR(q_c);
	Ok
}


/*
*	SETUP MATRIX FOR COMPUTING c_alpha
*	USES:
*		R_A, QtSQ
*	EFFECT:
*		BL, BD
*/
sof vs_transfer()
{
	SYMBAND sum;
	TMP_ALLOC_BAND(sum,MG*(NG-2),3*MG-1)

	form_band(sum, R_A, QtSQ, MG, NG);
	ZERO_BAND(BL, 3*MI-1, MI*(NI-2));
	Call(chol_sym_band, (BL, BD, sum, 3*MG-1, MG*(NG-2)))

	TMP_FREE_BAND(sum);
	Ok
}


/*
*	The routines below actually compute the rougness penalty:
*	alpha \int g''(t)^2 dt
*	They differ in their prerequisites and computation counts.
*/

/*
*	COMPUTE ROUGHNESS PENALTY
*	vs_smooth() MUST BE CALLED FIRST.
*	USES:
*		c_alpha AND alpha_inv
*	RETURNS:
*		rough			roughness term c' (R @ a) c
*/
double roughness()
{
	int m, n;
	double rough = 0.;
	if (!Flag_equispace)
		return post_penalty();

	for (m=MG; m >= 1; m--)
	{
		double ssq = 2. * c_alpha[(NG-3)*MG+m] * c_alpha[(NG-3)*MG+m];
		for (n=NG-3; n >= 1; n--)
		{
			ssq += c_alpha[(n-1)*MG+m] *
			(2. * c_alpha[(n-1)*MG+m] + c_alpha[(n  )*MG+m]);
		}
		rough += ssq * alpha_inv[m];
	}

	return rough * 1./3.;
}


/*
*	COMPUTE ROUGHNESS PENALTY FOR SPLINE COEFFICIENTS
*	IMMEDIATELY AFTER SMOOTHING
*	RETURNS:	\sum_m alpha_m c' R c = c_alpha' R_A c_alpha
*	USES:	c_alpha, R_A, MG, NG
*/
double post_penalty()
{
	double pen;
	VECT r_c;
	TMP_ALLOC_VECT(r_c,MG*(NG-2))

	mult_band_vect(r_c, R_A, c_alpha, 2*MG-1, MG*(NG-2));
	pen = DOT_VECTOR(r_c, c_alpha, MG*(NG-2));

	TMP_FREE_VECTOR(r_c);
	return pen;
}


/*
*	COMPUTE ROUGHNESS PENALTY FOR SPLINE COEFFICIENTS
*	vs_setup_alloc() and vs_setup_ti() MUST BE CALLED FIRST.
*	INPUT:
*		alpha[1..M]		smoothing parameters
*		yi[1..MxN]		points to be interpolated
*	OUTPUT:
*		pen			penalty term: c' (R @ a) c
*					= \sum_{m=1}^M alpha_m c_m' R c_m
*	RETURNS:			computed?
*	USES:	QtY, c_alpha, alpha_inv, R_A, BL, BD, MG, NG
*/
sof penalty(double *pen, ALPHA_SET alpha, VECT yi, int M, int N)
{
	VECT RAD, res;
	BAND RAL;

	Call(vs_setup_yi, (yi, M, N))
	Call(vs_setup_alpha, (alpha))

	TMP_ALLOC_VECT(RAD,MG*(NG-2))
	TMP_ALLOC_VECT(res,MG*(NG-2))
	TMP_ALLOC_BAND(RAL,MG*(NG-2),2*MG-1)

	Call(chol_sym_band, (RAL, RAD, R_A, 2*MG-1, MG*(NG-2)))
	solve_sym_band(res, RAL, RAD, QtY, 2*MG-1, MG*(NG-2), 0);
	*pen = DOT_VECTOR(res, QtY, MG*(NG-2));

	TMP_FREE_VECTOR(RAD);
	TMP_FREE_VECTOR(res);
	TMP_FREE_BAND(RAL);
	Ok
}
