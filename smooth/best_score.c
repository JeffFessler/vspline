/*********************** best_score.c ************************************
* SMOOTH VECTOR MEASUREMENTS WITH AUTOMATIC SCORE CRITERIA
* Jeff Fessler
* COPYRIGHT 1989, THE BOARD OF TRUSTEES OF THE
* LELAND STANFORD JUNIOR UNIVERSITY. ALL RIGHTS RESERVED.
*************************************************************************/
#include "smooth.h"

	static int MM, nn, score_id;
	static int Flag_print;


/*
*	RETURN PARTICULAR SCORE
*/
static OPT_FUNC score_func(ALPHA_SET alpha)
{
	double ur, cv, gcv;
	double score;

	if ( !vs_setup_alpha(alpha) )			ABORT("8.")
	if ( !vs_smooth((VECT) 0, 0))			ABORT("9.")
	if ( !cross_val(&ur, &cv, &gcv, score_id) )	ABORT("10.")

	if (score_id == JOBUR)			score = ur;
	else if (score_id == JOBCV)		score = cv;
	else if (score_id == JOBGCV)		score = gcv;
	else					ABORT("Invalid score id.")

	if (Flag_print) {
#if 0
		message0("cv ", cv);
#else
		printf("cv %g ", cv);
#endif

		PRINT_ALPHA(alpha, MM);
	}

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
double best_score(CMAT ys, ALPHA_SET alpha, VECT ti, CMAT yi, COVAR *covar,
ctruf flag_iid, cint M, cint N, cint l_score_id, cdouble tol, cint chat)
{
	OPT_ARG score;
	DIRSET dirset;
	int iter;

	// "SAVE" THE PASSED VARIABLES
	score_id = l_score_id;	MM = M;	nn = N;
	Flag_print = chat;

	// SETUP: INDEPENDENT OF ALPHA
	if (!vs_setup_alloc(MM, nn))		ABORT("1.")
	if (!vs_setup_ti(ti, nn))		ABORT("2.")
	if (!vs_setup_yi(yi, MM, nn))		ABORT("3.")
	if (!vs_setup_cov(covar, flag_iid))	ABORT("4.")
	if (!cv_setup_alloc())			ABORT("5.")
	if (!cv_setup_covar())			ABORT("6.")

	// THE ROUTINES BELOW ALL DEPEND ON ALPHA

	dirset = SETUP_DIRSET(MM);
	powell(alpha, dirset, MM, (OPT_ARG) tol, &iter, &score, score_func);

	if ( !vs_setup_alpha(alpha) )		ABORT("6.2")
	if ( !vs_smooth(ys, 1))			ABORT("6.5")

	vs_cleanup();
	cv_cleanup();

	return (double) score;
}
