/*********************** mle() *******************************************
*	Independent maximum likelihood estimates for several measurements
*	Note the input formats!
*************************************************************************/
#include "nonlin.h"

	void mle_ind();

	extern void nl_alloc(), nl_setup(), nl_init(), nl_cleanup();
	extern sof nl_iter();
	extern Nonlin *nl_info();

#define alpha		((ALPHA *) 0)
#define ti		((VECT) 0)
#define estfile		((char *) 0)
#define FLAG_MLE	True

/*
*	INPUT:
*		meas[1..N][1..L]	original noisy measurements
*		istates[1..M*N]		initial estimates
*		prob			problem dependent structure
*		L			length of each measurement vector
*		M			# of unknowns
*		N			# of measurement sets
*		P			# of offsets
*		W			maximum signal width
*	EFFECT:
*		ostates[1..M*N]		final estimates
*/

void mle_ind(ostates, meas, istates, prob, LL, MM, nnmeas, PP, WW)
PARAM ostates, istates;
MEAS meas;
Problem *prob;
int LL, MM, nnmeas, PP, WW;
{
	int nn;
	Nonlin *init;
	truf flag_mod = prob->modify != 0;
	PARAM par;

	nl_alloc(LL, MM, 1, PP, WW, FLAG_MLE, flag_mod);
	nl_setup(ti, 1);
	init = nl_info();

	for (nn=0; nn < nnmeas; nn++) {
		par = istates + nn*MM;
		nl_init(init, &par-1, meas+nn, prob);
		init->bad.roughness = 0.;
		if (!nl_iter(init, alpha, meas+nn, prob, estfile))
			message("WARNING: unsuccessful iterations", (double) nn);
		Bcopy(&init->par[1][1], &ostates[nn*MM+1], MM);
	}

	nl_cleanup();
	free_info(init, 1);
}
