/************************ sigjac_taurad() ********************************
* COMPUTE SIGNALS AND JACOBIANS OF PROJECTION OF A CIRCLE
* WITH UNIT DENSITY AND VARIABLE POSITION AND RADIUS
*
* PARAMETERS: (M = 2) (tau,rad)
* USES: fabs(), sqrt(), asin(), ceil(), floor();
*************************************************************************/
#include <math.h>
#include "nonlin.h"
#include "sigarg.h"


/*
* in
*	par: [1..M]		parameters
*	arg:			L, the image width
* out
*	sig: [1..1][1..W]	Computed signal
*	jac: [1..M][1..1][1..W]	Jacobians
*/

static sof sigjac_taurad(
SIGNAL sig,
JACOBIAN jac;
OFFSET *off,
VECT par,
Void *arg)
{
	int ii;
	ARGSTRUCT *as	= (ARGSTRUCT *) arg;
	int LL	= as->LL;
	int WW	= as->WW;
	float *Jtau	= &jac[1][1][1]; // Jacobians
	float *Jrad	= &jac[2][1][1];
	float *thesig	= &sig[1][1];

	double tau	= par[1];
	double rad	= par[2];
	double gain	= rad * rad;

	off->start[1]	= 1 + floor(tau - rad);
	off->count[1]	= ceil(tau + rad) - off->start[1] + 1;
	if (off->count[1] > WW || off->start[1] < 1 || off->start[1]+off->count[1]-1 > LL || rad <= 0.)
		return FAILURE;

	tau = tau - off->start[1] + 1;

	/* LOOP THROUGH THE NON-ZERO SIGNAL PIXELS			*/
	for (ii=1; ii <= off->count[1]; ii++, thesig++, Jtau++, Jrad++)
	{
		double d = ii - tau;

		double s0, s1, as1, as0;
		register double t0 = (d-1.) / rad;
		register double t1 = d / rad;

		if (t1 >=  1.)	t1 =  1.;
		if (t0 <= -1.)	t0 = -1.;
		s1	= sqrt(1 - t1*t1);
		s0	= sqrt(1 - t0*t0);
		as1	= asin(t1);
		as0	= asin(t0);

		*Jtau	= rad * 2. * (s0 - s1);
		*Jrad	= rad * 2. * (as1 - as0);
		*thesig	= gain * ( (t1 * s1 + as1) - (t0 * s0 + as0) );
	}

	Ok
}


/*
* SIGNALS AND JACOBIANS FOR ENTIRE DATA SET
*/

void all_sigjac_taurad(
SIGJAC_RETURN *out,
SIGNAL_SET sig,
JACOB_SET jac,
OFFSET_SET off,
PARAM_SET par,
Void *sigarg,
int MM, int NN, int PP, int WW)
{
	int nn;

	if (MM != 2)		ABORT("Wrong # of parameters!")

	*out = 0;
	ZERO_SIGNAL_SET(sig, NN, PP, WW);

	for (nn = NN; nn >= 1; nn--)
	{
		ZERO_JACOB(jac[nn], MM, PP, WW);
		if (!sigjac_taurad(sig[nn], jac[nn], &off[nn], par[nn], sigarg))
		{
			*out = nn;
			return;
		}
	}
}
