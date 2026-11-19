/************************ sigjac_edge() **********************************
* COMPUTE SIGNALS AND JACOBIANS
* FOR EDGE - No offsets.
*
* PARAMETERS: (M = 1) (position of edge)
* USES:	ceil()
*************************************************************************/
#include <math.h>
#include "nonlin.h"
#include "sigarg.h"

//	void all_sigjac_edge();
//	sof sigjac_edge();


/*
* in
*	par: [1..M]		parameters
*	arg:			simply L, the image width
* out
*	sig: [1..1][1..L]	Computed signal
*	jac: [1..1][1..1][1..L]	Jacobians
*	off: [1..1]		1
*/

static sof sigjac_edge(
SIGNAL sig,
JACOBIAN jac,
OFFSET *off,
VECT par,
Void *arg)
{
	int ii;
	ARGSTRUCT *as = (ARGSTRUCT *) arg;
	int LL = as->LL;
	int WW = as->WW;
	float *Jedge = &jac[1][1][1];
	float *Sig = &sig[1][1];
	float edge = par[1];

	off->start[1] = 1;
	off->count[1] = 1 + ceil(edge);
	if (off->count[1] > LL || off->count[1] > WW)
		return FAILURE;

	for (ii=1; ii <= edge; ii++, Sig++, Jedge++)
	{
		*Jedge = 0.0;
		*Sig = 1.0;
	}

	// now do final point
	if (ii - 1 <= edge)
	{
		*Jedge = 1.0;
		*Sig = edge - (ii-1);
	}

	Ok
}


/*
* Signals and jacobians for entire data set
*/

void all_sigjac_edge(out, sig, jac, off, par, sigarg, MM, NN, OO, WW)
SIGJAC_RETURN *out;
SIGNAL_SET sig;
JACOB_SET jac;
OFFSET_SET off;
PARAM_SET par;
Void *sigarg;
int MM, NN, OO, WW;
{
	int nn;

	if (MM != 1)	ABORT("Wrong # of parameters!")

	*out = 0;
        ZERO_SIGNAL_SET(sig, NN, OO, WW);

    for (nn = NN; nn >= 1; nn--)
    {
        ZERO_JACOB(jac[nn], MM, OO, WW);
	if (!sigjac_edge(sig[nn], jac[nn], &off[nn], par[nn], sigarg))
	{
		*out = nn;
		return;
	}
    }
}
