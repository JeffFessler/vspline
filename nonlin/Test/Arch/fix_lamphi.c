/***************************** fix_lamphi.c ******************************
*	UNWRAP ELLIPSE ORIENTATION
*	USES floor()
*************************************************************************/
#include <math.h>
#include "nonlin.h"

	void fix_lamphi();


/*
*	INPUT
*		pp[1..N][1..2]		lambda, phi
*	RETURNS
*/

void fix_lamphi(parset, MM, NN)
PARAM_SET parset;
int MM, NN;
{
	int nn;

	PARAM_TYPE phi = parset[1][5];
	PARAM_TYPE frac = phi - floor(phi);
	int pi = ((int) floor(phi)) % 180;
	phi = parset[1][5] = pi + frac;

return;
	for (nn=2; nn <= NN; nn++)
	{
		PARAM_TYPE new = parset[nn][5];
		frac = new - floor(new);
		pi = ((int) floor(new)) % 180;

		while (pi < phi)
			pi += 180;

		if (phi-(pi-180) < pi-phi)
			pi -= 180;
		phi = parset[nn][5] = pi + frac;
	}
}
