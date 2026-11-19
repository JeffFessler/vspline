/************************** interp.c *************************************
* CUBIC-SPLINE INTERPOLATION
* Jeff Fessler
* COPYRIGHT 1989, THE BOARD OF TRUSTEES OF THE
* LELAND STANFORD JUNIOR UNIVERSITY. ALL RIGHTS RESERVED.
*************************************************************************/
#include "smooth.h"


/*
* return i such that list[i] <= x < list[i+1]
*/
static int dumb_search(cdouble x, cVECT list)
{
	int ii;
	for (ii=2; ;ii++)
		if (x < list[ii])
			return ii-1;
}


/*
* linear interpolation between (ti[1],yi[1]) and (ti[2],yi[2]).
*/
static void simple_interp(VECT ys, cVECT ti, cVECT yi, cVECT xs, cint nx)
{
	cdouble dx = ti[2] - ti[1];
	cdouble dy = yi[2] - yi[1];

	for (int ix=1; ix <= nx; ix++)
		ys[ix] = yi[1] + (xs[ix] - ti[1])*dy/dx;
}


/*
*	INTERPOLATION
*	INPUT:
*		ti[1..N]		knots
*		yi[1..N]		function value at knots
*		xs[1..nx]		new sample points
*	OUTPUT:
*		ys[1..nx]		interpolated values
*	RETURN:				1 if OK, 0 if error
*/
sof cubic_interp(VECT ys, cVECT ti, cVECT yi, cVECT xs, cint nn, cint nx)
{
	BAND Qt, R, RL;
	VECT RD, cc;

	if (nn == 2) {
		simple_interp(ys, ti, yi, xs, nx);
		Ok
	}

	TMP_ALLOC_BAND(Qt, nn-2, 2);
	TMP_ALLOC_BAND(R, nn-2, 1);
	TMP_ALLOC_BAND(RL, nn-2, 1);
	TMP_ALLOC_VECT(RD, nn-2);
	TMP_ALLOC_VECT(cc, nn-2);

	// COMPUTE 'C' VECTOR
	Call(spline_mat, (R, Qt, ti, nn))
	form_qt_Y(cc, Qt, yi, 1, nn);
	Call(chol_sym_band, (RL, RD, R, 1, nn-2))
	solve_sym_band((VECT) 0, RL, RD, cc, 1, nn-2, 1);

	for (int ix=1; ix <= nx; ix++)
	{
		double xx = xs[ix];

		// FIX: WRONG END CONDITIONS
		if (xx >= ti[nn])
		{
			double hn1 = ti[nn] - ti[nn-1];
			double bn1= (yi[nn]-yi[nn-1])/hn1 + hn1/6.*cc[nn-2];
			ys[ix] = yi[nn] + bn1 * (xx - ti[nn]);
		}
		else if (xx <= ti[1])
		{
			double h1 = (ti[2]-ti[1]);
			double b1 = (yi[2]-yi[1])/h1 - h1/6.*cc[1];
			ys[ix] = yi[1] + b1*(xx - ti[1]);
		}
		else
		{
			// 1 <= it <= nn-1
			cint it = dumb_search(xx, ti);
			cdouble dt = xx - ti[it];
			cdouble hi = ti[it+1] - ti[it];
			cdouble ai = yi[it];
			cdouble ci = (it > 1) ? cc[it-1] : 0;
			cdouble ci1 = (it+1 < nn) ? cc[it] : 0;
			cdouble di = (ci1 - ci) / hi;
			cdouble bi = (yi[it+1]-ai)/hi - (ci1 + 2*ci) * hi / 6.;
			ys[ix] = ai + (bi + (ci/2. + di/6.*dt)*dt)*dt;
		}
	}

	TMP_FREE_BAND(Qt);
	TMP_FREE_BAND(R);
	TMP_FREE_BAND(RL);
	TMP_FREE_VECTOR(RD);
	TMP_FREE_VECTOR(cc);

	Ok
}
