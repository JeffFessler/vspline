/***************************** pseudo.c **********************************
*	Routines for forming pseudo measurements for multiple object recon.
*************************************************************************/
#include <stdio.h>
#include "nonlin.h"

	void meas_dec(), meas_inc();


/*
*	DECREMENT MEASUREMENTS BY ONE SIGNAL SET
*/

void meas_dec(meas, sig, off, nn, PP)
register MEAS meas;
register SIGNAL_SET sig;
register OFFSET_SET off;
int nn, PP;
{
	register int jj, pp, ww;

	for (jj=nn; jj >= 1; jj--)
	{
		for (pp=PP; pp >= 1; pp--)
		{
			MEAS_TYPE *pm = &meas[jj][off[jj].start[pp]];
			SIGNAL_TYPE *po = &sig[jj][pp][1];
			for (ww=off[jj].count[pp]; ww; ww--)
				*pm++ -= *po++;
		}
	}
}


/*
*	INCREMENT MEASURMENT BY ONE SIGNAL SET
*/

void meas_inc(meas, sig, off, nn, PP)
register MEAS meas;
register SIGNAL_SET sig;
register OFFSET_SET off;
int nn, PP;
{
	register int jj, pp, ww;

	for (jj=nn; jj >= 1; jj--)
	{
		for (pp=PP; pp >= 1; pp--)
		{
			MEAS_TYPE *pm = &meas[jj][off[jj].start[pp]];
			SIGNAL_TYPE *po = &sig[jj][pp][1];
			for (ww=off[jj].count[pp]; ww >= 1; ww--)
				*pm++ += *po++;
		}
	}
}
