/***************************** alloc.c ***********************************
* Allocation subroutines for nonlinear estimation.
*************************************************************************/
#include "nonlin.h"

#if 0
	JACOB_SET alloc_jacob_set();
	void free_jacob_set(), copy_jacob_set();
	OFFSET_SET alloc_offset_set();
	void free_offset_set(), copy_offset_set();
#endif


/*
* Allocate N JACOBIANS
* RETURNS: jac[1..N][1..M][1..P][1..W]
*/

JACOB_SET alloc_jacob_set(int nn, int MM, int PP, int WW)
{
	JACOB_SET jac;
	int n;

	Alloc0(jac, JACOBIAN, nn, sizeof(JACOBIAN))
	jac--;

	for (n=nn; n >= 1; n--)
	    if (!(jac[n] = ALLOC_JACOB(MM, PP, WW)))
		return (JACOB_SET) 0;

	return jac;
}


/*
*
*/

void free_jacob_set(JACOB_SET jac, int nn)
{
	int n;
	for (n=nn; n >= 1; n--)
	    FREE_MATRIX3(jac[n]);
	Free(jac + 1);
}


/*
* Copy a JACOB_SET
* src,dst[1..N][1..M][1..P][1..W]
*/

void copy_jacob_set(JACOB_SET src, JACOB_SET dst, int nn, int MM, int PP, int WW)
{
	int n;

    for (n=nn; n >= 1; n--)
	COPY_MATRIX3(src[n], dst[n], MM, PP, WW);
}


/*
* RETURN [1..N] offset structures,
* each with pointers for [1..P] pulses.
*/

OFFSET_SET alloc_offset_set(int nn, int PP)
{
	int n;
	OFFSET_SET off;
	Alloc0(off, OFFSET, nn, sizeof(OFFSET))
	off--;
	Alloc0(off[1].start, short, nn*PP, sizeof(*off[1].start))
	Alloc0(off[1].count, short, nn*PP, sizeof(*off[1].count))
	off[1].start--;
	off[1].count--;
	for (n=nn; n > 1; n--)
	{
		off[n].start = off[1].start + (n-1) * PP;
		off[n].count = off[1].count + (n-1) * PP;
	}

	return off;
}


/*
*/

void free_offset_set(OFFSET_SET off)
{
	Free(&off[1].start[1]);
	Free(&off[1].count[1]);
	Free(&off[1]);
}


/*
* Copy a OFFSET_SET
*/

void copy_offset_set(OFFSET_SET src, OFFSET_SET dst, int nn, int PP)
{
	Bcopy(&src[1].count[1], &dst[1].count[1], nn*PP);
	Bcopy(&src[1].start[1], &dst[1].start[1], nn*PP);
}
