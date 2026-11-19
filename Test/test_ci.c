/*********************** test_ci.c ***************************************
* Test cubic_interp().
* Read data, interpolate, print out.
* Jeff Fessler
* COPYRIGHT 1989, THE BOARD OF TRUSTEES OF THE
* LELAND STANFORD JUNIOR UNIVERSITY.  ALL RIGHTS RESERVED.
*************************************************************************/
#include "def,test.h"

#define Usage	{ fprintf(stderr, "Usage: %s Ni Nx\n", argv[0]); exit(-1);}

int main(int argc, char **argv)
{
	VECT ti, yi, xs, ys;
	int Ni, Nx;

	if (argc != 3)
		Usage

	Ni = atoi(argv[1]);
	Nx = atoi(argv[2]);

	TMP_ALLOC_VECT(ti, Ni);
	TMP_ALLOC_VECT(yi, Ni);
	TMP_ALLOC_VECT(xs, Nx);
	TMP_ALLOC_VECT(ys, Nx);

	SCAN_VECTOR(ti, Ni);
	SCAN_VECTOR(yi, Ni);
	SCAN_VECTOR(xs, Nx);

	cubic_interp(ys, ti, yi, xs, Ni, Nx);

	PRINT_VECTOR(ys, Nx);

	return 0;
}
