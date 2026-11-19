/***************************** debug.c ***********************************
* MATRIX AND VECTOR DEBUGGING AND IOAND IO  SUBROUTINES
* print_xxx	ascii to stdout
* scan_xxx	ascii from stdin
* USES:  sprintf(), fwrite(), fopen(), fclose(), fflush(), etc
* Jeff Fessler
* COPYRIGHT 1989, THE BOARD OF TRUSTEES OF THE
* LELAND STANFORD JUNIOR UNIVERSITY.  ALL RIGHTS RESERVED.
*************************************************************************/
#include "mexdep.h"


// DOUBLE PRECISION SUBROUTINES

/*
* vect[1..N]
*/
void print_vector_d(cdouble *vect, cint N)
{
	for (int n=1; n <= N; n++)
	{
		printf("%g ", *++vect);
		if (!(n % 10))	printf("\n");
	}
	if (N % 10)	printf("\n");
}


/*
* vect[1..N] from stdin
*/
sof scan_vector_d(double *vect, cint N)
{
	for (int n=1; n <= N; n++)
		if (1 != scanf("%lg", &vect[n]))
			Fail("scanf")
	Ok
}


// SINGLE PRECISION SUBROUTINES

/*
* vect[1..N]
*/
void print_vector_f(cfloat *vect, cint N)
{
	for (int n=1; n <= N; n++)
	{
		printf("%g ", *++vect);
		if (!(n % 10))	printf("\n");
	}
	if (N % 10)	printf("\n");
}

/*
* vect[1..N] from stdin
*/
sof scan_vector_f(float *vect, cint N)
{
	for (int n=1; n <= N; n++)
		if (1 != scanf("%g", &vect[n]))
			Fail("scanf")
	Ok
}
