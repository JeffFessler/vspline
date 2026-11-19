/***************************** ident.c ***************************
* Return identity matrix.
* Jeff Fessler
* COPYRIGHT 1989,  THE BOARD OF TRUSTEES OF THE
* LELAND STANFORD JUNIOR UNIVERSITY.  ALL RIGHTS RESERVED.
*****************************************************************/
#include "defs_vmatrix.h"


/*
* IDENTITY MATRIX
* OUTPUT:
*	ident[1..M][1..M]		M*M identity matrix
*/
float **identity_f(int M)
{
	float **ident;
	if (!(ident = ALLOC_MATRIX_F(M, M)))	ABORT("mem alloc")

	for (; M; M--)
		ident[M][M] = 1.;

	return ident;
}


/*
* IDENTITY MATRIX
* OUTPUT:
*	ident[1..M][1..M]		M*M identity matrix
*/
double **identity_d(int M)
{
	double **ident;
	if (!(ident = ALLOC_MATRIX_D(M, M)))	ABORT("mem alloc")

	for (; M; M--)
		ident[M][M] = 1.;

	return ident;
}
