/***************************** alloc.c ***********************************
*	SUBROUTINES FOR ALLOCATING AND FREEING MATRICES AND VECTORS
*	Jeff Fessler
*	COPYRIGHT 1989,  THE BOARD OF TRUSTEES OF THE
*	LELAND STANFORD JUNIOR UNIVERSITY.  ALL RIGHTS RESERVED.
*************************************************************************/
#include "defs-vspline.h"
#include "defs_vmatrix.h"

// DOUBLE PRECISION SUBROUTINES

/*
* ALLOCATE SPACE FOR A UNIT OFFSET VECTOR: V[1..N]
*/
double *alloc_vector_d(cint nn)
{
	double *p;
	Mem0(p, nn)
	return p-1;
}


/*
* ALLOCATE SPACE FOR A [1..rows][1..cols] MATRIX
* SUCH THAT MATRIX ELEMENTS ARE CONSECUTIVE
*/
double **alloc_matrix_d(int rows, cint cols)
{
	double **p;
	Mem0(p, rows)
	p--;
	if ( !(p[1] = alloc_vector_d(rows*cols)) )
		return NULL;
	for (; rows > 1; rows--)
		p[rows] = (rows-1)*cols + p[1];
	return p;
}

// SINGLE PRECISION SUBROUTINES


/*
* ALLOCATE SPACE FOR A UNIT OFFSET VECTOR: V[1..N]
*/
float *alloc_vector_f(cint nn)
{
	float *p;
	Mem0(p, nn)
	return p-1;
}


/*
* ALLOCATE SPACE FOR A [1..rows][1..cols] MATRIX
* SUCH THAT MATRIX ELEMENTS ARE CONSECUTIVE
*/
float **alloc_matrix_f(int rows, cint cols)
{
	float **p;
	Mem0(p, rows)
	p--;
	if ( !(p[1] = alloc_vector_f(rows*cols)) )
		return 0;
	for (; rows > 1; rows--)
		p[rows] = (rows-1)*cols + p[1];
	return p;
}


/*
* ALLOCATE SPACE FOR A [1..rows][1..cols][1..last] MATRIX
* SUCH THAT MATRIX ELEMENTS ARE CONSECTIVE
*/
float ***alloc_matrix3_f(cint rows, cint cols, cint last)
{
	int r;
	float ***p;
	Mem0(p, rows)
	p--;
	Mem0(p[1], rows*cols)
	p[1]--;
	for (r=rows; r > 1; r--)
		p[r] = (r-1)*cols + p[1];

	Mem0(p[1][1], rows*cols*last)
	p[1][1]--;

	for (r=rows*cols; r > 1; r--)
		p[1][r] = (r-1) * last + p[1][1];

	return p;
}
