// vector_d.c
// UTILITY VECTOR SUBROUTINES
// Jeff Fessler
// COPYRIGHT 1989, THE BOARD OF TRUSTEES OF THE
// LELAND STANFORD JUNIOR UNIVERSITY. ALL RIGHTS RESERVED.

#include "defs_vmatrix.h"

void vector_add_d(Vd out, cVd in1, cVd in2, int N)
{
	for (; N; N--)
		*(++out) = *(++in1) + *(++in2);
}

void vector_sub_d(Vd out, cVd in1, cVd in2, int N)
{
	for (; N; N--)
		*(++out) = *(++in1) - *(++in2);
}

void vector_inc_d(Vd in1, cVd in2, int N)
{
	for (; N; N--)
		*(++in1) += *(++in2);
}

void vector_dec_d(Vd in1, cVd in2, int N)
{
	for (; N; N--)
		*(++in1) -= *(++in2);
}

void scale_vector_d(Vd v, cdouble s, int N)
{
	for ( ; N; N--)
		*++v *= s;
}

void vector_f_to_d(Vd d, cVf f, int n)
{
	for (; n; n--)
		*++d = *++f;
}

void inc_scaled_vector_d(Vd v1, cdouble s, cVd v2, int n)
{
	for (; n; n--)
		*++v1 += s * *++v2;
}

double norm_vector_d(cVd v, int n)
{
	register double d = 0.;
	v++;
	for (; n; n--, v++)
		d += *v * *v;
	return d;
}

double dot_vector_d(cVd v1, cVd v2, int n)
{
	register double d = 0.;
	for (; n; n--)
		d += *++v1 * *++v2;
	return d;
}

// find min and max of a vector v1[1..n]
void min_max_vector_d(double *vmin, double *vmax, cVd v, int n)
{
	*vmin = *vmax = v[1];

	for (; n; n--) {
		++v;
		if (*v > *vmax)	*vmax = *v;
		if (*v < *vmin)	*vmin = *v;
	}
}
