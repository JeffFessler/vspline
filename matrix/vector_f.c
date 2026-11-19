// vector_f.c
// Precision UTILITY VECTOR SUBROUTINES
// Jeff Fessler
// COPYRIGHT 1989, THE BOARD OF TRUSTEES OF THE
// LELAND STANFORD JUNIOR UNIVERSITY. ALL RIGHTS RESERVED.

#include "defs_vmatrix.h"

void vector_add_f(Vf out, cVf in1, cVf in2, int N)
{
	for (; N; N--)
		*(++out) = *(++in1) + *(++in2);
}

void vector_sub_f(Vf out, cVf in1, cVf in2, int N)
{
	for (; N; N--)
		*(++out) = *(++in1) - *(++in2);
}

void vector_inc_f(Vf in1, cVf in2, int N)
{
	for (; N; N--)
		*(++in1) += *(++in2);
}

void vector_dec_f(Vf in1, cVf in2, int N)
{
	for (; N; N--)
		*(++in1) -= *(++in2);
}

void scale_vector_f(Vf v, double s, int N)
{
	for ( ; N; N--)
		*++v *= s;
}

void vector_d_to_f(Vf f, cVd d, int n)
{
	for (; n; n--)
		*++f = *++d;
}

void inc_scaled_vector_f(Vf v1, double s, cVf v2, int n)
{
	for (; n; n--)
		*++v1 += s * *++v2;
}

double norm_vector_f(Vf v, int n)
{
	register double d = 0.;
	v++;
	for (; n; n--, v++)
		d += *v * *v;
	return d;
}

double dot_vector_f(Vf v1, cVf v2, int n)
{
	register double d = 0.;
	for (; n; n--)
		d += *++v1 * *++v2;
	return d;
}

// find min and max of a vector v1[1..n]
void min_max_vector_f(double *vmin, double *vmax, cVf v, int n)
{
	*vmin = *vmax = v[1];

	for (; n; n--) {
		++v;
		if (*v > *vmax)	*vmax = *v;
		if (*v < *vmin)	*vmin = *v;
	}
}
