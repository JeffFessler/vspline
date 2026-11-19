// sub_mat.c
// UTILITY MATRIX SUBROUTINES
// Jeff Fessler
// COPYRIGHT 1989, THE BOARD OF TRUSTEES OF THE
// LELAND STANFORD JUNIOR UNIVERSITY. ALL RIGHTS RESERVED.

#ifdef NOTUSED
/*
* MULTIPLY TWO MATRICES (OR MATRIX AND VECTOR)
* INPUT:
*	in1[1..r1xc1]
*	in2[1..c1xc2]
* OUTPUT:
*	out[1..r1xc2]		out = in1 * in2
* FLOPS:	r1 c2 (2 c1 - 1)
*/

void cmat_mul(out
CMAT out,
cCMAT in1,
cint r1,
cint c1,
cint c2,
cCMAT in2)
{
	register DATA *t;

	for (int i = r1; i >= 1; i--)
	{
		for (int j = c2; j >= 1; j--)
		{
			t = &out[(j-1)*r1 + i];
			*t += in1[(c1-1)*r1 + i]
			    * in2[(j-1)*c1 + c1];

			for (int k = c1-1; k >= 1; k--)
				*t += in1[(k-1)*r1 + i]
				    * in2[(j-1)*c1 + k];
		}
	}
}


/*
* TRACE OF A MATRIX
* INPUT:
*	in[1..nxN]
* RETURNS:
*	trace(in)
* FLOPS: N-1
*/

DATA cmat_trace(
cCMAT in,
cint N)
{
	register DATA tr = in[N*N];

	for (int n = N-1; n >= 1; n--)
		tr += in[(n-1)*N + n];

	return tr;
}


/*
* TAKE THE "TRANSPOSE" OF A SET OF SIGNALS
* INPUT:
*	in[0..NM-1]
* OUTPUT:
*	out[0..MN-1]
*/

void cmat_transpose(
VECT out,
cVECT in,
cint N,
cint M)
{
	for (int n=N-1; n >= 0; n--)
	    for (int m=M-1; m >= 0; m--)
		out[n*M+m] = in[m*N+n];
}
#endif
