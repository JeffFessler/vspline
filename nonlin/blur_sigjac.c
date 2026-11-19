/************************ blur_sigjac() **********************************
* Blur signals and jacobians by a point spread function.
*
* PARAMETERS: (M = 2) (tau,rad)
* USES: fabs(), sqrt(), asin(), ceil(), floor();
*************************************************************************/
#include "nonlin.h"


/*
* dot product of two vectors
*/

static float dot_f(cfloat *f1, cfloat *f2, int n)
{
	register float sum = 0.;
	for (n--; n >= 0; n--)
		sum += *f1++ * *f2++;
	return sum;
}

/*
* Convolve: v1 <= v1 convolved with v2
*	v1[0..n1-1]
*	v2[0..n2-1]
*	work[0..n1+n2-1]
*/

static void blur_f(float *v1, int n1, float *v2, int n2, float *work)
{
	register int n;

	n = n1+n2-1;

	Bzero(work, n1+2*(n2-1));
	Bcopy(v1, work+n2-1, n1);
	Bzero(v1, n);

	for (n--; n>=0; n--)
		v1[n] = dot_f(work+n, v2, n2);
}


/*
* in
*	sig: [1..P][1..W]		Computed signal
*	jac: [1..M][1..P][1..W]		Jacobians
*	off:
*	psf[0..K-1]			Blur samples
*	K:				number of blur samples (odd!)
* out
*	sig: [1..P][1..W]		Blurred signal
*	jac: [1..M][1..P][1..W]		Blurred jacobians
*	off: [1..M][1..P][1..W]		Corrected offset
*/

sof blur_sigjac(
SIGNAL sig,
JACOBIAN jac,
OFFSET *off,
int MM, int PP, int WW, int LL,
float *psf,
int KK)
{
	int pp;
	SIGNAL_TYPE *tmp;
	TMP_ALLOC_VECTOR_F(tmp, WW+2*(KK-1));

	for (pp=PP; pp >= 1; pp--) {
		int mm;
		int count = off->count[pp];

		off->start[pp] -= (KK-1)/2;
		off->count[pp] += KK-1;

		// Ensure that the blurred size isn't too large
		if (off->count[pp] > WW || off->start[pp] < 1 || off->count[pp]+off->start[pp]-1 > LL)
			Fail("too large");

		// Blur signals and Jacobians
		blur_f(&sig[pp][1], count, psf, KK, tmp);
		for (mm=MM; mm >= 1; mm--)
			blur_f(&jac[mm][pp][1], count, psf, KK, tmp);
	}

	TMP_FREE_VECTOR(tmp);

	Ok
}
