/***************************** mypowell.c ********************************
*       Powell's method for nonlinear minimization
*       USES:
*       Jeff Fessler
*       COPYRIGHT 1989,  THE BOARD OF TRUSTEES OF THE
*       LELAND STANFORD JUNIOR UNIVERSITY.  ALL RIGHTS RESERVED.
*************************************************************************/
#include "defs_vmatrix.h"

#define MAX_ITER 200


/*
*	Minimize a function of MM variables
*	INPUT:
*	OUTPUT:
*/
void powell(float *p, float **xi, int MM,
float ftol, int *iter, float *fret, double (*func)(float *))
{
	int i,ibig,j;
	float t,fptt,fp,del;
	float *pt,*ptt,*xit;

	TMP_ALLOC_VECTOR_F(pt, MM);
	TMP_ALLOC_VECTOR_F(ptt, MM);
	TMP_ALLOC_VECTOR_F(xit, MM);

	/*
	*	Initial function evaluation
	*/
	*fret=(*func)(p);
	COPY_VECTOR(p, pt, MM);

	/*
	*	Iterate to find minimum
	*/
	for (*iter=1; ; (*iter)++) {
		fp=(*fret);
		ibig=0;
		del=0.0;
		for (i=1;i<=MM;i++) {
			for (j=1;j<=MM;j++) xit[j]=xi[j][i];
			fptt=(*fret);
			linmin(p,xit,MM,fret,func);
			if (fabs(fptt-(*fret)) > del) {
				del=fabs(fptt-(*fret));
				ibig=i;
			}
		}
		if (2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))) {
			TMP_FREE_VECTOR(xit);
			TMP_FREE_VECTOR(ptt);
			TMP_FREE_VECTOR(pt);
			return;
		}
		if (*iter == MAX_ITER) ABORT("Too many iterations in powell");

		for (j=1;j<=MM;j++) {
			ptt[j]=2.0*p[j]-pt[j];
			xit[j]=p[j]-pt[j];
			pt[j]=p[j];
		}

		fptt=(*func)(ptt);
		if (fptt < fp) {
			t=2.0*(fp-2.0*(*fret)+fptt)*Sqr(fp-(*fret)-del)-del*Sqr(fp-fptt);
			if (t < 0.0) {
				linmin(p,xit,MM,fret,func);
				for (j=1;j<=MM;j++) xi[j][ibig]=xit[j];
			}
		}
	}
}
