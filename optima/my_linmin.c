/***************************** mylinmin.c ********************************
*       One dimensional minimization
*       USES:
*       Jeff Fessler
*       COPYRIGHT 1989,  THE BOARD OF TRUSTEES OF THE
*       LELAND STANFORD JUNIOR UNIVERSITY.  ALL RIGHTS RESERVED.
*************************************************************************/
#include "defs_vmatrix.h"

#define TOL 2.0e-4

	static int MM_GLOBAL;
	double (*func_GLOBAL)(float *);
	float *pnt_GLOBAL;
	float *dir_GLOBAL;


/*
*
*/
static double f1dim(float x)
{
	float f,*xt;

	TMP_ALLOC_VECTOR_F(xt, MM_GLOBAL);
	COPY_VECTOR(pnt_GLOBAL, xt, MM_GLOBAL);
	inc_scaled_vector_f(xt, x, dir_GLOBAL, MM_GLOBAL);
	f = (*func_GLOBAL)(xt);

	TMP_FREE_VECTOR(xt);
	return f;
}


/*
*	Minimize a function of MM variables along one direction
*	INPUT:
*	OUTPUT:
*/
void linmin(float pnt[], float dir[], int MM, float *fret, double (*func)(float *))
{
	int j;
	float xx,xmin,fx,fb,fa,bx,ax;

	MM_GLOBAL  = MM;
	TMP_ALLOC_VECTOR_F(pnt_GLOBAL, MM);
	TMP_ALLOC_VECTOR_F(dir_GLOBAL, MM);

	/*
	*	Save globals
	*/
	func_GLOBAL = func;
	COPY_VECTOR(pnt, pnt_GLOBAL, MM);
	COPY_VECTOR(dir, dir_GLOBAL, MM);

	ax=0.0;
	xx=1.0;
	bx=2.0;
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
	*fret=brent(ax,xx,bx,f1dim,TOL,&xmin);
	for (j=1; j <= MM; j++) {
		dir[j] *= xmin;
		pnt[j] += dir[j];
	}

	TMP_FREE_VECTOR(dir_GLOBAL);
	TMP_FREE_VECTOR(pnt_GLOBAL);
}
