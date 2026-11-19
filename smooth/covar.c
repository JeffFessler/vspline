/***************************** covar.c ***********************************
* UTILITY SUBROUTINES FOR COVAR STRUCTURES
* Jeff Fessler
* COPYRIGHT 1989,  THE BOARD OF TRUSTEES OF THE
* LELAND STANFORD JUNIOR UNIVERSITY.  ALL RIGHTS RESERVED.
*************************************************************************/
#include "smooth.h"


/****************************** inv_to_covar() ***************************
* CONVERT INVERSE OF COVARAINCE MATRIX INTO COVAR STRUCTURE
* INPUT:
*	inv: CMAT[1..MxM]
* OUTPUT:
*	cov:
* RETURNS:	1 if OK, 0 if singular
*************************************************************************/
static sof inv_to_covar(COVAR *cov, CMAT inv, cint M)
{
	mat_to_sym_band(cov->inv, inv, M-1, M);
	if ( !chol_sym_band(cov->Linv, cov->Dinv, cov->inv, M-1, M) ) return Failure;
	invert_sym_band(cov->cov, cov->Linv, cov->Dinv, M-1, M);
	cov->trace = band_trace(cov->cov, M-1, M);
	Ok
}


/*
* ALLOCATE SPACE FOR A SET OF N COVARIANCE STRUCTURES
* WITH COVARIANCE SIZE MxM
*/
COVAR *alloc_covar(cint M, int N)
{
	COVAR *cov;
	Mem0(cov, N)
	cov--;

	for (; N >= 1; N--) {
		if ( !(cov[N].cov = ALLOC_MATRIX(M, M)) )	return NULL;
		if ( !(cov[N].cov2 = ALLOC_MATRIX(M, M)) )	return NULL;
		if ( !(cov[N].inv = ALLOC_MATRIX(M, M)) )	return NULL;
		if ( !(cov[N].Linv = ALLOC_MATRIX(M, M)) )	return NULL;
		if ( !(cov[N].Dinv = ALLOC_VECTOR(M)) )		return NULL;
	}
	return cov;
}


void covar_free(COVAR *cov, int N)
{
	for (; N >= 1; N--) {
		FREE_MATRIX(cov[N].cov);
		FREE_MATRIX(cov[N].cov2);
		FREE_MATRIX(cov[N].inv);
		FREE_MATRIX(cov[N].Linv);
		FREE_VECTOR(cov[N].Dinv);
	}
	Free(cov+1);
}


void copy_covar(cCOVAR *src, COVAR *dst, cint MM)
{
	COPY_MATRIX(src->cov, dst->cov, MM, MM);
	COPY_MATRIX(src->cov2, dst->cov2, MM, MM);
	COPY_MATRIX(src->inv, dst->inv, MM, MM);
	COPY_MATRIX(src->Linv, dst->Linv, MM, MM);
	COPY_VECTOR(src->Dinv, dst->Dinv, MM);
}


/*
* MULTIPLY covariance by scale - correcting others
*/
void scale_covar(COVAR *cov, cdouble scale, cint MM)
{
	SCALE_MATRIX(cov->cov, scale, MM, MM);
	SCALE_MATRIX(cov->cov2, scale*scale, MM, MM);
	SCALE_MATRIX(cov->inv, 1. / scale, MM, MM);
	SCALE_VECTOR(cov->Dinv, 1. / scale, MM);
}


/****************************** inv_to_covar_set() ***********************
* CONVERT SET OF COVARIANCE MATRIX INVERSES INTO COVAR STRUCTURE
* INPUT:
*	isigma[1..MxMx(Nc)]		covariance inverses
*	flag_iid			Nc = 1 or N ?
* RETURN:
*	cov[1..(Nc)]:			or 0 if problem
*************************************************************************/
COVAR *inv_to_covar_set(CMAT isigma, cint M, cint N, ctruf flag_iid)
{
	cint Nc = flag_iid ? 1 : N;
	COVAR *covar;

	if ( !(covar = alloc_covar(M, Nc)) )	return NULL;

	for (int nn=Nc; nn >= 1; nn--)
		if (!inv_to_covar(&covar[nn], isigma+(nn-1)*M*M, M)) return NULL;

	return covar;
}
