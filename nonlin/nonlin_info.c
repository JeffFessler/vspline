/***************************** nonlin_info.c *****************************
*	allocation subroutines for "info" structures
*************************************************************************/
#include "nonlin.h"


/*
* ALLOCATE Nonlin STRUCTURE
*/
Nonlin *alloc_info(int MM, int nn, int PP, int WW, ctruf flag_mod)
{
	Nonlin *info;

	Alloc0(info, Nonlin, 1, sizeof(Nonlin))
	if (!(info->par		= ALLOC_PARAM(nn, MM)))			ABORT("")
	if (!(info->xh		= ALLOC_PARAM(nn, MM)))			ABORT("")
	if (!(info->off		= ALLOC_OFFSET_SET(nn, PP)))		ABORT("")
	if (!(info->sig_mod	= ALLOC_SIGNAL_SET(nn, PP, WW)))	ABORT("")
	if (!(info->resid	= ALLOC_SIGNAL_SET(nn, PP, WW)))	ABORT("")
	if (!(info->jac_mod	= ALLOC_JACOB_SET(nn, MM, PP, WW)))	ABORT("")
	if (!(info->covar	= alloc_covar(MM, nn)))			ABORT("")

	/*
	* For the typical case with no "modify" function,
	* set sig_mod equal to sig_raw so that both pointers work
	*/
	if (flag_mod) {
	    if (!(info->sig_raw	= ALLOC_SIGNAL_SET(nn, PP, WW)))	ABORT("")
	    if (!(info->jac_raw	= ALLOC_JACOB_SET(nn, MM, PP, WW)))	ABORT("")
	}
	else {
		info->sig_raw = info->sig_mod;
		info->jac_raw = info->jac_mod;
	}

	return info;
}


/*
* FREE Nonlin STRUCTURE
*/
void free_info(Nonlin *info, int nn)
{
	FREE_PARAM(info->par);
	FREE_PARAM(info->xh);
	FREE_OFFSET_SET(info->off);
	FREE_SIGNAL_SET(info->resid);
	covar_free(info->covar, nn);

	if (info->sig_raw == info->sig_mod) {
		FREE_SIGNAL_SET(info->sig_raw);
	} else {
		FREE_SIGNAL_SET(info->sig_raw);
		FREE_SIGNAL_SET(info->sig_mod);
	}

	if (info->jac_raw == info->jac_mod) {
		FREE_JACOB_SET(info->jac_raw, nn);
	} else {
		FREE_JACOB_SET(info->jac_raw, nn);
		FREE_JACOB_SET(info->jac_mod, nn);
	}
	Free(info);
}


/*
* COPY ONE Nonlin STRUCTURE'S CONTENTS TO ANOTHER
*/
void copy_info(Nonlin *src, Nonlin *dst, int MM, int nn, int PP, int WW)
{
	int ii;

	dst->bad = src->bad;

	COPY_PARAM(src->par, dst->par, nn, MM);
	COPY_PARAM(src->xh, dst->xh, nn, MM);
	COPY_OFFSET_SET(src->off, dst->off, nn, PP);
	COPY_SIGNAL_SET(src->resid, dst->resid, nn, PP, WW);

	COPY_SIGNAL_SET(src->sig_mod, dst->sig_mod, nn, PP, WW);
	if (src->sig_raw != src->sig_mod)
		COPY_SIGNAL_SET(src->sig_raw, dst->sig_raw, nn, PP, WW);

	COPY_JACOB_SET(src->jac_mod, dst->jac_mod, nn, MM, PP, WW);
	if (src->jac_raw != src->jac_mod)
		COPY_JACOB_SET(src->jac_raw, dst->jac_raw, nn, MM, PP, WW);

	for (ii=nn; ii >= 1; ii--)
		COPY_MATRIX(src->covar[ii].inv, dst->covar[ii].inv, MM, MM);
}
