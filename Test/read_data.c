/*********************** read_data.c *************************************
*	FOR TESTING: READ Data TO BE SMOOTHED FROM stdin
*	FILE FORMAT:
*	N M then alpha1 ... alphaM then
*	if (!flag_equ): t1, ..., tN
*	for n=1 to N y1n y2n ... yMn
*	for n=1 to N isig11n .. isigMMn
*	Jeff Fessler
*	COPYRIGHT 1989,  THE BOARD OF TRUSTEES OF THE
*	LELAND STANFORD JUNIOR UNIVERSITY.  ALL RIGHTS RESERVED.
*************************************************************************/
#include "def,test.h"

	extern COVAR *inv_to_covar_set();


/*
*	INPUT:
*		fd->file_iid, file_equ
*/

void read_data(FILEDATA *fd)
{
	int MM, NN, nn, Nc;
	CMAT isig;

	if (2 != scanf("%d %d", &NN, &MM))	ABORT("Bad data file.\n")
	if (NN < 3 || MM < 1)			ABORT("Bad data file.\n")

	fd->N = NN;	fd->M = MM;
	Nc = fd->flag_iid ? 1 : NN;

	if (!fd->flag_equ)
		fd->ti = ALLOC_VECTOR(NN);
	else
		fd->ti = (VECT) 0;
	if (!(fd->alpha	= ALLOC_ALPHA(MM)))		ABORT("Mem 0")
	if (!(fd->yi	= ALLOC_VECTOR(MM*NN)))		ABORT("Mem 1")
	if (!(isig	= ALLOC_VECTOR(Nc*MM*MM)))	ABORT("Mem 2")

	SCAN_ALPHA(fd->alpha, MM);

	if (!fd->flag_equ)
		SCAN_VECTOR(fd->ti, NN);
	SCAN_VECTOR(fd->yi, NN*MM);

	for (nn = 1; nn <= Nc; nn++)
		SCAN_VECTOR(isig+(nn-1)*MM*MM, MM*MM);

	if (!(fd->covar	= inv_to_covar_set(isig, MM, Nc, fd->flag_iid)))
			ABORT("Mem 3")

	FREE_VECTOR(isig);
}
