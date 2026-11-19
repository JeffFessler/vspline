/*********************** test_gen.c **************************************
* Generate signals for later estimation
*************************************************************************/
#define Usage	{ fprintf(stderr, "Usage: %s info par measout [binary?]\n", argv[0]); exit(-1); }

#include <stdio.h>
#include "nonlin.h"
#include "sigarg.h"

	extern void SIGJACFUNC();


typedef struct
{
	PARAM_SET states;
	MEAS meas;
	int L, M, N, O, W;
} FILEDATA;


#define ARG_INFO	argv[1]
#define ARG_PAR		argv[2]
#define ARG_MEAS	argv[3]
#define ARG_BIN		argv[4]



/*
* READ Data FROM FILES
* FILE FORMAT:
* FROM info
*	L M N O W
* FROM par
*	for n=1 to N x1n ... xMn
*/

static void read_data(
FILEDATA *fd,
char *info,
char *init)
{
	FILE *fp;
	int LL, MM, NN, OO, WW;

	if (!(fp = fopen(info, "r")))				ABORT("rd -3");
	if (5 != fscanf(fp, "%d %d %d %d %d", &LL, &MM, &NN, &OO, &WW))	ABORT("rd -2");
	if (LL < MM || MM < 1 || NN < 3)			ABORT("rd -1");

	fd->L = LL;
	fd->M = MM;
	fd->N = NN;
	fd->O = OO;
	fd->W = WW;

	if (!(fd->meas	= ALLOC_MEAS(NN, LL)))		ABORT("rd 1");
	if (!(fd->states= ALLOC_PARAM(NN, MM)))		ABORT("rd 2");

	FSCAN_PARAM(fd->states, init, NN, MM);
}


int main(int argc, char *argv[])
{
	ARGSTRUCT as;
	FILEDATA fd;
	SIGNAL_SET sig;
	JACOB_SET jac;
	OFFSET_SET off;
	SIGJAC_RETURN out;
	int nn;
	int LL, MM, NN, OO, WW;

	if (argc < 4)
		Usage

	read_data(&fd, ARG_INFO, ARG_PAR);
	LL=fd.L; MM=fd.M; NN=fd.N; OO=fd.O; WW=fd.W;

	as.LL = LL;
	as.WW = WW;

	sig = ALLOC_SIGNAL_SET(NN, OO, WW);
	jac = ALLOC_JACOB_SET(NN, MM, OO, WW);
	off = ALLOC_OFFSET_SET(NN, OO);

	SIGJACFUNC(&out, sig, jac, off, fd.states, (Void *) &as, MM, NN, OO, WW);

	for (nn=1; nn <= NN; nn++)
	{
		int oo, ww;
		for (oo=OO; oo >= 1; oo--)
		{
			// COMPUTE resid[nn][oo]
			SIGNAL_TYPE *p1 = &fd.meas[nn][off[nn].start[oo]];
			SIGNAL_TYPE *p2 = &sig[nn][oo][1];
			for (ww=off[nn].count[oo]; ww >= 1; ww--)
				*p1++ = *p2++;
		}
	}

	if ((argc > 4) ? atoi(argv[4]) : 0)
		WRITE_MEAS(ARG_MEAS, fd.meas, LL, NN);
	else
		FPRINT_MEAS(ARG_MEAS, (int *) 0, fd.meas, LL, NN);
}
