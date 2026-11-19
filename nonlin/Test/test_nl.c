/*********************** test_nl.c ***************************************
* TEST nonlinear smoothing on SIGJACFUNC
* COMPILE TIME PARAMETERS:
* SIGJACFUNC:	name of signal and jacobian function
* FIX:		if DEF'd, then pre-clean parameters
* FIXFUNC:	function to pre-clean parameters
*************************************************************************/
#define USAGE	{ fprintf(stderr, "Usage: %s info init data [binary?] [est] [mle?]\n", argv[0]); exit(-1); }

#include <stdio.h>
#include "nonlin.h"
#include "sigarg.h"

//	extern void nonlin_est();
//	extern void SIGJACFUNC();

#ifdef FIX
	extern void FIXFUNC();
	static void (*fixfunc)() = FIXFUNC;
#else
#define FIXFUNC	0
#endif


typedef struct
{
	PARAM_SET states;
	MEAS meas;
	ALPHA_SET alpha;
	VECT ti;
	int L, M, N, O, W;
} FILEDATA;

	static FILEDATA fd;
	static Problem prob;
	static char *estfile = (char *) 0;


/*
* FOR TESTING: READ Data FROM FILES
* FILE FORMAT:
* FROM root.info
* L M N O W then alpha1 ... alphaM
* FROM root.init
* for n=1 to N x1n ... xMn
* FROM root.data
* for n=1 to N y1n y2n ... yLn
*/

static void read_data(
FILEDATA *fd,
char *info,
char *init,
char *data,
truf flag_bin,
truf flag_mle)
{
	FILE *fp;
	register int m;
	int LL, MM, NN, OO, WW;

	FOPEN(fp, info, "r")
	if (5 != fscanf(fp, "%d %d %d %d %d", &LL, &MM, &NN, &OO, &WW))	ABORT("rd -2");
	if (LL < MM || MM < 1 || NN < 3)			ABORT("rd -1");

	fd->L = LL;
	fd->M = MM;
	fd->N = NN;
	fd->O = OO;
	fd->W = WW;

	if (!(fd->meas	= ALLOC_MEAS(NN, LL)))		ABORT("rd 1");
	if (!(fd->states= ALLOC_PARAM(NN, MM)))		ABORT("rd 2");

	if (flag_mle)
		fd->alpha = 0;
	else
	{
		if (!(fd->alpha	= ALLOC_ALPHA(MM)))		ABORT("rd 0");
	    for (m = 1; m <= MM; m++)
		if (1 != fscanf(fp, "%g", &(fd->alpha)[m]))	ABORT("rd 3");
	}
	FCLOSE(fp);

	FSCAN_PARAM(fd->states, init, NN, MM);

	if (flag_bin)	READ_MEAS(fd->meas, data, NN, LL);
	else		FSCAN_MEAS(fd->meas, data, NN, LL);
}


int main(int argc, char **argv)
{
	ARGSTRUCT as;
	double cv;
	int dum;
	truf flag_bin, flag_mle;

	if (argc < 4)
		USAGE;

	system("date");

	prob.sigjac = SIGJACFUNC;
	prob.sigarg = (Void *) &as;
	prob.print = TRUE;
	prob.borrow = FALSE;
	prob.modify = 0;
	prob.modarg = (Void *) 0;

	fd.ti = (VECT) 0;
	if (argc > 5)
		estfile = argv[5];
	flag_bin = (argc > 4) ?  atoi(argv[4]) : FALSE;
	flag_mle = (argc > 6) ?  atoi(argv[6]) : FALSE;

	read_data(&fd, argv[1], argv[2], argv[3], flag_bin, flag_mle);

	if (fd.alpha)
		PRINT_ALPHA(fd.alpha, fd.M);

	as.LL = fd.L;
	as.WW = fd.W;
	nonlin_est(&cv, fd.states, fd.alpha, fd.ti, fd.meas, fd.states,
		&prob, JOBCV, fd.L, fd.M, fd.N, fd.O, fd.W, estfile);

	dum = 1;
	FPRINT_PARAM("final", &dum, fd.states, fd.N, fd.M);

	system("date");
	return 0;
}
