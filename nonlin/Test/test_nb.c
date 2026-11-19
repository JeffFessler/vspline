/*********************** test_nb.c ***************************************
*	TEST nonlinear smoothing on SIGJACFUNC WITH BEST SCORE
*	USES:
*		atof()
*************************************************************************/
#define USAGE	{ fprintf(stderr, "Usage: %s info init data [binary?] [tol]\n", argv[0]); exit(-1); }

#include <stdio.h>
#include <math.h>
#include "nonlin.h"
#include "sigarg.h"

	static void read_data();

	extern double nonlin_best();
	extern void SIGJACFUNC();


typedef struct
{
	PARAM_SET states;
	MEAS meas;
	ALPHA_SET alpha;
	VECT ti;
	int L, M, N, O, W;
} FILEDATA;


#define TOL	0.01


main(argc, argv)
int argc;
char *argv[];
{
	double tol;
	ARGSTRUCT as;
	FILEDATA fd;
	double best;
	int dum;
	Problem prob;

	prob.sigjac = SIGJACFUNC;
	prob.sigarg = (Void *) &as;
	prob.borrow = FALSE;
	prob.print = TRUE;
	prob.modify = 0;
	prob.modarg = (Void *) 0;

	if (argc < 4)
		USAGE;
	system("date");

	fd.ti = (VECT) 0;

	tol = (argc > 5) ? atof(argv[5]) : TOL;

	read_data(&fd, argv[1], argv[2], argv[3], argc > 4);

	PRINT_ALPHA(fd.alpha, fd.M);

	as.LL = fd.L;
	as.WW = fd.W;
	best = nonlin_best(fd.states, fd.alpha, fd.ti, fd.meas, fd.states, &prob, JOBCV, tol, fd.L, fd.M, fd.N, fd.O, fd.W);

	dum = 0;
	PRINT_ALPHA(fd.alpha, fd.M);
	message("best cv", best);
	FPRINT_PARAM("finalbest", &dum, fd.states, fd.N, fd.M);
	system("date");
}


/*
*	FOR TESTING: READ Data FROM FILES
*	FILE FORMAT:
*	FROM root.info
*	L M N O W then alpha1 ... alphaM
*	FROM root.init
*	for n=1 to N x1n ... xMn
*	FROM root.data
*	for n=1 to N y1n y2n ... yLn
*/

static void read_data(fd, info, init, data, flag_bin)
FILEDATA *fd;
char *info, *init, *data;
int flag_bin;
{
	FILE *fp;
	register int m;
	int LL, MM, NN, OO, WW;

	if (!(fp = fopen(info, "r")))					ABORT("rd -3");
	if (5 != fscanf(fp, "%d %d %d %d %d", &LL, &MM, &NN, &OO, &WW))	ABORT("rd -2");
	if (LL < MM || MM < 1 || NN < 3)				ABORT("rd -1");

	fd->L = LL;
	fd->M = MM;
	fd->N = NN;
	fd->O = OO;
	fd->W = WW;

	if (!(fd->alpha	= ALLOC_ALPHA(MM)))				ABORT("rd 0");
	if (!(fd->meas	= ALLOC_MEAS(NN, LL)))				ABORT("rd 1");
	if (!(fd->states= ALLOC_PARAM(NN, MM)))				ABORT("rd 2");

	for (m = 1; m <= MM; m++)
		if (1 != fscanf(fp, "%g", &(fd->alpha)[m]))		ABORT("rd 3");
	FCLOSE(fp)

	FSCAN_PARAM(fd->states, init, NN, MM);

	if (flag_bin)	READ_MEAS(fd->meas, data, NN, LL);
	else		FSCAN_MEAS(fd->meas, data, NN, LL);
}
