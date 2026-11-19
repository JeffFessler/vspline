/*********************** test_vs.c ***************************************
* Read data from file (see read_data.c), smooth it, compute scores.
* USES atof()
* Jeff Fessler
* COPYRIGHT 1989, THE BOARD OF TRUSTEES OF THE
* LELAND STANFORD JUNIOR UNIVERSITY.  ALL RIGHTS RESERVED.
*************************************************************************/
#include "def,test.h"

#define USAGE	{ fprintf(stderr, "Usage: %s iid equi [score tol]\n", argv[0]); exit(-1);}

#define TOL	0.01

int main(int argc, char **argv)
{
	FILEDATA fd;
	double ur, cv, gcv, penin, penout, rough;
	CMAT ys;
	truf flag_best;

	if (argc < 3)
		USAGE;
	flag_best = argc > 3;
	fd.flag_iid = atoi(argv[1]);
	fd.flag_equ = atoi(argv[2]);
	fd.ti = (VECT) 0;

	read_data(&fd);

	PRINT_ALPHA(fd.alpha, fd.M);
	if (!fd.flag_equ)
		PRINT_VECTOR(fd.ti, fd.N);
	PRINT_VECTOR(fd.yi, fd.N*fd.M);

	PRINT_BAND(fd.covar[1].inv, fd.M, fd.M-1);

	if (!(ys = ALLOC_VECTOR(fd.N*fd.M)))    ABORT("Mem 99")

    if (flag_best)
    {
	int score_id = atoi(argv[3]);
	double tol = (argc > 4) ? atof(argv[4]) : TOL;
	double bs = best_score(ys, fd.alpha, fd.ti, fd.yi, fd.covar, fd.flag_iid, fd.M, fd.N, score_id, tol, True);

	message0("best score ", bs);
	PRINT_ALPHA(fd.alpha, fd.M);
	PRINT_VECTOR(ys, fd.M*fd.N);
    }
    else
    {
	if ( !vect_smooth(ys, &ur, &cv, &gcv, &penin, &penout, &rough,
		fd.alpha, fd.ti, fd.yi, fd.covar, fd.flag_iid, fd.M, fd.N) )
		ABORT("E 7")

	PRINT_VECTOR(ys, fd.N*fd.M);
	message("ur score", ur);
	message("cv score", cv);
	message("gcv score", gcv);
	message("roughness", rough);
	message("penalty yi (in)", penin);
	message("penalty ys (out)", penout);
    }

	return 0;
}
