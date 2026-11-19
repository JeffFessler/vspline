/*********************** test_ss.c ***************************************
* TEST best_indiv()
* Jeff Fessler
* COPYRIGHT 1989, THE BOARD OF TRUSTEES OF THE
* LELAND STANFORD JUNIOR UNIVERSITY.  ALL RIGHTS RESERVED.
*************************************************************************/
#include "def,test.h"

#define USAGE	{ fprintf(stderr, "Usage: %s score_id\n", argv[0]); exit(-1);}

	extern double best_indiv();

#define TOL	0.01


int main(int argc, char **argv)
{
	FILEDATA fd;
	VECT trace;
	CMAT ys;
	int score_id;

	if (argc < 2)
		USAGE;
	score_id = atoi(argv[1]);
	fd.flag_iid = 1;
	fd.flag_equ = 0;
	fd.ti = (VECT) 0;

	read_data(&fd);

	if (!(ys = ALLOC_VECTOR(fd.N*fd.M)))	ABORT("Mem 98")
	if (!(trace = ALLOC_VECTOR(fd.M)))	ABORT("Mem 99")

	if (!best_indiv(ys, trace, fd.alpha, fd.ti, fd.yi, score_id, TOL, fd.M, fd.N))
		ABORT("Problem in best_ind")

	PRINT_ALPHA(fd.alpha, fd.M);
	PRINT_VECTOR(trace, fd.M);

	return 0;
}
