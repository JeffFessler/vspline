// matrix_io_d.c
// Double Precision IO SUBROUTINES (mex independent)
// write_xxx	binary to file
// fprint_xxx	ascii to file.ext (or file if !ext)
// fscan_xxx	ascii from file
// read_xxx	binary from file
// USES:  sprintf(), fwrite(), fopen(), fclose(), fflush(), etc
// Jeff Fessler
// COPYRIGHT 1989, THE BOARD OF TRUSTEES OF THE
// LELAND STANFORD JUNIOR UNIVERSITY. ALL RIGHTS RESERVED.

#include "defs-vspline.h"

#define TYPE double
#define cTYPE cdouble

/*
* mat[1..N][1..M]
*/
void print_matrix_d(cTYPE **mat, cint N, cint M)
{
	for (int n=1; n <= N; n++)
		print_vector_d(mat[n], M);
}


/*
* vect[1..N]
*/
sof write_vector_d(cchar *name, cTYPE *vect, cint N)
{
	FILE *fp;
	Fopen0w(fp, name)
	Fwrite0(&vect[1], sizeof(*vect), N, fp)
	Fclose0(fp)
	Ok
}


/*
* mat[1..N][1..M] -> file "name"(.ext)
*/
sof fprint_matrix_d(cchar *name, cint *ext, cTYPE **mat, cint N, cint M)
{
	FILE *fp;
	char fname[256];

	if (ext)	sprintf(fname, "%s.%02d", name, *ext);
	else		sprintf(fname, "%s", name);

	Fopen0w(fp, fname)
	for (int n=1; n <= N; n++)
	{
		for (int m=1; m <= M; m++)
			if (EOF == fprintf(fp, "%g ", mat[n][m]))
				Fail("fprintf")
		if (EOF == fprintf(fp, "\n"))
			Fail("fprintf")
	}
	Fclose0(fp)
	Ok
}


/*
* vect[1..N] from file "name"
*/

sof fscan_vector_d(TYPE *vect, cchar *name, cint N)
{
	FILE *fp;
	Fopen0r(fp, name)
	for (int n=1; n <= N; n++)
		if (1 != fscanf(fp, "%lg", &vect[n]))	Fail("fscanf")
	Fclose0(fp)
	Ok
}


/*
* vect[1..N] from file "name" (binary)
*/

sof read_vector_d(TYPE *vect, cchar *name, cint N)
{
	FILE *fp;

	Fopen0r(fp, name)
	Fread0(&vect[1], sizeof(*vect), N, fp)
	Fclose0(fp)
	Ok
}
