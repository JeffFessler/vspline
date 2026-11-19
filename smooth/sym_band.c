/****************************** sym_band.c *******************************
*	SUBROUTINES FOR SYMMETRIC BAND MATRIX MANIPULATIONS
*	SEE band.h FOR BAND MATRIX STORAGE FORMAT.
*	NOTE: BAND(N-1,N) MEANS THESE ROUTINES WORK FOR "FULL" MATRICES.
*	Jeff Fessler
*	COPYRIGHT 1989,  THE BOARD OF TRUSTEES OF THE
*	LELAND STANFORD JUNIOR UNIVERSITY.  ALL RIGHTS RESERVED.
*************************************************************************/
#include "smooth.h"


/*
*	SUBTRACT n POINT InnER PRODUCT OF p1, p2 FROM x
*	INPUT:
*		x
*		p1[0..n-1]
*		p2[0..n-1]
*	OUTPUT:
*		x		x -=  <p1, p2>
*/
static void dotdiff(DATA *x, DATA *p1, DATA *p2, int n)
{
	for (; n > 0; n--)
		*x -= *p1++ * *p2++;
}


/*
*	SUBTRACT n POINT InnER PRODUCT OF p1, p2, p3
*	FROM init
*	pi[0..n-1]
*/
static double ddot123(double sum, VECT p1, VECT p2, VECT p3, int n)
{
	for (; n > 0; n--, p1++, p2++, p3++)
		sum -= *p1 * *p2 * *p3;

	return sum;
}


static truf istiny(DATA x)
{
	return (x == 0.);
}


/**************************** kron_diag() *******************************
*	TENSOR PRODUCT OF SYMMETRIC BAND MATRIX WITH DIAGONAL MATRIX
*	INPUT:
*		B SYMBAND(J,N)
*		D[1..M]			DIAGONAL OF DIAGONAL MATRIX
*	OUTPUT:
*		K SYMBAND((J+1)M-1,NM)	K = B @ Diag(D)
*	FLOPS:	M[(N-J-1)(J+1) + (J+2)(J+1)/2] = M[N(J+1) - J(J+1)/2]
*		IF J=1, N<=N-2, THEN  N(2M) - 5M
************************************************************************/
void kron_diag(SYMBAND K, cSYMBAND B, cint J, cint N, cVECT D, cint M)
{
	cint nk = (J+1) * M;

	for (int m=M; m >= 1; m--)		/* UP THE DIAGONAL */
	{
	    for (int i=J+1; i >= 1; i--)
		for (int j=i; j >= 1; j--)
			K[(i-1)*M+m][(j-i)*M+nk] = B[i][j-i+J+1] * D[m];

	    for (int i=J+2; i <= N; i++)
		for (int j=i-J; j <= i; j++)
			K[(i-1)*M+m][(j-i)*M+nk] = B[i][j-i+J+1] * D[m];
	}
}


#ifdef NOTUSED
/**************************** kron_sym_band() ***************************
*	KRONECKER TENSOR PRODUCT OF SYMMETRIC BAND MATRIX WITH SYMMETRIC
*	FULL MATRIX.
*	INPUT:
*		B SYMBAND(J,N)
*		A[1..MM]	M*M MATRIX STORED COLUMN-WISE.
*	OUTPUT:
*		K SYMBAND((J+1)M-1,NM)
*	FLOPS:	MM[(N-J-1)J + J(J+1)/2] +NM((M-1)/2 = N[MMJ+M(M-1)/2]-MM/2(JJ+J)
*		OLD: M*M*[(J+1)*(J+2)/2 + (N-J)*(J+1)]
*		IF J=1, N<=N-2, THEN  N(3/2MM+M/2)-4MM-M
************************************************************************/
void kron_sym_band(
SYMBAND K,
cSYMBAND B,
cint J,
cint N,
cCMAT A,
cint M)
{
	register int i, j, m, n, nk;

	nk = (J+1) * M;

    for (m=M; m >= 1; m--)	/* DOWN THE MATRIX A	*/
    {
	/* THE BELOW DIAGONAL BLOCKS	*/
	for (n=M; n >= 1; n--)	/* ACROSS  THE MATRIX A	*/
	{
	    for (i=J+1; i >= 1; i--)
		for (j=i-1; j >= 1; j--)
			K[(i-1)*M+m][(j-i)*M+(n-m)+nk] = B[i][j-i+J+1] * A[(n-1)*M+m];

	    for (i=J+2; i <= N; i++)
		for (j=i-J; j < i; j++)
			K[(i-1)*M+m][(j-i)*M+(n-m)+nk] = B[i][j-i+J+1] * A[(n-1)*M+m];
	}

	/* THE ON DIAGONAL BLOCKS	*/
	for (n=m; n >= 1; n--)
	    for (i=N; i >= 1; i--)
			K[(i-1)*M+m][        (n-m)+nk] = B[i][    J+1] * A[(n-1)*M+m];
    }
}
#endif // NOTUSED


/****************************** chol_sym_band() **************************
*	CHOLESKY DECOMPOSITION OF A nxN SYMMETRIC BAND MATRIX B (2*J+1)
*	J MUST BE AT LEAST 1 FOR THIS TO WORK
*	INPUT:
*		B: SYMBAND(J,N)		-> L * diag(D) * trans(L)
*	OUTPUT:
*		L: BAND(J, N)		UNIT LOWER TRIANGULAR
*		D[1..N]
*	RETURNS:
*		1 IF SUCCESFUL, 0 IF THE MATRIX IS SINGULAR.
*	FLOPS:	(N-1)(3/2JJ+5/2J) - [JJJ+JJ/2-3/2J] = N(3/2JJ+5/2J)-J(J+1)^2
*		J==N-1 => N(N-1)(N-2)/2
*		IF J=3M-1 AND N<=M(N-2), THEN N(27/2MMM -3/2MM -M)-54MMM+12MM+2M
*************************************************************************/
sof chol_sym_band(BAND L, VECT D, cSYMBAND B, cint J, cint N)
{
	register int n, k, Jn, nJ;

	// LOOP DOWN THROUGH ROWS OF BAND MATRIX STRUCTURE

	D[1] = B[1][J+1];
	if (istiny(D[1]))	return Failure;

	if (N == 1)
		return 1;

    for (n=2; n <= J+1; n++)
    {
	Jn = J-n;

	L[n][Jn+2] = B[n][Jn+2] / D[1];

	for (k=2; k < n; k++)
		L[n][k+Jn+1] = ddot123((double) B[n][k+Jn+1],
			&D[1], &L[n][Jn+2], &L[k][J-k+2], k-1) / D[ k ];

	D[n] = ddot123((double) B[n][J+1], &D[1], &L[n][Jn+2], &L[n][Jn+2], n-1);
	if (istiny(D[n]))	return Failure;
    }

    for (; n <= N; n++)
    {
	nJ = n-J;

	L[n][1] = B[n][1] / D[ nJ ];

	for (k=2; k <= J; k++)
		L[n][k] = ddot123((double) B[n][k],
			&D[nJ], &L[n][1], &L[nJ+k-1][J-k+2], k-1) / D[ nJ+k-1 ];

	D[n] = ddot123((double) B[n][J+1], &D[nJ], &L[n][1], &L[n][1], J);
	if (istiny(D[n]))	return Failure;
    }

	Ok
}


/****************************** solve_sym_band() *************************
* Solve (L D L') x = z
* SEE Golub and Van Loan Alg. 5.3-2,3
* in
*	L: BAND(J, N)
*	D[1..N]
*	z[1..N]	caution: can be modified!
*	job:	if 1, return in z, if 0 return in x
* out
*	x[1..n]
* flops:	N+2[J(J-1)+2J(N-J)] = N(1+4J) -2JJ -2J
*		J==N-1 => 2N^2 - N
*		IF J=3M-1, N<=M(N-2)  THEN N(12MM-3M) -42MM +12M
*************************************************************************/
void solve_sym_band(VECT x, cBAND L, cVECT D,
VECT z, // trick
cint J, cint N, cint job)
{
	register int i, j;

	if (job)
		x = z;

	// NOW x IS REALLY u	(SOLVING Lu = z)
	else
		COPY_VECTOR(z, x, N);

	if (N == 1)
	{
		x[1] /= D[1];
		return;
	}

	for (i=2; i <= J; i++)
	{
		// x[i] = z[i];
		dotdiff(&x[i], &L[i][J-i+2], &x[1], i-1);
	}
	for (   ; i <= N; i++)
	{
		// x[i] = z[i];
		dotdiff(&x[i], &L[i][1], &x[i-J], J);
	}

	// NOW x IS REALLY v	(SOLVING Dv = u)
	for (i=1; i <= N; i++)
		x[i] /= D[i];

	// NOW IS REALLY x ON LHS, v ON RHS (SOLVING L' x = v)
	// x[N] = x[N];

	for (i=N-1; i >= N-J; i--)
		for (j=i+1; j <= N; j++)
			x[i] -= L[j][i-j+J+1] * x[j];
	for (; i >= 1; i--)
		for (j=i+1; j <= i+J; j++)
			x[i] -= L[j][i-j+J+1] * x[j];
}



/****************************** invert_sym_band() ************************
*	FIND CENTRAL INVERSE OF A SYMMETRIC BAND MATRIX B = L * D * L'
*	See Hutchinson and deHoog, Numerische Mathematik, 1985.
*	IB = inv(B)
*	INPUT:
*		L: SYMBAND(J,N)
*		D[1..N]
*	OUTPUT:
*		IB: SYMBAND(J,N)	central bands of inv(B)
*	FLOPS: 2JJN + ...
*	IF J=3M-1, N->M(N-2) THEN 18MMMN+...
*************************************************************************/
void invert_sym_band(SYMBAND IB, cSYMBAND L, cVECT D, cint J, cint N)
{
	register int n, l, k, nJ;

	IB[N][J+1] = 1. / D[N];

    for (n=N-1 ; n >= 1; n--)
    {
	nJ = (J < N-n) ? J : N-n;

	for (l=1; l <= nJ; l++)
	{
		IB[n+l][-l+J+1] = 0.;		/* WASTES A FLOP...*/
		for (k=1; k <= l; k++)
			IB[n+l][-l+J+1] -= L[n+k][-k+J+1] * IB[n+l][k-l+J+1];
		for (   ; k <= nJ; k++)
			IB[n+l][-l+J+1] -= L[n+k][-k+J+1] * IB[n+k][l-k+J+1];
	}

	IB[n][J+1] = 1. / D[n];

	for (l=1; l <= nJ; l++)
		IB[n][J+1] -= L[n+l][-l+J+1] * IB[n+l][-l+J+1];
    }
}
