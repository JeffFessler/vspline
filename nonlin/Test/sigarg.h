/***************************** sigarg.h **********************************
*	Definitions for sigjac_*.c
*************************************************************************/

typedef struct
{
	int LL;
	int WW;

	int n_psf;		/* 0 means no blurring			*/
	float *psf;		/* [0..n_psf-1] psf samples		*/

	int *trans_type;	/* [1..M]  each one of TRANS_XXX below	*/
				/* (null indicates no transformations)	*/
} ARGSTRUCT;


#define TRANS_NONE	0
#define TRANS_LOG	1
#define TRANS_MAX	2
#define TRANS_RANGE	3


/*
*	logarithmic transform for non-negative variables
*	p is the original value, x is the internal value
*/
#define LOG_TRANS(p)		log(p)
#define LOG_UNTRANS(x)		exp(x)
#define LOG_DERFAC(p)		(p)


/*
*	transform for [a,b] bounded variables
*/
#define RANGE_TRANS(p,a,b)	log( ((p)-(a)) / ((b)-(p)) )
#define RANGE_UNTRANS(x,a,b)	( (a) + ((b)-(a)) / (1.+exp(-x)) )
#define RANGE_DERFAC(p,a,b)	( ((p)-(a)) / ((b)-(a)) * ((b)-(p)) )


/*
*	transform for [0,m] bound variables
*/
#define MAX_TRANS(p,m)		RANGE_TRANS((p),0,(m))
#define MAX_UNTRANS(x,m)	RANGE_UNTRANS((x),0,(m))
#define MAX_DERFAC(p,m)		RANGE_DERFAC((p),0,(m))
