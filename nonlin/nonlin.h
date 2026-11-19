/***************************** nonlin.h **********************************
*	Definitions for nonlinear vector spline smoothing routines.
*************************************************************************/
#include "smooth.h"

#define MEAS_TYPE	float
typedef MEAS_TYPE **MEAS;
#define	NORM_MEAS(m,r,c)	norm_vector_f(m[1],(r)*(c))
#define	FSCAN_MEAS		FSCAN_MATRIX_F
#define	FPRINT_MEAS		FPRINT_MATRIX_F
#define	WRITE_MEAS		WRITE_MATRIX_F
#define	READ_MEAS		READ_MATRIX_F
#define	COPY_MEAS		COPY_MATRIX
#define	ALLOC_MEAS		ALLOC_MATRIX_F
#define	FREE_MEAS		FREE_MATRIX
#define	TMP_ALLOC_MEAS		TMP_ALLOC_MATRIX_F
#define	TMP_FREE_MEAS		TMP_FREE_MATRIX


#define PARAM_TYPE	double
typedef PARAM_TYPE *PARAM;		/* [1..M]			*/
typedef PARAM_TYPE **PARAM_SET;		/* [1..N][1..M]			*/
#define ALLOC_PARAM		ALLOC_MATRIX_D
#define FREE_PARAM		FREE_MATRIX
#define COPY_PARAM		COPY_MATRIX
#define FPRINT_PARAM		FPRINT_MATRIX_D
#define FSCAN_PARAM		FSCAN_MATRIX_D
#define PARAMSTR		"%lg"
#define PARAM_TO_INFO(info, param, m, n)	COPY_MATRIX((param), info->par, (m), (n))
#define INFO_TO_PARAM(param, info, m, n)	COPY_MATRIX(info->par, (param), (m), (n))


#define SIGNAL_TYPE	float
typedef SIGNAL_TYPE **SIGNAL;		/* [1..P][1..W]			*/
typedef SIGNAL *SIGNAL_SET;		/* [1..N][1..P][1..W]		*/

#define ALLOC_SIGNAL_SET	ALLOC_MATRIX3_F
#define COPY_SIGNAL_SET		COPY_MATRIX3
#define ZERO_SIGNAL_SET		ZERO_MATRIX3
#define FREE_SIGNAL_SET		FREE_MATRIX3


#define JACOB_TYPE	float
typedef JACOB_TYPE ***JACOBIAN;		/* [1..M][1..P][1..W]		*/
typedef JACOBIAN *JACOB_SET;		/* [1..N][1..M][1..P][1..W]	*/

#define ALLOC_JACOB_SET		alloc_jacob_set
#define COPY_JACOB_SET		copy_jacob_set
#define FREE_JACOB_SET		free_jacob_set
#define ALLOC_JACOB		ALLOC_MATRIX3_F
#define COPY_JACOB		COPY_MATRIX3
#define ZERO_JACOB		ZERO_MATRIX3
#define FREE_JACOB		FREE_MATRIX3
#define DOT_JACOB		DOT_VECTOR_F


/*	THIS STRUCTURE EXPLOITS THE SPARSENESS OF ARTERIAL PROJECTIONS	*/
typedef struct
{
	short *start;			/* [1..P] where pulse begins	*/
	short *count;			/* [1..P] length of pulse	*/
} OFFSET;

typedef OFFSET *OFFSET_SET;		/* [1..N]			*/

	extern JACOB_SET alloc_jacob_set(int, int, int, int);
	extern OFFSET_SET alloc_offset_set(int, int);
	extern void free_jacob_set(JACOB_SET, int), free_offset_set(OFFSET_SET);
	extern void copy_jacob_set(JACOB_SET, JACOB_SET, int, int, int, int);
	extern void copy_offset_set(OFFSET_SET, OFFSET_SET, int, int);


#define ALLOC_OFFSET_SET	alloc_offset_set
#define COPY_OFFSET_SET		copy_offset_set
#define ZERO_OFFSET_SET
#define FREE_OFFSET_SET		free_offset_set


/*	This structure is used to carry around problem dependent info	*/

typedef struct
{
	/* REQUIRED	*/
	void (*sigjac)(void);		/* Function returning sig & jac	*/
	truf borrow;			/* ok to borrow jacobians?	*/

	/* RECOMMENDED	*/
	Void *sigarg;			/* arguments to that function	*/
	truf print;			/* print debugging statements?	*/

	/* OPTIONAL	*/
	void (*modify)(void);		/* modify signals (from overlap)*/
	Void *modarg;			/* args to the modify function	*/
} Problem;

typedef int SIGJAC_RETURN;


/*
*	The "badness" criterion of an estimate
*	|Y-S|   = (       |Y|(S==0)       ) + |Y-S|(S!=0)
*	|Y-S|   = |Y|       - |Y|(S!=0)     + |Y-S|(S!=0)
*	norm_resid = meas_norm + norm_loc_meas + norm_loc_resid
*/
typedef struct
{
	double meas_norm;		/* MEASUREMENT NORM |Y|		*/
	double norm_loc_meas;		/* - |Y|(S!=0)			*/
	double norm_loc_resid;		/* |Y-S|(S!=0)			*/
	double roughness;		/* ROUGHNESS PENALTY TERM	*/
} Badness;

#define BADDIFF(a,b)	(	(a.meas_norm-b.meas_norm) + \
				(a.norm_loc_meas-b.norm_loc_meas) + \
				(a.norm_loc_resid-b.norm_loc_resid) + \
				(a.roughness-b.roughness) )


/*
*	This structure carries around the "state" of an estimate
*/
typedef struct
{
	/*
	*	"DATA" precision matrices
	*/
	OMATRIX par;			/* [1..N][1..M]			*/
	OMATRIX xh;			/* [1..N][1..M]			*/
	COVAR *covar;			/* [1..N]			*/

	/*
	*	float precision matrices
	*/
	OFFSET_SET off;			/* pulse signal offsets		*/
	SIGNAL_SET resid;		/* (meas - sig_mod)		*/
	SIGNAL_SET sig_raw;		/* raw signals			*/
	SIGNAL_SET sig_mod;		/* modified signals		*/
	JACOB_SET jac_raw;		/* raw jacobians		*/
	JACOB_SET jac_mod;		/* modified jacobians		*/

	Badness bad;			/* see above			*/
} Nonlin;

extern Nonlin *alloc_info(int MM, int nn, int PP, int WW, ctruf flag_mod);
extern void free_info(Nonlin *info, int nn);
extern void copy_info(Nonlin *src, Nonlin *dst, int MM, int nn, int PP, int WW);

extern sof blur_sigjac(SIGNAL sig, JACOBIAN jac, OFFSET *off,
	int MM, int PP, int WW, int LL, float *psf, int KK);
