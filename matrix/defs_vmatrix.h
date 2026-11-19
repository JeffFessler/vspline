// matrix.h
// General purpose vector and matrix utilites.
// I use the "Numerical Recipes in C" convention for allocating
// vectors and matrices, so
// A[a..b][c..d] is a matrix with b-a+1 rows and d-c+1 columns.
// A[1..nxm]: an nxm "column" matrix stored 1,1 2,1 ... n,1 ... n,m.
// However, not only are they allocated as a doubly subscripted array,
// i.e., mat[1..N][1..M], but the data is stored sequentially
// as well.  So &mat[1][1] is a pointer to nxM elements stored
// in the order:
// mat[1][1], mat[1][2], ..., mat[1][M], mat[2][1], ... mat[N][M].
// Thus we can easily operate on the individual rows of the matrix,
// or on the entire matrix.
// Jeff Fessler
// COPYRIGHT 1989, THE BOARD OF TRUSTEES OF THE
// LELAND STANFORD JUNIOR UNIVERSITY. ALL RIGHTS RESERVED.

#include "defs-vspline.h"

typedef double *Vd;
typedef cdouble *cVd;
typedef float *Vf;
typedef cfloat *cVf;

// alloc.c
extern double *alloc_vector_d(cint);
extern double **alloc_matrix_d(int, cint);
extern float *alloc_vector_f(int);
extern float **alloc_matrix_f(int, cint);
extern float ***alloc_matrix3_f(cint, cint, cint);

// debug.c
// see defs-vspline.h

// ident.c
extern float **identity_f(int);
extern double **identity_d(int);

// vector_d.c
extern void vector_add_d(Vd out, cVd in1, cVd in2, cint N);
extern void vector_sub_d(Vd out, cVd in1, cVd in2, int N);
extern void vector_inc_d(Vd in1, cVd in2, int N);
extern void vector_dec_d(Vd in1, cVd in2, int N);
extern void scale_vector_d(Vd v, cdouble s, int N);
extern void vector_f_to_d(Vd d, cVf f, int n);
extern void inc_scaled_vector_d(Vd v1, cdouble s, cVd v2, int n);
extern double norm_vector_d(cVd v, int n);
extern double dot_vector_d(cVd v1, cVd v2, int n);
extern void min_max_vector_d(double *vmin, double *vmax, cVd v, int n);

// vector_f.c
extern void vector_add_f(Vf out, cVf in1, cVf in2, int N);
extern void vector_sub_f(Vf out, cVf in1, cVf in2, int N);
extern void vector_inc_f(Vf in1, cVf in2, int N);
extern void vector_dec_f(Vf in1, cVf in2, int N);
extern void scale_vector_f(Vf v, double s, int N);
extern void vector_d_to_f(Vf d, cVd f, int n);
extern void inc_scaled_vector_f(Vf v1, double s, cVf v2, int n);
extern double norm_vector_f(Vf v, int n);
extern double dot_vector_f(Vf v1, cVf v2, int n);
extern void min_max_vector_f(double *vmin, double *vmax, cVf v, int n);


// precision independent
#define COPY_VECTOR(src,dst,n)		Bcopy(&src[1],&dst[1],n)
#define COPY_MATRIX(src,dst,r,c)	COPY_VECTOR(src[1],dst[1],(r)*(c))
#define COPY_MATRIX3(src,dst,r,c,l)	COPY_VECTOR(src[1][1],dst[1][1],(r)*(c)*(l))
#define ZERO_VECTOR(v,n)		Bzero(&v[1],n)
#define ZERO_MATRIX(m,r,c)		ZERO_VECTOR(m[1],(r)*(c))
#define ZERO_MATRIX3(m,r,c,l)		ZERO_VECTOR(m[1][1],(r)*(c)*(l))
#define FREE_VECTOR(v)			{ Free(&v[1]); }
#define FREE_MATRIX(m)			{ Free(&m[1][1]); Free(&m[1]); }
#define FREE_MATRIX3(m)			{ Free(&m[1][1][1]); Free(&m[1][1]); Free(&m[1]); }

// double precision
#define ALLOC_VECTOR_D			alloc_vector_d
#define ALLOC_MATRIX_D			alloc_matrix_d
#define ALLOC_MATRIX3_D			alloc_matrix3_d
#define SCALE_VECTOR_D			scale_vector_d
#define NORM_VECTOR_D			norm_vector_d
#define DOT_VECTOR_D			dot_vector_d
#define INC_SCALED_VECTOR_D		inc_scaled_vector_d
#define ADD_VECTOR_D			vector_add_d
#define INC_VECTOR_D			vector_inc_d
#define SUB_VECTOR_D			vector_sub_d
#define DEC_VECTOR_D			vector_dec_d

#define SCAN_VECTOR_D			scan_vector_d
#define READ_VECTOR_D			read_vector_d
#define FSCAN_VECTOR_D			fscan_vector_d
#define PRINT_VECTOR_D			print_vector_d
#define PRINT_MATRIX_D			print_matrix_d
#define WRITE_VECTOR_D			write_vector_d
#define FPRINT_MATRIX_D			fprint_matrix_d

#define SCALE_MATRIX_D(m,s,r,c)		SCALE_VECTOR_D(m[1],s,(r)*(c))
#define SCAN_MATRIX_D(m,r,c)		SCAN_VECTOR_D(m[1],(r)*(c))
#define READ_MATRIX_D(m,f,r,c)		READ_VECTOR_D(m[1],f,(r)*(c))
#define FSCAN_MATRIX_D(m,f,r,c)		FSCAN_VECTOR_D(m[1],f,(r)*(c))
#define WRITE_MATRIX_D(name,m,r,c)	WRITE_VECTOR_D(name,m[1],(r)*(c))

// single precision
#define ALLOC_VECTOR_F			alloc_vector_f
#define ALLOC_MATRIX_F			alloc_matrix_f
#define ALLOC_MATRIX3_F			alloc_matrix3_f
#define SCALE_VECTOR_F			scale_vector_f
#define NORM_VECTOR_F                   norm_vector_f
#define DOT_VECTOR_F                    dot_vector_f
#define INC_SCALED_VECTOR_F             inc_scaled_vector_f
#define ADD_VECTOR_F                    vector_add_f
#define INC_VECTOR_F                    vector_inc_f
#define SUB_VECTOR_F                    vector_sub_f
#define DEC_VECTOR_F                    vector_dec_f

#define SCAN_VECTOR_F			scan_vector_f
#define READ_VECTOR_F			read_vector_f
#define FSCAN_VECTOR_F			fscan_vector_f
#define PRINT_VECTOR_F			print_vector_f
#define PRINT_MATRIX_F			print_matrix_f
#define WRITE_VECTOR_F			write_vector_f
#define FPRINT_MATRIX_F			fprint_matrix_f

#define SCALE_MATRIX_F(m,s,r,c)		SCALE_VECTOR_F(m[1],s,(r)*(c))
#define SCAN_MATRIX_F(m,r,c)		SCAN_VECTOR_F(m[1],(r)*(c))
#define READ_MATRIX_F(m,f,r,c)		READ_VECTOR_F(m[1],f,(r)*(c))
#define FSCAN_MATRIX_F(m,f,r,c)		FSCAN_VECTOR_F(m[1],f,(r)*(c))
#define WRITE_MATRIX_F(name,m,r,c)	WRITE_VECTOR_F(name,m[1],(r)*(c))


// macros for allocating and freeing temporary work space

// precision independent
#define TMP_FREE_VECTOR(v)	TMP_Free(&v[1])
#define TMP_FREE_MATRIX(m)	{ TMP_Free(&m[1][1]); TMP_Free(&m[1]); }
#define TMP_FREE_MATRIX3(m)	{ TMP_Free(&m[1][1][1]); TMP_Free(&m[1][1]); TMP_Free(&m[1]); }

// double precision
#define TMP_ALLOC_VECTOR_D(v,n) {v = (double *) TMP_Alloc((n),sizeof(*v)); \
				if (!v) ABORT("alloc failed"); \
				Bzero(v, (n)); \
				v--; \
				}
#define TMP_ALLOC_MATRIX_D(m,r,c)	{ int i; \
				m = (double **) TMP_Alloc((r),sizeof(*m)); \
				if (!m) ABORT("alloc failed"); \
				m--; \
				TMP_ALLOC_VECTOR_D(m[1], (r)*(c)) \
				for (i=(r); i > 1; i--) \
					m[i] = (i-1)*(c) + m[1]; \
				}

// single precision
#define TMP_ALLOC_VECTOR_F(v,n) {v = (float *) TMP_Alloc((n),sizeof(*v)); \
				if (!v) ABORT("alloc failed"); \
				Bzero(v, (n)); \
				v--; \
				}
#define TMP_ALLOC_MATRIX_F(m,r,c)	{ int i; \
				m = (float **) TMP_Alloc((r),sizeof(*m)); \
				if (!m)	ABORT("alloc failed"); \
				m--; \
				TMP_ALLOC_VECTOR_F(m[1], (r)*(c)) \
				for (i=(r); i > 1; i--) \
					m[i] = (i-1)*(c) + m[1]; \
				}
#define TMP_ALLOC_MATRIX3_F(m,r,c,l)	{ int i; \
				m = (float ***) TMP_Alloc((r),sizeof(*m)); \
				if (!m)		ABORT("alloc failed"); \
				m--; \
				m[1] = (float **) TMP_Alloc((r)*(c),sizeof(**m)); \
				if (!m[1])	ABORT("alloc failed"); \
				m[1]--; \
				TMP_ALLOC_VECTOR_F(m[1][1], (r)*(c)*(l)); \
				for (i=(r); i > 1; i--) \
					m[i] = (i-1)*(c) + m[1]; \
				for (i=(r)*(c); i > 1; i--) \
					m[1][i] = (i-1)*(l) + m[1][1]; \
				}
