// decs-smooth.h


// band_util.c
extern double band_trace(cSYMBAND A, int J, int M);
extern double trace_prod(cSYMBAND A, cSYMBAND B, cint M);
extern void mult_band_vect(VECT y, cSYMBAND A, cVECT x, cint J, cint N);
extern void mult_full_vect(VECT v, cSYMBAND A, cVECT z, cint M);
extern void mult_full_band(SYMBAND C, cSYMBAND A, cSYMBAND B, cint M);
extern void sym_band_to_mat(CMAT A, cSYMBAND B, cint J, cint N);
extern void mat_to_sym_band(SYMBAND B,
cCMAT A, cint J, cint N);

// best_indiv.c
extern double best_indiv(VECT ys, VECT trace, ALPHA_SET alpha, VECT ti,
	VECT yi, cint l_score_id, cdouble tol, cint MM, cint nn);

// best_mse.c
extern double best_mse(CMAT l_ys, CMAT l_yr, ALPHA_SET alpha,
	VECT ti, CMAT yi, COVAR * covar,
	ctruf flag_iid, cint M, cint N, cdouble tol);

// best_score.c
extern double best_score(CMAT ys, ALPHA_SET alpha,
	VECT ti, CMAT yi, COVAR *covar,
	ctruf flag_iid, cint M, cint N, cint l_score_id,
	cdouble tol, cint chat);

// covar.c
extern COVAR *alloc_covar(cint M, int N);
extern void covar_free(COVAR *cov, int N);
extern void copy_covar(cCOVAR *src, COVAR *dst, cint MM);
extern void scale_covar(COVAR *cov, cdouble scale, cint MM);
extern COVAR *inv_to_covar_set(CMAT isigma, cint M, cint N, ctruf flag_iid);

// cv_score.c
extern sof cv_setup_alloc(void);
extern sof cv_setup_covar(void);
extern void cv_cleanup(void);
extern sof cross_val(double *ur, double *cv, double *gcv, cint job);
extern sof trace_i_a(double *trace, ALPHA_SET alpha);

// cv_sub.c
extern void select_nth_qbiqt(SYMBAND F, cSYMBAND IB, cSYMBAND Qt,
	cint M, cint N, cint n);
extern double rss_term(cVECT err, cSYMBAND inv, cint M);
extern double ur_term(cSYMBAND Fn, cSYMBAND cov2, cint M);
extern double cv_term(cVECT err, cSYMBAND FL, cVECT FD, cSYMBAND inv, cint M);
extern double gcv_term(cSYMBAND Fn, cSYMBAND cov, cint M);

// example1.c
extern sof vspline(CMAT ys, ALPHA_SET alpha, VECT ti, CMAT yi,
	CMAT isigma, ctruf flag_iid, cint MM, cint nn);
extern sof vect_smooth(CMAT ys, double *ur, double *cv, double *gcv,
	double *penin, double *penout, double *rough, ALPHA_SET alpha,
	VECT ti, CMAT yi, COVAR *covar, truf flag_iid, cint MM, cint nn);

// interp.c
extern sof cubic_interp(VECT ys, cVECT ti, cVECT yi, cVECT xs, cint nn, cint nx);

// sym_band.c
extern void kron_diag(SYMBAND K, cSYMBAND B, cint J, cint N, cVECT D, cint M);
extern sof chol_sym_band(BAND L, VECT D, cSYMBAND B, cint J, cint N);
extern void solve_sym_band(VECT x, cBAND L, cVECT D, VECT z,
	cint J, cint N, cint job);
extern void invert_sym_band(SYMBAND IB, cSYMBAND L, cVECT D, cint J, cint N);

// vs_init.c
extern sof spline_mat(SYMBAND R, BAND Qt, cVECT ti, cint N);
extern void form_qt_Y(VECT QtY, cBAND Qt, cVECT Y, cint M, cint N);
extern void form_qt_s_q(SYMBAND QtSQ, cBAND Qt, cCOVAR *covar,
	cint M, cint N, ctruf flag_iid);

// vs_proc.c
extern void form_band(SYMBAND sum, cSYMBAND R_A, cSYMBAND QtSQ, cint M, cint N);
extern void form_q_c(VECT q_c, cBAND Qt, cVECT c, cint M, cint N);
extern void form_sig_q_c(OMATRIX sig_q_c, cCOVAR *covar, cCMAT q_c,
	cint M, cint N, ctruf flag_iid);

extern sof vs_setup_alloc(int M, int N);
extern sof vs_setup_ti(VECT ti, int N);
extern sof vs_setup_yi(VECT yi, int M, int N);
extern sof vs_setup_cov(COVAR *covar, truf flag_iid);
extern void vs_cleanup(void);
extern sof vs_setup_alpha(ALPHA_SET alpha);
extern sof vs_smooth(CMAT result, int job);
extern sof vs_transfer(void);
extern double roughness(void);
extern double post_penalty(void);
extern sof penalty(double *pen, ALPHA_SET alpha, VECT yi, int M, int N);
