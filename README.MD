#Introduction

+	Spatt includes some 'Large Deviation' samples.  Migrated this over to win64
to play with those samples
+	Included latest gsl at this time.

#Library folder's structure

* `\gsl-2.1-src` - 
* `\gsl-2.1-src\VC14` -   (open libgsl.sln)
* `\spatt-1.2.2.t3` -      (open spatt-1.2.2.t3.sln)

#Build Instructions

+	Build libgsl.sln first. Review that configuration is for a x64 build. This creates a new directory with required include files
\gsl-2.1-src\gsl ( this subdirectory is created during this build)
+	for gsl-2.1-src, I did not run full regression tests
+	Only test was whether it worked with spatt-1.2.2.t3
+	Build spatt-1.2.2.t3.sln
+	Used VC14 (Visual Studio 2015)
+	VC14, Update 1 works fine with relative paths.Previous release was buggy.

#Changes to original code. Spatt

+	In spatt-1.2.2.t3\src  are several file names, renamed with prefix 'main'.
	Each of these have a int main function.  These are test files.
	In VC14, Solution Explorer, you can include or exclude as to test.
	
#Changes to original code. gsl

Made minor changes for build.. you can run a diff against an original gsl 2.1
(1) copied fp-win.c to subdirectory gsl-2.1-src\ieee-utils.
(2) moved structs into file gsl-2.1-src\multifit\gsl_multifit.h
(3) commented out this code in gsl-2.1-src\sys\gsl_sys.h
(4) added include's to files
(5) added macro to cache.h

(2.1) Details .. moved struct lmder_state_t from multifit\lmder.c to multifit\gsl_multifit.h
and added #include <gsl/gsl_multifit.h> to multifit\lmder.c
typedef struct
  {
    size_t iter;
    double xnorm;
    double fnorm;
    double delta;
    double par;
    gsl_matrix *J;             // Jacobian matrix 
    gsl_matrix *r;             // R matrix in J = Q R P^T 
    gsl_vector *tau;
    gsl_vector *diag;          // scaling matrix D = diag(d1,...,dp) 
    gsl_vector *qtf;           // Q^T f 
    gsl_vector *newton;
    gsl_vector *gradient;      // gradient g = J^T f 
    gsl_vector *x_trial;       // trial step x + dx 
    gsl_vector *f_trial;       // trial function f(x + dx) 
    gsl_vector *df;
    gsl_vector *sdiag;
    gsl_vector *rptdx;
    const gsl_vector *weights; // data weights 
    gsl_vector *w;
    gsl_vector *work1;
    gsl_permutation * perm;
  }
lmder_state_t;


(2.2) Details .. moved from multifit\lmniel.c to multifit\gsl_multifit.h
typedef struct
{
  gsl_matrix *A;             /* J^T J */
  gsl_matrix *A_copy;        /* copy of J^T J */
  gsl_matrix *J;             /* Jacobian J(x) */
  gsl_vector *diag;          /* D = diag(J^T J) */
  gsl_vector *rhs;           /* rhs vector = -g = -J^T f */
  gsl_vector *x_trial;       /* trial parameter vector */
  gsl_vector *f_trial;       /* trial function vector */
  gsl_vector *work;          /* workspace length p */
  long nu;                   /* nu */
  double mu;                 /* LM damping parameter mu */
  double tau;                /* initial scale factor for mu */
} lmniel_state_t;


(3) Details... commented out this code in gsl-2.1-src\sys\gsl_sys.h
/* nick
double gsl_log1p (const double x);
double gsl_expm1 (const double x);
double gsl_hypot (const double x, const double y);
double gsl_hypot3 (const double x, const double y, const double z);
double gsl_acosh (const double x);
double gsl_asinh (const double x);
double gsl_atanh (const double x);
*/

(4) added additional include's to these files.
\multifit\lmiterate.c
\multifit\lmmisc.c
\multifit\lmpar.c
\multifit\lmset.c
\multifit\lmutil.c
\multifit\qrsolv.c
\multifit\gsl_multifit.c


(5) added macro to cache.h
#define RETURN_IF_NULL(x) if (!x) { return ; }
