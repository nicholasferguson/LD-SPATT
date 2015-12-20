/* min/test.c
*
* Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or (at
* your option) any later version.
*
* This program is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/
#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include "ldstat.h"
#define finite gsl_finite

const gsl_min_fminimizer_type  * gsl_min_fminimizer_goldensection;
const gsl_min_fminimizer_type  * gsl_min_fminimizer_brent;

gsl_function create_function(double(*f)(double, void *));

void
test_f_e(const gsl_min_fminimizer_type * T,
	const char * description, gsl_function *f,
	double lower_bound, double minimum, double upper_bound,
	double correct_minimum);

void
test_f(const gsl_min_fminimizer_type * T,
	const char * description, gsl_function *f,
	double lower_bound, double middle, double upper_bound,
	double correct_minimum);

int
test_bracket(const char * description, gsl_function *f, double lower_bound,
	double upper_bound, unsigned int max);

double f_cos(double x, void * p);
double func1(double x, void * p);
double func2(double x, void * p);
double func3(double x, void * p);
double func4(double x, void * p);

#define SAFE_FUNC_CALL(f, x, yp) \
do { \
  *yp = GSL_FN_EVAL(f,x); \
  if (!finite(*yp)) \
    GSL_ERROR("computed function value is infinite or NaN", GSL_EBADFUNC); \
} while (0)




/* stopping parameters */

const double EPSABS = 0.001;
const double EPSREL = 0.001;

double rate();

const unsigned int MAX_ITERATIONS = 100;

void my_error_handler(const char *reason, const char *file,
	int line, int err);

#define WITHIN_TOL(a, b, epsrel, epsabs) \
 (fabs((a) - (b)) < (epsrel) * GSL_MIN(fabs(a), fabs(b)) + (epsabs))

int
main(void)
{
	gsl_function F_cos, F_func1, F_func2, F_func3, F_func4;

	const gsl_min_fminimizer_type * fminimizer[4];
	const gsl_min_fminimizer_type ** T;

	gsl_ieee_env_setup();

	rate();

	fminimizer[0] = gsl_min_fminimizer_goldensection;
	fminimizer[1] = gsl_min_fminimizer_brent;
	fminimizer[2] = 0;

	F_cos = create_function(f_cos);
	F_func1 = create_function(func1);
	F_func2 = create_function(func2);
	F_func3 = create_function(func3);
	F_func4 = create_function(func4);

	gsl_set_error_handler(&my_error_handler);

	for (T = fminimizer; *T != 0; T++)
	{
		test_f(*T, "cos(x) [0 (3) 6]", &F_cos, 0.0, 3.0, 6.0, M_PI);
		test_f(*T, "x^4 - 1 [-3 (-1) 17]", &F_func1, -3.0, -1.0, 17.0, 0.0);
		test_f(*T, "sqrt(|x|) [-2 (-1) 1.5]", &F_func2, -2.0, -1.0, 1.5, 0.0);
		test_f(*T, "func3(x) [-2 (3) 4]", &F_func3, -2.0, 3.0, 4.0, 1.0);
		test_f(*T, "func4(x) [0 (0.782) 1]", &F_func4, 0, 0.782, 1.0, 0.8);

		test_f_e(*T, "invalid range check [4, 0]", &F_cos, 4.0, 3.0, 0.0, M_PI);
		test_f_e(*T, "invalid range check [1, 1]", &F_cos, 1.0, 1.0, 1.0, M_PI);
		test_f_e(*T, "invalid range check [-1, 1]", &F_cos, -1.0, 0.0, 1.0, M_PI);
	}
	test_bracket("cos(x) [1,2]", &F_cos, 1.0, 2.0, 15);
	test_bracket("sqrt(|x|) [-1,0]", &F_func2, -1.0, 0.0, 15);
	test_bracket("sqrt(|x|) [-1,-0.6]", &F_func2, -1.0, -0.6, 15);
	test_bracket("sqrt(|x|) [-1,1]", &F_func2, -1.0, 1.0, 15);

	gsl_test_summary();
	return 0;
}

void
test_f(const gsl_min_fminimizer_type * T,
	const char * description, gsl_function *f,
	double lower_bound, double middle, double upper_bound,
	double correct_minimum)
{
	int status;
	size_t iterations = 0;
	double m, a, b;
	double x_lower, x_upper;
	gsl_min_fminimizer * s;

	x_lower = lower_bound;
	x_upper = upper_bound;

	s = gsl_min_fminimizer_alloc(T);
	gsl_min_fminimizer_set(s, f, middle, x_lower, x_upper);

	do
	{
		iterations++;

		status = gsl_min_fminimizer_iterate(s);

		m = gsl_min_fminimizer_x_minimum(s);
		a = gsl_min_fminimizer_x_lower(s);
		b = gsl_min_fminimizer_x_upper(s);

#ifdef DEBUG
		printf("%.12f %.18f %.12f %.18f %.12f %.18f\n",
			a, GSL_FN_EVAL(f, a), m, GSL_FN_EVAL(f, m), b, GSL_FN_EVAL(f, b));
#endif

		if (a > b)
			gsl_test(GSL_FAILURE, "interval is invalid (%g,%g)", a, b);

		if (m < a || m > b)
			gsl_test(GSL_FAILURE, "m lies outside interval %g (%g,%g)", m, a, b);

		if (status) break;

		status = gsl_min_test_interval(a, b, EPSABS, EPSREL);
	} while (status == GSL_CONTINUE && iterations < MAX_ITERATIONS);

	gsl_test(status, "%s, %s (%g obs vs %g expected) ",
		gsl_min_fminimizer_name(s), description,
		gsl_min_fminimizer_x_minimum(s), correct_minimum);

	/* check the validity of the returned result */

	if (!WITHIN_TOL(m, correct_minimum, EPSREL, EPSABS))
	{
		gsl_test(GSL_FAILURE, "incorrect precision (%g obs vs %g expected)",
			m, correct_minimum);
	}

	gsl_min_fminimizer_free(s);

}

void
test_f_e(const gsl_min_fminimizer_type * T,
	const char * description, gsl_function *f,
	double lower_bound, double middle, double upper_bound,
	double correct_minimum)
{
	int status;
	size_t iterations = 0;
	double x_lower, x_upper;
	double a, b;
	gsl_min_fminimizer * s;

	x_lower = lower_bound;
	x_upper = upper_bound;

	s = gsl_min_fminimizer_alloc(T);
	status = gsl_min_fminimizer_set(s, f, middle, x_lower, x_upper);

	if (status != GSL_SUCCESS)
	{
		gsl_min_fminimizer_free(s);
		gsl_test(status == GSL_SUCCESS, "%s, %s", T->name, description);
		return;
	}

	do
	{
		iterations++;
		gsl_min_fminimizer_iterate(s);
		a = gsl_min_fminimizer_x_lower(s);
		b = gsl_min_fminimizer_x_upper(s);

		status = gsl_min_test_interval(a, b, EPSABS, EPSREL);
	} while (status == GSL_CONTINUE && iterations < MAX_ITERATIONS);

	gsl_test(!status, "%s, %s", gsl_min_fminimizer_name(s), description,
		gsl_min_fminimizer_x_minimum(s) - correct_minimum);

	gsl_min_fminimizer_free(s);
}

void
my_error_handler(const char *reason, const char *file, int line, int err)
{
	if (0)
		printf("(caught [%s:%d: %s (%d)])\n", file, line, reason, err);
}

int
test_bracket(const char * description, gsl_function *f, double lower_bound,
	double upper_bound, unsigned int max)
{
	int status;
	double x_lower, x_upper;
	double f_upper, f_lower, f_minimum;
	double x_minimum;

	x_lower = lower_bound;
	x_upper = upper_bound;
	SAFE_FUNC_CALL(f, x_lower, &f_lower);
	SAFE_FUNC_CALL(f, x_upper, &f_upper);
	status = gsl_min_find_bracket(f, &x_minimum, &f_minimum, &x_lower, &f_lower, &x_upper, &f_upper, max);
	gsl_test(status, "%s, interval: [%g,%g], values: (%g,%g), minimum at: %g, value: %g",
		description, x_lower, x_upper, f_lower, f_upper, x_minimum, f_minimum);
	return status;
}

gsl_function create_function(double(*f)(double, void *))
{
	gsl_function F;
	F.function = f;
	F.params = 0;
	return F;
}

double
f_cos(double x, void * p)
{
	p = 0;  /* avoid warning about unused parameter */
	return cos(x);
}

/* f(x) = x^4 - 1 */
/* minimum at x = 0 */

double
func1(double x, void * p)
{
	p = 0;  /* avoid warning about unused parameter */
	return pow(x, 4.0) - 1;
}

/* f(x) = sqrt(|x|) */
/* minimum at x = 0 */

double
func2(double x, void * p)
{
	p = 0;  /* avoid warning about unused parameter */
	return sqrt(fabs(x));
}


/* f(x) = 1 for x < 1 and -exp(-x) for x >= 1 */
/* minimum at x = 1 */

double
func3(double x, void * p)
{
	p = 0;  /* avoid warning about unused parameter */

	if (x < 1)
		return 1;
	else
		return -exp(-x);
}

/* f(x) = x - 30/(1+1e5*(x-0.8)**2) */
/* minimum near x = 0.8 */

double
func4(double x, void * p)
{
	p = 0;  /* avoid warning about unused parameter */

	return x - 30.0 / (1.0 + 1e5 * pow(x - 0.8, 2.0));
}
double fn1(double x, void * params)
{
	return cos(x) + 1.0;
}
double rate() {
	//  double ldstat::sparse_rate(){

	double ax, bx, cx, fa, fb, fc;
	int status;
	int iter = 0, max_iter = 100;
//	const gsl_min_fminimizer_type *T;
	gsl_min_fminimizer_type T = { "test", 10 };
	gsl_min_fminimizer *s;
	gsl_function F;


	//F.function = &ldstat_sparse_F;
	F.function = &fn1;
	F.params = 0;

		/* brackting the miniminum */
		ax = -1.0; bx = -2.0; cx = -2.5;
		fb = 1.0; fc = 2.0; fa = 2.0;
	
		/* brent minimizer */
	//	T = new gsl_min_fminimizer_type;

		s = gsl_min_fminimizer_alloc(&T);
		gsl_min_fminimizer_set_with_values(
				s,  //gsl_min_fminimizer
				&F,  //gsl_function
				bx,  //x_minimum  x_minimum < x_upper || x_minimum > x_lower
				fb, //f_minimum
				cx,  //x_lower must be lt x_upper
				fc,  //f_lower  (f_minimum < f_lower || f_minimum < f_upper)
				ax, //x_upper
				fa  //f_upper
			);
		/* main loop */
		do
		{
			iter++;
			status = gsl_min_fminimizer_iterate(s);

			bx = gsl_min_fminimizer_x_minimum(s);
			ax = gsl_min_fminimizer_x_lower(s);
			cx = gsl_min_fminimizer_x_upper(s);
			fb = gsl_min_fminimizer_f_minimum(s);
			fa = gsl_min_fminimizer_f_lower(s);
			fc = gsl_min_fminimizer_f_upper(s);
			//printf("[%e,%e,%e] with function values = [%e,%e,%e]\n",ax,bx,cx,fa,fb,fc);	
			//	  printf("rel errors: %e\t%e\n",(fa-fb)/fb,(fc-fb)/fb);
			status = gsl_min_test_interval(ax, cx, 0.001, 0.0);

			//  if (status == GSL_SUCCESS)
			//    printf ("Converged:\n");
			//  
			// printf("%i\tF(%e)=%e\n",iter,bx,fb);
		} while (status == GSL_CONTINUE && iter < max_iter);
	


	return fb;
};
