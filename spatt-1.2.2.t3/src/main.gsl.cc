#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
/*
#define finite gsl_finite

const gsl_min_fminimizer_type  * gsl_min_fminimizer_goldensection;
const gsl_min_fminimizer_type  * gsl_min_fminimizer_brent;

double fn1 (double x, void * params)
{
  return cos(x) + 1.0;
}

int
main (void)
{
  int status;
  int iter = 0, max_iter = 100;
  // gsl_min_fminimizer_type *T;
  gsl_min_fminimizer *s;
  double m = 2.0, m_expected = M_PI;
  double a = 0.0, b = 6.0;
  gsl_function F;

  F.function = &fn1;
  F.params = 0;

  gsl_min_fminimizer_type T = { "test", 10 };

  s = gsl_min_fminimizer_alloc (&T);
  gsl_min_fminimizer_set (s, &F, m, a, b);

  printf ("using %s method\n",
          gsl_min_fminimizer_name (s));

  printf ("%5s [%9s, %9s] %9s %10s %9s\n",
          "iter", "lower", "upper", "min",
          "err", "err(est)");

  printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
          iter, a, b,
          m, m - m_expected, b - a);

  do
    {
      iter++;
      status = gsl_min_fminimizer_iterate (s);

      m = gsl_min_fminimizer_x_minimum (s);
      a = gsl_min_fminimizer_x_lower (s);
      b = gsl_min_fminimizer_x_upper (s);

      status 
        = gsl_min_test_interval (a, b, 0.001, 0.0);

      if (status == GSL_SUCCESS)
        printf ("Converged:\n");

      printf ("%5d [%.7f, %.7f] "
              "%.7f %+.7f %.7f\n",
              iter, a, b,
              m, m - m_expected, b - a);
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_min_fminimizer_free (s);

  return status;
}
*/
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>

const gsl_min_fminimizer_type  * gsl_min_fminimizer_goldensection;
const gsl_min_fminimizer_type  * gsl_min_fminimizer_brent;

double fn1(double x, void * params)
{
	return cos(x) + 1.0;
}

int
main(void)
{
	int status;
	int iter = 0, max_iter = 100;
	const gsl_min_fminimizer_type *T;
	gsl_min_fminimizer *s;
	double m = 2.0, m_expected = M_PI;
	double a = 0.0, b = 6.0;
	gsl_function F;

	F.function = &fn1;
	F.params = 0;

	T = new gsl_min_fminimizer_type();
	s = gsl_min_fminimizer_alloc(T);
	gsl_min_fminimizer_set(s, &F, m, a, b);

	printf("using %s method\n",
		gsl_min_fminimizer_name(s));

	printf("%5s [%9s, %9s] %9s %10s %9s\n",
		"iter", "lower", "upper", "min",
		"err", "err(est)");

	printf("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
		iter, a, b,
		m, m - m_expected, b - a);

	do
	{
		iter++;
		status = gsl_min_fminimizer_iterate(s);

		m = gsl_min_fminimizer_minimum(s);
		a = gsl_min_fminimizer_x_lower(s);
		b = gsl_min_fminimizer_x_upper(s);

		status
			= gsl_min_test_interval(a, b, 0.001, 0.0);

		if (status == GSL_SUCCESS)
			printf("Converged:\n");

		printf("%5d [%.7f, %.7f] "
			"%.7f %.7f %+.7f %.7f\n",
			iter, a, b,
			m, m_expected, m - m_expected, b - a);
	} while (status == GSL_CONTINUE && iter < max_iter);

	return status;
}
/*
bash$ ./a.out
0 [0.0000000, 6.0000000] 2.0000000 -1.1415927 6.0000000
1 [2.0000000, 6.0000000] 3.2758640 +0.1342713 4.0000000
2 [2.0000000, 3.2831929] 3.2758640 +0.1342713 1.2831929
3 [2.8689068, 3.2831929] 3.2758640 +0.1342713 0.4142862
4 [2.8689068, 3.2831929] 3.2758640 +0.1342713 0.4142862
5 [2.8689068, 3.2758640] 3.1460585 +0.0044658 0.4069572
6 [3.1346075, 3.2758640] 3.1460585 +0.0044658 0.1412565
7 [3.1346075, 3.1874620] 3.1460585 +0.0044658 0.0528545
8 [3.1346075, 3.1460585] 3.1460585 +0.0044658 0.0114510
9 [3.1346075, 3.1460585] 3.1424060 +0.0008133 0.0114510
10 [3.1346075, 3.1424060] 3.1415885 -0.0000041 0.0077985
Converged:
11 [3.1415885, 3.1424060] 3.1415927 -0.0000000 0.0008175
*/