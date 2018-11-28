#include <stdio.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

//  f1 =  -dt1 + sqrt(((x-s2_x)^2)/(s2_v00^2)+((y-s2_y)^2)/(s2_v90^2)) - sqrt(((x-s1_x)^2)/(s1_v00^2)+((y-s1_y)^2)/(s1_v90^2));
//  f2 =  -dt2 + sqrt(((x-s3_x)^2)/(s3_v00^2)+((y-s3_y)^2)/(s3_v90^2)) - sqrt(((x-s1_x)^2)/(s1_v00^2)+((y-s1_y)^2)/(s1_v90^2));

struct mparams
  {
    double dt1;
    double dt2;
    double s1_x;
    double s1_y;
    double s2_x;
    double s2_y;
    double s3_x;
    double s3_y;
    double s1_v00;
    double s1_v90;
    double s2_v00;
    double s2_v90;
    double s3_v00;
    double s3_v90;
  };

int mmea_f (const gsl_vector * x, void *params,
              gsl_vector * f)
{
  double dt1 = ((struct mparams *) params)->dt1;
  double dt2 = ((struct mparams *) params)->dt2;
  double s1_x = ((struct mparams *) params)->s1_x;
  double s1_y = ((struct mparams *) params)->s1_y;
  double s2_x = ((struct mparams *) params)->s2_x;
  double s2_y = ((struct mparams *) params)->s2_y;
  double s3_x = ((struct mparams *) params)->s3_x;
  double s3_y = ((struct mparams *) params)->s3_y;
  double s1_v00 = ((struct mparams *) params)->s1_v00;
  double s1_v90 = ((struct mparams *) params)->s1_v90;
  double s2_v00 = ((struct mparams *) params)->s2_v00;
  double s2_v90 = ((struct mparams *) params)->s2_v90;
  double s3_v00 = ((struct mparams *) params)->s3_v00;
  double s3_v90 = ((struct mparams *) params)->s3_v90;

  const double x0 = gsl_vector_get (x, 0);
  const double x1 = gsl_vector_get (x, 1);

  //  f1 =  -dt1 + sqrt(((x-s2_x)^2)/(s2_v00^2)+((y-s2_y)^2)/(s2_v90^2)) - sqrt(((x-s1_x)^2)/(s1_v00^2)+((y-s1_y)^2)/(s1_v90^2));
  const double y0 = -dt1 + sqrt( ((x0-s2_x)*(x0-s2_x))/((s2_v00)*(s2_v00)) + ((x1-s2_y)*(x1-s2_y))/((s2_v90)*(s2_v90)) ) - sqrt(((x0-s1_x)*(x0-s1_x))/((s1_v00)*(s1_v00))+((x1-s1_y)*(x1-s1_y))/((s1_v90)*(s1_v90)));
  //  f2 =  -dt2 + sqrt(((x-s3_x)^2)/(s3_v00^2)+((y-s3_y)^2)/(s3_v90^2)) - sqrt(((x-s1_x)^2)/(s1_v00^2)+((y-s1_y)^2)/(s1_v90^2));
  const double y1 = -dt2 + sqrt( ((x0-s3_x)*(x0-s3_x))/((s3_v00)*(s3_v00)) + ((x1-s3_y)*(x1-s3_y))/((s3_v90)*(s3_v90)) ) - sqrt(((x0-s1_x)*(x0-s1_x))/((s1_v00)*(s1_v00))+((x1-s1_y)*(x1-s1_y))/((s1_v90)*(s1_v90)));

  gsl_vector_set (f, 0, y0);
  gsl_vector_set (f, 1, y1);

  return GSL_SUCCESS;
}
struct rparams
  {
    double a;
    double b;
  };

int rosenbrock_f (const gsl_vector * x, void *params,
              gsl_vector * f)
{
  double a = ((struct rparams *) params)->a;
  double b = ((struct rparams *) params)->b;

  const double x0 = gsl_vector_get (x, 0);
  const double x1 = gsl_vector_get (x, 1);

  const double y0 = a * (1 - x0);
  const double y1 = b * (x1 - x0 * x0);

  gsl_vector_set (f, 0, y0);
  gsl_vector_set (f, 1, y1);

  return GSL_SUCCESS;
}

int rosenbrock_df (const gsl_vector * x, void *params,
               gsl_matrix * J)
{
  const double a = ((struct rparams *) params)->a;
  const double b = ((struct rparams *) params)->b;

  const double x0 = gsl_vector_get (x, 0);

  const double df00 = -a;
  const double df01 = 0;
  const double df10 = -2 * b  * x0;
  const double df11 = b;

  gsl_matrix_set (J, 0, 0, df00);
  gsl_matrix_set (J, 0, 1, df01);
  gsl_matrix_set (J, 1, 0, df10);
  gsl_matrix_set (J, 1, 1, df11);

  return GSL_SUCCESS;
}

int rosenbrock_fdf (const gsl_vector * x, void *params,
                gsl_vector * f, gsl_matrix * J)
{
  rosenbrock_f (x, params, f);
  rosenbrock_df (x, params, J);

  return GSL_SUCCESS;
}

int print_state_f (size_t iter, gsl_multiroot_fsolver * s)
{
  printf ("iter = %3u x = % .3f % .3f "
          "f(x) = % .3e % .3e\n",
          iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->f, 0),
          gsl_vector_get (s->f, 1));
}

int print_state_d (size_t iter, gsl_multiroot_fdfsolver * s)
{
  printf ("iter = %3u x = % .3f % .3f "
          "f(x) = % .3e % .3e\n",
          iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->f, 0),
          gsl_vector_get (s->f, 1));
}

/*
 * Bessel Function Example
 */
/*
int main (void)
{
  double z = 5.0;
  double y = gsl_sf_bessel_J0 (z);
  printf ("J0(%g) = %.18e\n", z, y);
  return 0;
}
*/

/*
 *  MMEA Solver/Algorithm Example
 */ 
int main (void)
{
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;

  int status;
  size_t i, iter = 0;

  const size_t n = 2;
  struct mparams{
    double dt1;
    double dt2;
    double s1_x;
    double s1_y;
    double s2_x;
    double s2_y;
    double s3_x;
    double s3_y;
    double s1_v00;
    double s1_v90;
    double s2_v00;
    double s2_v90;
    double s3_v00;
    double s3_v90;
  };
  struct mparams p = {1.0, 10.0};
  gsl_multiroot_function f = {&mmea_f, n, &p};

  double x_init[2] = {-1.0, -1.0};
  gsl_vector *x = gsl_vector_alloc (n);

  gsl_vector_set (x, 0, x_init[0]);
  gsl_vector_set (x, 1, x_init[1]);

  T = gsl_multiroot_fsolver_hybrids;
  s = gsl_multiroot_fsolver_alloc (T, 2);
  gsl_multiroot_fsolver_set (s, &f, x);

  print_state_f (iter, s);

  do
    {
      iter++;
      status = gsl_multiroot_fsolver_iterate (s);

      print_state_f (iter, s);

      if (status)   /* check if solver is stuck */
        break;

      status = 
        gsl_multiroot_test_residual (s->f, 1e-7);
    }
  while (status == GSL_CONTINUE && iter < 1000);

  printf ("status = %s\n", gsl_strerror (status));

  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x);
  return 0;
}

/*
 *  Rosenbrock Example
 */
/*
int main (void)
{
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;

  int status;
  size_t i, iter = 0;

  const size_t n = 2;
  struct rparams p = {1.0, 10.0};
  gsl_multiroot_function f = {&rosenbrock_f, n, &p};

  double x_init[2] = {-10.0, -5.0};
  gsl_vector *x = gsl_vector_alloc (n);

  gsl_vector_set (x, 0, x_init[0]);
  gsl_vector_set (x, 1, x_init[1]);

  T = gsl_multiroot_fsolver_hybrids;
  s = gsl_multiroot_fsolver_alloc (T, 2);
  gsl_multiroot_fsolver_set (s, &f, x);

  print_state_f (iter, s);

  do
    {
      iter++;
      status = gsl_multiroot_fsolver_iterate (s);

      print_state_f (iter, s);

      if (status)
        break;

      status = 
        gsl_multiroot_test_residual (s->f, 1e-7);
    }
  while (status == GSL_CONTINUE && iter < 1000);

  printf ("status = %s\n", gsl_strerror (status));

  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x);
  return 0;
}
*/

/*
 *  Rosenbrock Example Devirative
 */
/*
int main(){
  const gsl_multiroot_fdfsolver_type *T;
  gsl_multiroot_fdfsolver *s;

  int status;
  size_t i, iter = 0;

  const size_t n = 2;
  struct rparams p = {1.0, 10.0};
  gsl_multiroot_function_fdf f = {&rosenbrock_f,
                                  &rosenbrock_df,
                                  &rosenbrock_fdf,
                                  n, &p};

  double x_init[2] = {-10.0, -5.0};
  gsl_vector *x = gsl_vector_alloc (n);

  gsl_vector_set (x, 0, x_init[0]);
  gsl_vector_set (x, 1, x_init[1]);

  T = gsl_multiroot_fdfsolver_gnewton;
  s = gsl_multiroot_fdfsolver_alloc (T, n);
  gsl_multiroot_fdfsolver_set (s, &f, x);

  print_state_d (iter, s);

  do
    {
      iter++;

      status = gsl_multiroot_fdfsolver_iterate (s);

      print_state_d (iter, s);

      if (status)
        break;

      status = gsl_multiroot_test_residual (s->f, 1e-7);
    }
  while (status == GSL_CONTINUE && iter < 1000);

  printf ("status = %s\n", gsl_strerror (status));

  gsl_multiroot_fdfsolver_free (s);
  gsl_vector_free (x);
  return 0;
}
*/
