#include "wrapper_gsl.h"
#include "cisstICPNormalsAlgorithm_Kent.h"
#include <stdio.h>

#include <assert.h>
#undef NDEBUG       // enable debug in release mode

// define to print iteration details
//#define GSL_VERBOSE


// Global variables
//  (referenced from global functions)
//  suffix prevents collision with other wrappers
cisstICPNormalsAlgorithm_Kent *Kent_gsl = NULL;


double kent_f( const gsl_vector *x, void *params )
{
  static vct6 vctx;
  vctx.Assign( 
    gsl_vector_get(x, 0),
    gsl_vector_get(x, 1),
    gsl_vector_get(x, 2),
    gsl_vector_get(x, 3),
    gsl_vector_get(x, 4),
    gsl_vector_get(x, 5)
  );
 
  return Kent_gsl->CostFunctionValue( vctx );
}

void kent_df( const gsl_vector *x, void *params, gsl_vector *g )
{
  static vct6   vctx;
  static vct6   vctg;

  vctx.Assign( 
    gsl_vector_get(x, 0),
    gsl_vector_get(x, 1),
    gsl_vector_get(x, 2),
    gsl_vector_get(x, 3),
    gsl_vector_get(x, 4),
    gsl_vector_get(x, 5)
  );

  Kent_gsl->CostFunctionGradient( vctx, vctg );

  gsl_vector_set(g, 0, vctg[0]);
  gsl_vector_set(g, 1, vctg[1]);
  gsl_vector_set(g, 2, vctg[2]);
  gsl_vector_set(g, 3, vctg[3]);
  gsl_vector_set(g, 4, vctg[4]);
  gsl_vector_set(g, 5, vctg[5]);
}

void kent_fdf( const gsl_vector *x, void *params, double *f, gsl_vector *g )
{
  *f = kent_f( x, params ); 
  kent_df( x, params, g );
}



// Constructor
wrapper_gsl::wrapper_gsl( cisstICPNormalsAlgorithm_Kent *kent )
: T(NULL),
  min(NULL),
  maxIter( 50 ),
  //tol_df( 1.0e-6 ),
  thresh_gradient( 1.0e-3 )
{
  Kent_gsl = kent;  // initialize global variable

  double lineSearchTol = 0.1;

  // gsl minimizer
  T = gsl_multimin_fdfminimizer_vector_bfgs2;
  min = gsl_multimin_fdfminimizer_alloc(T, 6);
  if (min == NULL)
  {
    std::cout << "ERROR: create GSL minimizer failed" << std::endl;
    assert(0);
  }
  
  // gsl function
  func.n = 6;
  func.params = NULL;
  func.f = &kent_f;
  func.df = &kent_df;
  func.fdf = &kent_fdf;

  // gsl vector
  x_gsl = gsl_vector_alloc(6);
}

// Destructor
wrapper_gsl::~wrapper_gsl( )
{
  gsl_multimin_fdfminimizer_free( min );
  gsl_vector_free( x_gsl );
}

vct6 wrapper_gsl::ComputeRegistration( const vct6 &x0 )
{
  int status;
  unsigned int iter;
  double initStepSize, lineSearchTol;

  gsl_vector_set(x_gsl, 0, x0[0]);
  gsl_vector_set(x_gsl, 1, x0[1]);
  gsl_vector_set(x_gsl, 2, x0[2]);
  gsl_vector_set(x_gsl, 3, x0[3]);
  gsl_vector_set(x_gsl, 4, x0[4]);
  gsl_vector_set(x_gsl, 5, x0[5]);

  initStepSize = 0.1;
  lineSearchTol = 0.1;
  iter = 0;

  gsl_multimin_fdfminimizer_set( min, &func, x_gsl, initStepSize, lineSearchTol );

  //double df;
  //double f0 = Kent_gsl->CostFunctionValue( x0 );
  //double thresh_df = tol_df * f0;
  //double thresh_df = tol_df * Kent_gsl->k_sum;
  //std::cout << "thresh_df: " << thresh_df << std::endl;
  //double fprev = f0;

#ifdef GSL_VERBOSE
  printf( "iter 0: f: %4.5f  x: %.3f %.3f %.3f %.3f %.3f %.3f\n", 
          f0,
          gsl_vector_get(x_gsl, 0), 
          gsl_vector_get(x_gsl, 1), 
          gsl_vector_get(x_gsl, 2), 
          gsl_vector_get(x_gsl, 3), 
          gsl_vector_get(x_gsl, 4), 
          gsl_vector_get(x_gsl, 5)
          );
#endif
  do
  {
    iter++;
    status = gsl_multimin_fdfminimizer_iterate( min );
    if (status) { break; }  // unable to improve on current estimate

    //df = fprev - min->f;

    status = gsl_multimin_test_gradient( min->gradient, thresh_gradient );
    //if (status == GSL_SUCCESS) {printf ("Minimum found");}

#ifdef GSL_VERBOSE
    printf( "iter %u: f: %4.5f  x: %.3f %.3f %.3f %.3f %.3f %.3f\n", 
            iter,
            min->f,
            gsl_vector_get(min->x, 0), 
            gsl_vector_get(min->x, 1), 
            gsl_vector_get(min->x, 2), 
            gsl_vector_get(min->x, 3), 
            gsl_vector_get(min->x, 4), 
            gsl_vector_get(min->x, 5)
            );
#endif
  }  
  while ( status == GSL_CONTINUE && iter < maxIter );
  //while ( df > thresh_df && iter < maxIter);

  return vct6( 
    gsl_vector_get(min->x, 0),
    gsl_vector_get(min->x, 1),
    gsl_vector_get(min->x, 2),
    gsl_vector_get(min->x, 3),
    gsl_vector_get(min->x, 4),
    gsl_vector_get(min->x, 5)
  );    

  // This function resets the minimizer s to use the current point as a new starting point. 
  //
  //  Function: int   gsl_multimin_fdfminimizer_restart (gsl_multimin_fdfminimizer * s)
  //
  // The minimizer maintains a current best estimate of the minimum at all times. 
  // This information can be accessed with the following auxiliary functions:
  //
  //  Function: gsl_vector *  gsl_multimin_fdfminimizer_x (const gsl_multimin_fdfminimizer * s)  
  //  Function: double        gsl_multimin_fdfminimizer_minimum (const gsl_multimin_fdfminimizer * s)
  //  Function: gsl_vector *  gsl_multimin_fdfminimizer_gradient (const gsl_multimin_fdfminimizer * s)
  //
  // Termination Conditions:
  //  Function: int   gsl_multimin_test_gradient (const gsl_vector * g, double epsabs)
  //    Returns: GSL_SUCCESS or GSL_CONTINUE
}