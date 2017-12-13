#ifndef _wrapper_gsl_h
#define _wrapper_gsl_h

#include "cisstVector.h"
#include "./gsl/gsl_multimin.h"
//#include "cisstICPNormalsAlgorithm_Kent.h"

class cisstICPNormalsAlgorithm_Kent;  // forward decleration

class wrapper_gsl
{
public:

  // -- Variables -- //

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *min;
  gsl_multimin_function_fdf func;

  gsl_vector *x_gsl;

  double tol_df, thresh_gradient;
  unsigned int maxIter;


  // -- Methods -- //

  // constructor
  wrapper_gsl( cisstICPNormalsAlgorithm_Kent *kent );
  ~wrapper_gsl();

  vct6  ComputeRegistration( const vct6 &x0 );
};

#endif