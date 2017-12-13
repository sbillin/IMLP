// ****************************************************************************
//
//    Copyright (c) 2014, Seth Billings, Russell Taylor, Johns Hopkins University
//    All rights reserved.
//
//    Redistribution and use in source and binary forms, with or without
//    modification, are permitted provided that the following conditions are
//    met:
//
//    1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//    2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//
//    3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
//
//    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
//    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
//    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
//    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
//    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//  
// ****************************************************************************
#include "wrapper_dlib.h"
#include "algDirICP_GIMLOP.h"

#include <assert.h>
#undef NDEBUG       // enable debug in release mode

// define to print iteration details
#define DLIB_VERBOSE


//--- Globals ---//

namespace
{
  // Global variables
  //  (referenced from global functions)
  algDirICP_GIMLOP *Kent_dlib = NULL;

  // Global functions
  //  (needed for function pointers)
  double fValue( const wrapper_dlib::dlib_vector &x_dlib )
  {
    static vct6 x;
    x.Assign( x_dlib(0), x_dlib(1), x_dlib(2), 
              x_dlib(3), x_dlib(4), x_dlib(5) );
 
    return Kent_dlib->CostFunctionValue( x );
  }

  wrapper_dlib::dlib_vector fDerivative( const wrapper_dlib::dlib_vector &x_dlib )
  {
    static vct6   x;
    static vct6   g;
    static wrapper_dlib::dlib_vector  g_dlib(6);    // 6-element vector

    x.Assign( x_dlib(0), x_dlib(1), x_dlib(2), 
              x_dlib(3), x_dlib(4), x_dlib(5) );

    Kent_dlib->CostFunctionGradient( x, g );
    g_dlib(0) = g[0];
    g_dlib(1) = g[1];
    g_dlib(2) = g[2];
    g_dlib(3) = g[3];
    g_dlib(4) = g[4];
    g_dlib(5) = g[5];

    return g_dlib;
  }
} // local namespace



//--- Non-Globals ---//

// Constructor
wrapper_dlib::wrapper_dlib()  // algDirICP_GIMLOP *kent )
  : maxIter(40), //maxIter( 20 ),
  //tol_df( 1.0e-6 ),
  gradientNormThresh( 1.0e-3 )
{
  Kent_dlib = NULL;     // initialize global variable
  //Kent_dlib = kent;  // initialize global variable
}


// passing a pointer to the algorithm is necessary if this library is being
//  used with multiple Kent algorithms simultaneously (even if single threaded)
vct6 wrapper_dlib::ComputeRegistration(const vct6 &x0, algDirICP_GIMLOP *kent)
{
  // initialize global pointer to algorithm
  Kent_dlib = kent;
  
  dlib_vector x_dlib(6);  // 6-element vector

  //double f0 = Kent_gsl->CostFunctionValue( x0 );
  //double thresh_df = tol_df * f0;
  //double thresh_df = tol_df * Kent_gsl->k_sum;
  //std::cout << "thresh_df: " << thresh_df << std::endl;

  try
  {
    // Now we use the find_min() function to find the minimum point.  The first argument
    // to this routine is the search strategy we want to use.  The second argument is the 
    // stopping strategy.
    //   objective_delta_stop_strategy:  stop if df < threshold

    // The other arguments to find_min() are the function to be minimized, its derivative, 
    // then the starting point, and the last is an acceptable minimum value.  
    // That is, if the algorithm finds any inputs that give an output value less than
    // this then it will stop immediately.  Usually you supply a number smaller than 
    // the actual global minimum.

    x_dlib(0) = x0[0];
    x_dlib(1) = x0[1];
    x_dlib(2) = x0[2];
    x_dlib(3) = x0[3];
    x_dlib(4) = x0[4];
    x_dlib(5) = x0[5];

    //std::cout << "Difference between analytic derivative and numerical approximation of derivative: " 
    //  << dlib::length(dlib::derivative(fValue)(x_dlib) - fDerivative(x_dlib)) << std::endl;


    dlib::find_min( dlib::bfgs_search_strategy(),
#ifdef DLIB_VERBOSE
      //dlib::objective_delta_stop_strategy( Tol_df,maxIter ).be_verbose(),
      dlib::gradient_norm_stop_strategy( gradientNormThresh,maxIter ).be_verbose(),
#else
      //dlib::objective_delta_stop_strategy( Tol_df,maxIter ),
      dlib::gradient_norm_stop_strategy( gradientNormThresh,maxIter ),
#endif
      fValue, fDerivative,
      x_dlib, -1.0);

    //find_min_using_approximate_derivatives( bfgs_search_strategy(),
    //                                        objective_delta_stop_strategy(Tol_df).be_verbose(),
    //                                        dlib_fValue, 
    //                                        x0_dlib, -1.0);
  }
  catch (std::exception& e)
  {
    std::cout << "DLIB EXCEPTION: " << e.what() << std::endl;
    assert(0);
  }

  Kent_dlib = NULL;

  return vct6( x_dlib(0), x_dlib(1), x_dlib(2), 
               x_dlib(3), x_dlib(4), x_dlib(5) );
}