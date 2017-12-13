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
#ifndef _algDirICP_GIMLOP_h
#define _algDirICP_GIMLOP_h

#include "algDirICP.h"
#include "DirPDTree_Mesh.h"
#include "TriangleClosestPointSolver.h"

#include "Ellipsoid_OBB_Intersection_Solver.h"
#include "wrapper_dlib.h"

// gsl not compiled for 64-bit
//#include "wrapper_gsl.h"


//#define KENT_POS_ISOTROPIC
//#define SAVE_MATCHES
//#define ENABLE_DEBUG_FILE


class algDirICP_GIMLOP : public algDirICP, public algDirPDTree
{
  // This algorithm implements the Kent + Gaussian negative
  //   log likelihood cost function for position and orientation based
  //   registration. Derived classes implement different variations of
  //   this algorithm depending on whether the target shape is represented
  //   by a mesh or a point cloud.
  //
  //    -loglik[ ... ]
  //

protected:

#ifdef SAVE_MATCHES
  bool MatchIsotropic;
  std::string matchesFile, isoMatchesFile;
#endif

#ifdef ENABLE_DEBUG_FILE
  std::ofstream debugStream;
#endif

  //--- Algorithm Parameters ---//

public:

  // noise parameter estimation methods
  enum  PARAM_EST_TYPE{ PARAMS_FIXED, PARAMS_DYN_THRESH, PARAMS_FIXED_1stITER_POS };
  PARAM_EST_TYPE paramEstMethod;

  // for reporting purposes
  double meanSigma2, meanK, meanE, meanE_msmt;

protected:  

  Ellipsoid_OBB_Intersection_Solver IntersectionSolver;
  //wrapper_gsl gsl;
  wrapper_dlib dlib;

  // Kent parameters for each sample
  //  effective noise model
  vctDynamicVector<vct3x2> L;
  vctDoubleVec k;
  vctDoubleVec B;
  //  measurement noise model (static)
  vctDoubleVec k_msmt;
  vctDoubleVec E_msmt;  // eccentricity
  //  match uncertainty model (dynamic)
  double k_est;
  double meanK_msmt;

  // Gaussian parameters for each sample
  //  effective noise model
  vctDynamicVector<vct3x3> M;         // noise covariance of sample positions
  vctDynamicVector<vct3x3> invM;      // inverse noise covariance of sample positions
  vctDynamicVector<vct3x3> N;         // decomposition of inv(M) = N'N = R*D^2*R'
  vctDynamicVector<vct3x3> invN;
  vctDoubleVec Dmin;        // sqrt of smallest eigenvalue of inv(M)
  vctDoubleVec Emin;        // smallest eigenvalue of inv(M)
  //  measurement noise model (static)
  vctDynamicVector<vct3x3> M_msmt;         // noise covariance of sample positions
  vctDynamicVector<vct3x3> invM_msmt;      // inverse noise covariance of sample positions
  vctDynamicVector<vct3x3> N_msmt;         // decomposition of inv(M) = N'N = R*D^2*R'
  vctDynamicVector<vct3x3> invN_msmt;
  vctDoubleVec Dmin_msmt;   // sqrt of smallest eigenvalue of inv(M)
  vctDoubleVec Emin_msmt;   // smallest eigenvalue of inv(M)
  //  match uncertainty model (dynamic)
  double traceM_est;
  double meanTraceM_msmt;

  // noise parameters for xfmd samples
  vctDynamicVector<vct3x2> R_L;
  vctDynamicVector<vct3x3> R_invM_Rt;
  vctDynamicVector<vct3x3> N_Rt;
  vctDynamicVector<vct3x3> inv_N_Rt;   // inv(N*Rt) = R*invN

  // Variables for current sample point undergoing a match search
  //  these are set once for each sample point by the pre-match function
  //  and used whenever a tree search routine is called for that sample thereafter.
  //  These do not change relative to the node being searched => storing them here
  //  allows using the same PD tree search structure without adding new
  //  function arguments.
  vctFixedSizeMatrix<double, 3, 2> sampleR_L;
  double sampleK, sampleB;
  vct3x3 sampleR_InvM_Rt;
  vct3x3 sampleN_Rt;
  vct3x3 sample_inv_N_Rt;
  double sampleDmin;
  double sampleEmin;

  // Optimizer calculations common to both cost function and gradient
  vct6 x_prev;
  vct3 a, t;
  vctRot3 Ra;
  //vctDynamicVector<vct3> Yp_RaXp_t;
  vctDynamicVector<vct3> Yp_t;
  vctDynamicVector<vct3> Rat_Yp_RaXp_t;
  vctDynamicVector<vct3> invM_Rat_Yp_RaXp_t;
  vctDynamicVector<vct3> RaXn;
  vctDynamicVector<vct3x2> RaRL;


  //--- Algorithm Methods ---//

public:

  // constructor
  algDirICP_GIMLOP(
    DirPDTreeBase *pDirTree,
    vctDynamicVector<vct3> &samplePts,
    vctDynamicVector<vct3> &sampleNorms,
    const vctDynamicVector<double> &argK,
    const vctDynamicVector<double> &argE,
    const vctDynamicVector<vct3x2> &argL,
    const vctDynamicVector<vct3x3> &M,
    PARAM_EST_TYPE paramEst = PARAMS_FIXED);

  // destructor
  virtual ~algDirICP_GIMLOP() {}

  // dlib routines
  void    UpdateOptimizerCalculations(const vct6 &x);
  void    CostFunctionGradient(const vct6 &x, vct6 &g);
  double  CostFunctionValue(const vct6 &x);

  inline double MatchError(
    const vct3 &Xp, const vct3 &Xn,
    const vct3 &Yp, const vct3 &Yn,
    double k, double B, const vct3x2 &L,
    const vct3x3 &invM);

  void SetSamples(
    const vctDynamicVector<vct3> &samplePts,
    const vctDynamicVector<vct3> &sampleNorms,
    const vctDynamicVector<double> &argK,
    const vctDynamicVector<double> &argE,
    const vctDynamicVector<vct3x2> &argL,
    const vctDynamicVector<vct3x3> &argM,
    PARAM_EST_TYPE paramEst = PARAMS_FIXED);

  void SetNoiseModel(
    const vctDynamicVector<double> &argK,
    const vctDynamicVector<double> &argE,
    const vctDynamicVector<vct3x2> &argL,
    const vctDynamicVector<vct3x3> &M,
    PARAM_EST_TYPE paramEst);

  void UpdateNoiseModel_DynamicEstimates();
  void UpdateNoiseModel_SamplesXfmd(vctFrm3 &Freg);

  double  ComputeRpos();
  void    ComputeMatchUncertaintyEstimates();

  void ComputeCovDecompositions(
    const vctDynamicVector<vct3x3> &M, vctDynamicVector<vct3x3> &invM,
    vctDynamicVector<vct3x3> &N, vctDynamicVector<vct3x3> &invN,
    vctDoubleVec &Dmin, vctDoubleVec &Emin);
  void ComputeCovDecomposition_NonIter(
    const vct3x3 &M, vct3x3 &Minv, vct3x3 &N, vct3x3 &Ninv,
    double &Dmin, double &Emin);
  void ComputeCovDecomposition_SVD(
    const vct3x3 &M, vct3x3 &Minv, vct3x3 &N, vct3x3 &Ninv, 
    double &Dmin, double &Emin );

  // standard virtual routines
  virtual void SamplePreMatch(unsigned int sampleIndex);


  //--- ICP Interface Methods ---//

  void          ICP_InitializeParameters(vctFrm3 &FGuess);
  //void          ICP_UpdateParameters_PostMatch();
  void          ICP_UpdateParameters_PostRegister(vctFrm3 &Freg);
  vctFrm3       ICP_RegisterMatches();
  double        ICP_EvaluateErrorFunction();
  unsigned int  ICP_FilterMatches();


  //--- PD Tree Interface Methods ---//

  int  NodeMightBeCloser(
    const vct3 &v, const vct3 &n,
    DirPDTreeNode const *node,
    double ErrorBound);

  virtual double FindClosestPointOnDatum(const vct3 &v, const vct3 &n,
    vct3 &closest, vct3 &closestNorm,
    int datum) = 0;

  virtual int DatumMightBeCloser(const vct3 &v, const vct3 &n,
    int datum,
    double ErrorBound) = 0;

};
#endif
