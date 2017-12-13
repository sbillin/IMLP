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
#ifndef _algDirICP_PIMLOP_h
#define _algDirICP_PIMLOP_h

#include "algDirICP.h"
#include "DirPDTree_Mesh.h"
#include "algDirPDTree_vonMisesPrj.h"
#include "algDirICP_PIMLOP_dlibWrapper.h"
//#include "algPDTree_MLP.h"


class algDirICP_PIMLOP : public algDirICP //, public algDirPDTree_vonMisesPrj
{
  //
  // This algorithm implements the projected von-Mises + generalized Gaussian
  // noise model for ICP.
  //
  // Derived classes implement different variations of
  // this algorithm depending on whether the target shape is represented
  // by a mesh or a point cloud.
  //

  //--- Algorithm Parameters ---//

public:

  //std::ofstream debugStream;
  //Ellipsoid_OBB_Intersection_Solver IntersectionSolver;
  
  algDirICP_PIMLOP_dlibWrapper dlib;

  bool dynamicParamEst;

  
  vctDynamicVector<vct2>    sampleNorms2d;  // 2D (in-plane) sample orientations (Xpln)
  vctDynamicVector<vctRot3> Rx_pln;         // xfms in-plane orientation to x-local 3d space
                                            //  i.e. sampleNorm3d = Rx_pln * [Xpln; 0]
  vctDynamicVector<vctRot3> Ry_pln;         // Ry_pln[i] = Rreg * Rx_pln[i]

  //vctDynamicVector<vct3> matchNorms3d;    // Y3d = Ry_pln' * Yn;
  //vctDynamicVector<vct2> matchNorms2d;    // match orientations transformed and projected to the 
  //                                        //  local 2d planar coordinates of each sample orientation


  // von-Mises noise-model params for each sample
  //  NOTE: points having no orientation data may be used by
  //        setting k = 0 for those points
  vctDynamicVector<double> k;
  //double k_sum;

  // Gaussian noise-model params for non-transformed samples
  //  NOTE: M = Mx since target is assumed to have zero error
  vctDynamicVector<vct3x3> M;         // noise covariance of sample positions
  vctDynamicVector<vct3x3> invM;      // inverse noise covariance of sample positions
  vctDynamicVector<vct3x3> N;         // decomposition of inv(M) = N'N = R*D^2*R'
  vctDynamicVector<vct3x3> invN;
  vctDynamicVector<double> Dmin;      // sqrt of smallest eigenvalue of inv(M)
  //vctDynamicVector<double> Emin;      // smallest eigenvalue of inv(M)

  // noise parameters for xfmd samples
  vctDynamicVector<vct3x3> R_invM_Rt;
  vctDynamicVector<vct3x3> N_Rt;
  vctDynamicVector<vct3x3> inv_N_Rt;   // inv(N*Rt) = R*invN

  // Optimizer calculations common to both cost function and gradient
  vct6 x_prev;
  vct3 a, t;
  vctRot3 Ra;
  //vctDynamicVector<vct3> Yp_RaXp_t;
  vctDynamicVector<vct3> Yp_t;
  vctDynamicVector<vct3> Rat_Yp_RaXp_t;
  vctDynamicVector<vct3> invM_Rat_Yp_RaXp_t;
  vctDynamicVector<vct2> Yprj;  
  vctDynamicVector<vct2> Ypln;
  vctDynamicVector<double> Ynorm;
  //vctDynamicVector<vct2> RtRatYn;
  //vctDynamicVector<vct3> RaXn;  
  //vctDynamicVector<vctFixedSizeMatrix<double, 3, 2>> RaRL;


  //--- Algorithm Methods ---//

public:

  // constructor
  algDirICP_PIMLOP(
    DirPDTreeBase *pDirTree,
    const vctDynamicVector<vct3> &samplePts,
    const vctDynamicVector<vct2> &sampleNorms2d,
    const vctDynamicVector<vctRot3> &Rx_pln,
    const vctDynamicVector<double> &sample_k,
    const vctDynamicVector<vct3x3> &sample_M);

  // destructor
  virtual ~algDirICP_PIMLOP() {}

  void SetSamples(
    const vctDynamicVector<vct3> &samplePts,
    const vctDynamicVector<vct2> &sampleNorms2d,
    const vctDynamicVector<vctRot3> &Rx_pln,
    const vctDynamicVector<double> &sample_k,
    const vctDynamicVector<vct3x3> &sample_M);

  //void UpdateSampleXfmPositions(const vctFrm3 &F);
  void UpdateNoiseModel_SamplesXfmd(vctFrm3 &Freg);
  void ComputeNoiseModelDecompositions();

  void    UpdateOptimizerCalculations(const vct6 &x);
  void    CostFunctionGradient(const vct6 &x, vct6 &g);
  double  CostFunctionValue(const vct6 &x);

  // Xp    ~ transformed source point
  // Xpln  ~ non-transformed source orientation in 2d plane coordinates
  //         where Xn (3d orientation in local coors) = Rx_pln * [Xpln; 0]
  // Yp    ~ target point
  // Yn    ~ target 3d orientation
  // Ry_pln  ~ Yprj = Pxy( Ry_pln' * Yn)  (projection of Yn to plane containing Xpln)
  // k     ~ orientation concentration
  // invM  ~ positional covariance for transformed sample
  inline double MatchError(
    const vct3 &Xp, const vct2 &Xpln,
    const vct3 &Yp, const vct3 &Yn,
    const vctRot3 &Ry_pln,
    double k, const vct3x3 &invM);

  // standard ICP algorithm virtual routines
  virtual void SamplePreMatch(unsigned int sampleIndex) = 0;


  //--- ICP Interface Methods ---//

  void          ICP_InitializeParameters(vctFrm3 &FGuess);
  //void          ICP_UpdateParameters_PostMatch();
  void          ICP_UpdateParameters_PostRegister(vctFrm3 &Freg);
  vctFrm3       ICP_RegisterMatches();
  double        ICP_EvaluateErrorFunction();
  unsigned int  ICP_FilterMatches();

};
#endif
