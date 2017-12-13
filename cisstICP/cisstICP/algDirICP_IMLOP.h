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
#ifndef _algDirICP_IMLOP_h
#define _algDirICP_IMLOP_h

#include "algDirICP.h"
#include "DirPDTree_Mesh.h"


// Debug Modes
//#define TEST_STD_ICP

class algDirICP_IMLOP : public algDirICP, public algDirPDTree
{
  //
  // This algorithm implements the von-Mises Fisher + Gaussian negative
  //   log likelihood cost function for position and orientation based
  //   registration. Derived classes implement different variations of
  //   this algorithm depending on whether the target shape is represented
  //   by a mesh or a point cloud.
  //
  //    -loglik[ C * exp( k*dot(Ny,Nx) - B*||Y-X||^2 ) ]
  //        where C  =  product of normalizations terms
  //                    for Fisher and 3D Gaussian distributions
  //              B  =  1/(2*sigma2)
  //
  //    cost    -k*Sum_i(dot(Ny,Nx)) + B*Sum_i(||Y-X||^2) - nSamples*log(C)
  //            
  //        NOTE:  This cost may be < 0 for very concentrated distn's since
  //                a pdf value may be > 1 (even though a probability value
  //                may not).
  //
  //  Covariance tree methods use the following modified cost:
  //
  //    cost    k*(1-dot(Ny,Nx)) + B*||Y-X||^2
  //    
  //        Note:  This modified cost is used so that cost is always >= 0.
  //               Note that the logC term does not affect the relative
  //               cost of one match vs another.
  //


  //--- Algorithm Parameters ---//

public:

  // noise model
  double k;       // Fisher dist'n
  double B;       // Gaussian dist'n
  double sigma2;  //  '' B = 1/(2*sigma2)

protected:

  double k_init, sigma2_init;
  double wRpos;   // weighting of positions for estimating k in each iteration

  bool   dynamicParamEst;  // if true, k & sigma2 are dynamically estimated from the match residuals
  double k_factor;         // k = k * k_factor; used to reduce the growth of k

  const double threshK;     // perfect match thresholds
  const double threshB;

  double R, Rnorm, Rpos;
  double circSD;


  // Match Filtering
  //
  //  Note: using buffers and buffer references for match filtering enables
  //        resizing the set of "good" points without altering memory allocation
  //

  unsigned int  nPosOutliers;       // # outliers by position
  unsigned int  nNormOutliers;      // # outliers by orientation that are not also position outliers

  //vctDynamicVector<double>  SquaredDistFromMean;

  // good samples (outliers removed)
  unsigned int  nGoodSamples;
  vctDynamicVector<vct3>        goodSamplePtsBuf;  
  vctDynamicVector<vct3>        goodMatchPtsBuf;
  vctDynamicVectorRef<vct3>     goodSamplePts;
  vctDynamicVectorRef<vct3>     goodMatchPts;

  vctDynamicVector<vct3>        goodSampleNormsBuf;
  vctDynamicVector<vct3>        goodMatchNormsBuf;
  vctDynamicVectorRef<vct3>     goodSampleNorms;
  vctDynamicVectorRef<vct3>     goodMatchNorms;

  //// summary error statistics
  //double gSumSqrDist_PostMatch;   // good sample distances
  //double oSumSqrDist_PostMatch;   // thresholded distances on outlier samples (helpful for smoothing cost function)
  //double gSumNormProducts_PostMatch;  // sum of norm products (good samples)
  //double oSumNormProducts_PostMatch;  // sum of norm products (outliers, thresholded)





  //--- Algorithm Methods ---//

public:

  // constructor
  algDirICP_IMLOP(
    DirPDTreeBase *pDirTree,
    vctDynamicVector<vct3> &samplePts,
    vctDynamicVector<vct3> &sampleNorms,
    double kinit = 0.0, double sigma2init = 1.0, double wRpos = 0.5,
    double kfactor = 1.0,
    bool dynamicallyEstParams = true)
    : algDirICP(pDirTree, samplePts, sampleNorms),
    algDirPDTree(pDirTree),
    k_init(kinit), sigma2_init(sigma2init), wRpos(wRpos),
    k_factor(kfactor),
    dynamicParamEst(dynamicallyEstParams),
    threshK(1.0e5), threshB(1.0e4)
  {
    // Ensure SetSamples function of this derived class gets called
    SetSamples(samplePts, sampleNorms);
  }

  // destructor
  virtual ~algDirICP_IMLOP() {}


  void    UpdateNoiseModel(double sumSqrDist, double sumNormProducts);
  double  ComputeRpos();

  void SetNoiseModel(
    double initK, double initSigma2, double wRpos, bool dynamicallyEstParams);

  void SetSamples(
    const vctDynamicVector<vct3> &argSamplePts,
    const vctDynamicVector<vct3> &argSampleNorms);


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
