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
#ifndef _algICP_StdICP_h
#define _algICP_StdICP_h

#include "algICP.h"


class algICP_StdICP 
  : public algICP //, public algPDTree_CP
{
  //
  // This is a base class for the standard ICP algorithm
  //
  //  Minimize Cost:  Sum_i( ||Yi,T(Xi)||^2 )
  //
  //    T:   Rigid body transform [R,t] wrt which the cost function is optimized
  //         (i.e. the values we wish to solve)
  //    Xi:  sample point i
  //    Yi:  model point closest to sample Xi
  //
  // --------------------------
  //
  //  -LogLikelihood function:   -loglik[ C * exp( B*||Y-T(X)||^2 ) ]
  //       where C  =  normalization term for 3D Gaussian distributions
  //             B  =  1/(2*sigma2)
  //
  //   =>  cost function:    Sum_i(||Yi-T(Xi)||^2)
  //


  //--- Algorithm Parameters ---//

public:

  //
  //  Note: using buffers and buffer references for match filtering enables
  //        resizing the set of "good" points without altering memory allocation
  //

  // good samples (outliers removed)
  unsigned int  nGoodSamples;

  vctDynamicVectorRef<vct3>     goodSamplePts;
  vctDynamicVectorRef<vct3>     goodMatchPts;

  vctDynamicVector<vct3>        goodSamplePtsBuf;
  vctDynamicVector<vct3>        goodMatchPtsBuf;


  //// error statistics
  //double gSumSqrDist_PostMatch;   // good sample distances
  //double oSumSqrDist_PostMatch;   // thresholded distances on outlier samples (helpful for smoothing cost function)


  //--- Algorithm Methods ---//

public:

  // constructor
  algICP_StdICP(PDTreeBase *pTree, vctDynamicVector<vct3> &samplePts)
    : algICP(pTree, samplePts)
  {
    SetSamples(samplePts);
  }  

  // destructor
  virtual ~algICP_StdICP() {}

  void SetSamples(vctDynamicVector<vct3> &samplePts)
  {
    algICP::SetSamples(samplePts);

    // size sample buffers
    goodSamplePtsBuf.SetSize(nSamples);
    goodMatchPtsBuf.SetSize(nSamples);
  }

  virtual void  ComputeMatchStatistics(double &Avg, double &StdDev);


  //--- ICP Interface Methods ---//

  void          ICP_InitializeParameters(vctFrm3 &FGuess);
  vctFrm3       ICP_RegisterMatches();
  double        ICP_EvaluateErrorFunction();
  unsigned int  ICP_FilterMatches();

  //void  ICP_UpdateParameters_PostMatch();
  //void  ICP_UpdateParameters_PostRegister(vctFrm3 &Freg);
  //void  ICP_ComputeMatches();
  //std::vector<cisstICP::Callback> ICP_GetIterationCallbacks();

};
#endif
