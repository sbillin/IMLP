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

#include "alg2D_DirICP_StdICP.h"
#include "DirPDTree2DNode.h"
#include "cisstICP.h"

#include "RegisterP2P.h"


#define EPS  1e-12

double alg2D_DirICP_StdICP::ICP_EvaluateErrorFunction()
{
  // Cost Function = RMS (root mean square) error of the non-outlier matches

  vct2 residual;
  double sumSqrDist = 0.0;

  for (unsigned int s = 0; s < nGoodSamples; s++)
  {
    residual = goodSamplePts.Element(s) - goodMatchPts.Element(s);
    sumSqrDist += residual.NormSquare();
  }

  return sqrt(sumSqrDist / nSamples);

  ////
  ////  Note: add thresholded outlier error as well for
  ////        smoother cost function
  ////
  //double SqrErr = gSumSqrDist_PostMatch + oSumSqrDist_PostMatch;
  //return sqrt(SqrErr / nSamples);
}

vctFrm2 alg2D_DirICP_StdICP::ICP_RegisterMatches()
{
  CISST_THROW("Not Implemented");
}

void alg2D_DirICP_StdICP::ICP_InitializeParameters(vctFrm2 &FGuess)
{
  // initialize base class
  alg2D_DirICP::ICP_InitializeParameters(FGuess);

  // for reporting purposes
  errFuncPosWeight = 1.0;
  errFuncNormWeight = 0.0;
}

unsigned int alg2D_DirICP_StdICP::ICP_FilterMatches()
{

#if 1
  // Method 1: No Outlier Detection
  nGoodSamples = 0;
  nOutliers = 0;

  // use all samples
  for (unsigned int s = 0; s < nSamples; s++)
  { // copy all samples to buffer
    goodSamplePtsBuf.Element(s) = samplePts.Element(s);
    goodMatchPtsBuf.Element(s) = matchPts.Element(s);
    
    goodSampleNormsBuf.Element(s) = sampleNorms.Element(s);
    goodMatchNormsBuf.Element(s) = matchNorms.Element(s);
  }
  nGoodSamples = nSamples;

  // Non-destructively resize good sample reference vectors
  goodSamplePts.SetRef(goodSamplePtsBuf, 0, nGoodSamples);
  goodMatchPts.SetRef(goodMatchPtsBuf, 0, nGoodSamples);

  goodSampleNorms.SetRef(goodSampleNormsBuf, 0, nGoodSamples);
  goodMatchNorms.SetRef(goodMatchNormsBuf, 0, nGoodSamples);

  //gSumSqrDist_PostMatch = sumSqrDist_PostMatch;
  //oSumSqrDist_PostMatch = 0.0;

  return nOutliers;
#endif

}
