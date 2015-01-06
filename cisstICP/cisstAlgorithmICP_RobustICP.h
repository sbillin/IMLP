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
#ifndef _cisstAlgorithmICP_RobustICP_h
#define _cisstAlgorithmICP_RobustICP_h

#include "cisstAlgorithmICP_StdICP.h"


class cisstAlgorithmICP_RobustICP : public cisstAlgorithmICP_StdICP //, public cisstAlgorithmCovTree_CP
{
  //
  // This is a base class for a robust ICP algorithm, adding outlier detection to
  //  standard ICP
  //
  //  Implements the outlier detection technique described in:
  //   Zhang, "Ierative Point Matching for Registration of Free-Form Curves and Surfaces", IJCV 1994
  //


  //--- Algorithm Parameters ---//

public:

  double D, D0max, DImax;
  double epsilon;
  double distAvg, distSD;

  bool bFirstIter_Matches;

  vctDynamicVector<double> matchDist;

  // for Round 1 of match filtering
  unsigned int nFilteredSamples;
  vctDynamicVectorRef<vct3>     filterSamplePts;
  vctDynamicVectorRef<vct3>     filterMatchPts;
  vctDynamicVectorRef<double>   filterMatchDist;

  vctDynamicVector<vct3>        filterSamplePtsBuf;
  vctDynamicVector<vct3>        filterMatchPtsBuf;
  vctDynamicVector<double>      filterMatchDistBuf;


  //--- Algorithm Methods ---//

public:

  // constructor
  cisstAlgorithmICP_RobustICP(
    cisstCovTreeBase *pTree, 
    vctDynamicVector<vct3> &samplePts,
    double D, double D0max)
    : cisstAlgorithmICP_StdICP(pTree, samplePts),
    D(D),
    D0max(D0max)
  {}

  double ComputeEpsilon(vctDynamicVector<double> &sampleDist);
  void printHistogram(
    vctDynamicVector<unsigned int> bins,
    unsigned int peakBin, unsigned int valleyBin,
    double minDist, double maxDist,
    double binWidth);

  //--- ICP Interface Methods ---//

  void          ICP_InitializeParameters(vctFrm3 &FGuess);
  void          ICP_UpdateParameters_PostMatch();
  unsigned int  ICP_FilterMatches();

  virtual std::vector<cisstICP::Callback> ICP_GetIterationCallbacks();

  //void    ICP_RegisterMatches(vctFrm3 &Freg);
  //double  ICP_EvaluateErrorFunction();
  //void  ICP_UpdateParameters_PostRegister(vctFrm3 &Freg);
  //void  ICP_ComputeMatches();

};
#endif
