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
#ifndef _alg2D_DirICP_h
#define _alg2D_DirICP_h

#include <cisstVector.h>
#include "DirPDTree2DBase.h"
#include "cisstException.h"

// debug
//#define SaveMatchesToFile
//#define ValidatePDTreeSearch


class alg2D_DirICP
{
  //
  // This is the base class for a family of ICP algorithms
  //

  //--- Standard Algorithm Parameters ---//

protected:

  DirPDTree2DBase  *pDirTree;   // target shape

public:

  unsigned int  nSamples;
  unsigned int  nOutliers;

  vctDynamicVector<vct2>  samplePts;      // source shape
  vctDynamicVector<vct2>  samplePtsXfmd;
  vctDynamicVector<vct2>  matchPts;

  vctDynamicVector<vct2>  sampleNorms;
  vctDynamicVector<vct2>  sampleNormsXfmd;
  vctDynamicVector<vct2>  matchNorms;

  vctDynamicVector<int>   matchDatums;
  vctDoubleVec            matchErrors;

  int   minNodesSearched, maxNodesSearched, avgNodesSearched;

  // These are merely variables for reporting run-time status to the terminal
  double errFuncNormWeight;
  double errFuncPosWeight;


  //// Match Error Statistics
  ////  NOTE: these are not computed by default
  ////        derived classes must call the appropriate match
  ////        error routines to compute these if desired
  ////

  //vctDynamicVector<vct2>  residuals_PostMatch;    // Pclosest - Psample
  //vctDoubleVec            sqrDist_PostMatch;      // ||Pclosest - Psample||^2
  //vctDoubleVec            dist_PostMatch;         // ||Pclosest - Psample||

  //vctDynamicVector<vct2>  residuals_PostRegister; // Pclosest - Psample
  //vctDoubleVec            sqrDist_PostRegister;   // ||Pclosest - Psample||^2
  //vctDoubleVec            dist_PostRegister;      // ||Pclosest - Psample||

  //vctDoubleVec  normProducts_PostMatch;     // dot(Nclosest,Nsamp)
  //vctDoubleVec  normProducts_PostRegister;

  //double  matchErrorAvg_PostMatch;
  //double  matchDistAvg_PostMatch;
  //double  sumSqrDist_PostMatch;

  //double  matchErrorAvg_PostRegister;
  //double  matchDistAvg_PostRegister;
  //double  sumSqrDist_PostRegister;

  //double  sumNormProducts_PostMatch;      // Sum_i( dot(Nclosest,Nsamp) )
  //double  sumNormProducts_PostRegister;


  //--- Standard Algorithm Methods ---//

public:

  // constructor
  alg2D_DirICP(
    DirPDTree2DBase *pDirTree,
    vctDynamicVector<vct2> &samplePts,
    vctDynamicVector<vct2> &sampleNorms);

  // destructor
  virtual ~alg2D_DirICP() {}

  void  ComputeCircErrorStatistics(double sumNormProducts, double &R, double &circSD);

protected:

  virtual void  SetSamples(
    vctDynamicVector<vct2> &argSamplePts, 
    vctDynamicVector<vct2> &argSampleNorms);

  virtual void  UpdateSampleXfmPositions(const vctFrm2 &F);

  virtual void  SamplePreMatch(unsigned int sampleIndex) {};
  virtual void  SamplePostMatch(unsigned int sampleIndex) {};

  //virtual void  ComputeErrors_PostMatch();
  //virtual void  ComputeErrors_PostRegister();


  //--- ICP Interface Methods ---//

public:

  // Initialize algorithm parameters before ICP begins
  virtual void    ICP_InitializeParameters(vctFrm2 &FGuess);

  // Update the algorithm parameters that depend on matching errors
  virtual void    ICP_UpdateParameters_PostMatch();

  // Update the algorithm parameters that depend on registration errors
  virtual void    ICP_UpdateParameters_PostRegister(vctFrm2 &Freg);

  virtual void    ICP_ComputeMatches();
  virtual unsigned int ICP_FilterMatches();   // handle outliers

  virtual vctFrm2 ICP_RegisterMatches() = 0;
  virtual double  ICP_EvaluateErrorFunction() = 0;

  //virtual std::vector<cisstICP::Callback> ICP_GetIterationCallbacks();

};

#endif
