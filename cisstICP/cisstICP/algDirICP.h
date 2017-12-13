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
#ifndef _algDirICP_h
#define _algDirICP_h

#include <cisstVector.h>
#include "algICP.h"
#include "DirPDTreeBase.h"

// debug
//#define SaveMatchesToFile
//#define ValidatePDTreeSearch


class algDirICP : public algICP
{
  //
  // This is the base class for a family of ICP algorithms
  //

  //--- Standard Algorithm Parameters ---//

public:

  DirPDTreeBase  *pDirTree;   // target shape

  vctDynamicVector<vct3>  sampleNorms;    // source shape
  vctDynamicVector<vct3>  sampleNormsXfmd;
  vctDynamicVector<vct3>  matchNorms;

  // these variables are for displaying run-time status on the terminal
  double errFuncNormWeight;
  double errFuncPosWeight;


  // Match Error Statistics

  //vctDoubleVec  normProducts_PostMatch;     // dot(Nclosest,Nsamp)
  //vctDoubleVec  normProducts_PostRegister;

  //double  sumNormProducts_PostMatch;      // Sum_i( dot(Nclosest,Nsamp) )
  //double  sumNormProducts_PostRegister;

  //double  dThetaMin, dThetaMax, dThetaAvg;

  //double  R_PostMatch;          // mean of norm products
  //double  circSD_PostMatch;     // circular standard deviation of theta (angular match error, radians)  
  //double  R_PostReg;
  //double  circSD_PostReg;


  //--- Standard Algorithm Methods ---//

public:

  // constructor
  algDirICP(
    DirPDTreeBase *pDirTree,
    const vctDynamicVector<vct3> &samplePts,
    const vctDynamicVector<vct3> &sampleNorms);

  // destructor
  virtual ~algDirICP() {}

  // R       ~  mean of norm products
  // circSD  ~  circular standard deviation of theta (angular match error in radians)
  void  ComputeCircErrorStatistics(double sumNormProducts, double &R, double &circSD);

  virtual void ComputeMatchStatistics(
    double &PosAvg, double &PosStdDev,
    double &AngAvg, double &AngStdDev);

  virtual void  SetSamples(
    const vctDynamicVector<vct3> &argSamplePts,
    const vctDynamicVector<vct3> &argSampleNorms);

protected:

  virtual void  UpdateSampleXfmPositions(const vctFrm3 &F);

  virtual void  SamplePreMatch(unsigned int sampleIndex) {};
  virtual void  SamplePostMatch(unsigned int sampleIndex) {};

  //virtual void  ComputeErrors_PostMatch();
  //virtual void  ComputeErrors_PostRegister();


  //--- ICP Interface Methods ---//

public:

  // Initialize algorithm parameters before ICP begins
  virtual void  ICP_InitializeParameters(vctFrm3 &FGuess);

  // Update the algorithm parameters that depend on matching errors
  virtual void  ICP_UpdateParameters_PostMatch();

  // Update the algorithm parameters that depend on registration errors
  virtual void  ICP_UpdateParameters_PostRegister(vctFrm3 &Freg);

  virtual void    ICP_ComputeMatches();
  virtual unsigned int ICP_FilterMatches();   // handle outliers

  virtual vctFrm3 ICP_RegisterMatches() = 0;
  virtual double  ICP_EvaluateErrorFunction() = 0;

  virtual std::vector<cisstICP::Callback> ICP_GetIterationCallbacks();

};

#endif
