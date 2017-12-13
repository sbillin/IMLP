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

#ifndef _algICP_h
#define _algICP_h

#include <cisstVector.h>
#include "PDTreeBase.h"
#include "cisstICP.h"   // for callbacks

// debug
//#define SaveMatchesToFile
//#define ValidatePDTreeSearch

#ifdef ValidatePDTreeSearch
//#define ValidateByEuclideanDist
extern std::ofstream validFS;
extern vct3 validPoint;
extern vct3 validNorm;
extern int  validDatum;
extern double validDist;
extern double validAng;
extern double validError;
extern double searchError;
extern unsigned int numValidDatums;
extern unsigned int numInvalidDatums;
extern double validPercent;
extern double doubleEps;
extern int validIter;
#endif

#ifdef SaveMatchesToFile
extern std::string saveMatchesDir;
extern int saveMatchesIter;
#endif


class algICP
{
  //
  // This is the base class for a family of ICP algorithms
  //

  //--- Standard Algorithm Parameters ---//

public:

  PDTreeBase  *pTree;   // target shape

  unsigned int  nSamples;
  unsigned int  nOutliers;

  vctDynamicVector<vct3>  samplePts;    // source shape
  vctDynamicVector<vct3>  samplePtsXfmd;
  vctDynamicVector<vct3>  matchPts;
  vctDynamicVector<int>   matchDatums;
  vctDoubleVec            matchErrors;

  int minNodesSearched, maxNodesSearched, avgNodesSearched;

  // current registration
  vctFrm3 Freg;

  //// TODO: put these only in the derived classes that actually use them
  //// Match Error Statistics
  ////  NOTE: these are not computed by default
  ////        derived classes must call the appropriate match
  ////        error routines to compute these if desired
  ////

  //vctDynamicVector<vct3>  residuals_PostMatch;    // Pclosest - Psample
  //vctDoubleVec            sqrDist_PostMatch;      // ||Pclosest - Psample||^2
  //vctDoubleVec            dist_PostMatch;         // ||Pclosest - Psample||

  //vctDynamicVector<vct3>  residuals_PostRegister; // Pclosest - Psample
  //vctDoubleVec            sqrDist_PostRegister;   // ||Pclosest - Psample||^2
  //vctDoubleVec            dist_PostRegister;      // ||Pclosest - Psample||

  ////double  matchErrorAvg_PostMatch;
  //double  matchDistAvg_PostMatch;
  //double  sumSqrDist_PostMatch;

  ////double  matchErrorAvg_PostRegister;
  //double  matchDistAvg_PostRegister;
  //double  sumSqrDist_PostRegister;


  //--- Standard Algorithm Methods ---//

public:

  // constructor
  algICP(PDTreeBase *pTree, const vctDynamicVector<vct3> &samplePts);

  // destructor
  virtual ~algICP() {}

  virtual void  SetSamples(const vctDynamicVector<vct3> &argSamplePts);

  virtual void ComputeMatchStatistics(double &Avg, double &StdDev);

protected:

  virtual void  UpdateSampleXfmPositions(const vctFrm3 &F);

  virtual void  SamplePreMatch(unsigned int sampleIndex) {};
  virtual void  SamplePostMatch(unsigned int sampleIndex) {};

  //// these are convenience routines that are not called by the base algorithm
  //virtual void  ComputeErrors_PostMatch();
  //virtual void  ComputeErrors_PostRegister();


  //--- ICP Interface Methods ---//

public:

  // Initialize algorithm parameters before ICP begins
  virtual void    ICP_InitializeParameters( vctFrm3 &FGuess );

  // Update the algorithm parameters that depend on matching errors
  virtual void    ICP_UpdateParameters_PostMatch();

  // Update the algorithm parameters that depend on registration errors
  virtual void    ICP_UpdateParameters_PostRegister( vctFrm3 &Freg );

  virtual void    ICP_ComputeMatches();
  virtual unsigned int ICP_FilterMatches();   // handle outliers

  virtual vctFrm3 ICP_RegisterMatches() = 0;
  virtual double  ICP_EvaluateErrorFunction() = 0;

  // enables an algorithm to request termination based on
  //  algorithm-specific criteria
  virtual bool    ICP_Terminate( vctFrm3 &Freg ) { return false; }

  virtual std::vector<cisstICP::Callback> ICP_GetIterationCallbacks();

};

#endif
