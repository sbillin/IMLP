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
#ifndef _ParameterizedTest_H
#define _ParameterizedTest_H

#include <stdio.h>
#include <iostream>
#include <vector>

#include "cisstAlgorithmICP.h"
#include "cisstMesh.h"

//#include "cpd.h"

//#define MINIMIZE_OUTPUT_FILES
//#define DISABLE_OUTPUT_FILES

struct CustomTestParams
{
  virtual std::string toString() { return ""; };
};

struct TestParameters
{
  typedef void(*GenerateSamplesFuncType)(
    cisstMesh &mesh,
    TestParameters &params,
    unsigned int randSeed, unsigned int &randSeqPos,
    std::ifstream &randnStream,
    vctDynamicVector<vct3> &samples,
    vctDynamicVector<vct3> &sampleNorms,
    vctDynamicVector<vct3> &noisySamples,
    unsigned int trialNumber);
    //std::string *pOutputDir);

  // point to function that generates samples and produces sample noise
  GenerateSamplesFuncType GenerateSamplesFunc;

  std::string TestDescription;
  std::string AlgorithmDescription;
  std::string inputDir;
  std::string outputBaseDir;
  std::string outputDataDir;
  std::string outputCommonDir;
  std::string sourceMeshFile;
  std::string targetMeshFile;

  cisstICP::Options opt_ICP;
  //CPD::Options opt_CPD;

  cisstAlgorithmICP *pAlg;
  enum AlgType { AlgTypeNone, StdICP, IMLP, RobustICP }; // , CPD };
  AlgType algType;

  CustomTestParams *pCustomTestParams;

  bool bEnableIterationCallbacks;

  unsigned int numTrials;
  unsigned int nSamples;
  unsigned int nValidationSamples;

  double minOffsetAng, maxOffsetAng;
  double minOffsetPos, maxOffsetPos;
  unsigned int randSeed_Samples, randSeed_Xfms, randSeed_Validation, randSeed_MeshNoise;

  double noiseMeshSD;

  // default constructor
  TestParameters()
    : GenerateSamplesFunc(NULL),
    algType(AlgTypeNone),
    numTrials(0),
    nSamples(0),
    nValidationSamples(0),
    minOffsetAng(0.0),
    maxOffsetAng(0.0),
    minOffsetPos(0.0),
    maxOffsetPos(0.0),
    randSeed_Samples(0),
    randSeed_Xfms(0),
    randSeed_Validation(0),
    randSeed_MeshNoise(0),
    noiseMeshSD(0.0),
    bEnableIterationCallbacks(true)
  {};

  std::string toString() {
    std::stringstream ss;
    ss
      << "Test Params:" << std::endl
      << " TestType:  " << TestDescription << std::endl
      << " Algorithm: " << AlgorithmDescription << std::endl
      << " inputDir:        " << inputDir << std::endl
      << " outputBaseDir:   " << outputBaseDir << std::endl
      << " outputDataDir:   " << outputDataDir << std::endl
      << " outputCommonDir: " << outputCommonDir << std::endl
      << " sourceMeshFile:  " << sourceMeshFile << std::endl
      << " targetMeshFile:  " << targetMeshFile << std::endl
      << " algType: " << algType << std::endl
      << " numTrials: " << numTrials << std::endl
      << " nSamples: " << nSamples << std::endl
      << " nValidationSamples: " << nValidationSamples << std::endl
      << " randSeed_Samples: " << randSeed_Samples << std::endl
      << " randSeed_Xfms: " << randSeed_Xfms << std::endl
      << " randSeed_Validation: " << randSeed_Validation << std::endl
      << " randSeed_MeshNoise: " << randSeed_MeshNoise << std::endl
      << " Noise_MeshSD: " << noiseMeshSD << std::endl
      << " minOffsetAng: " << minOffsetAng << std::endl
      << " maxOffsetAng: " << maxOffsetAng << std::endl
      << " minOffsetPos: " << minOffsetPos << std::endl
      << " maxOffsetPos: " << maxOffsetPos << std::endl
      << " bEnableIterationCallbacks: " << bEnableIterationCallbacks << std::endl
      << pCustomTestParams->toString()
      << std::endl
      << opt_ICP.toString()
      //<< std::endl
      //<< opt_CPD.toString()
      ;
    return ss.str();
  }
};

void Run_ParameterizedTest(
  cisstMesh &mesh,
  cisstCovTreeBase *pTree,
  TestParameters params);

#endif 
