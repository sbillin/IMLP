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
#ifndef _ParameterizedTest_PointCloud_PartialOverlap_H
#define _ParameterizedTest_PointCloud_PartialOverlap_H

#include <stdio.h>
#include <iostream>
#include <vector>

#include <windows.h>

#include <cisstOSAbstraction.h>

#include "cisstMesh.h"
#include "cisstPointCloud.h"

#include "cisstCovTree_PointCloud.h"
#include "cisstAlgorithmICP_IMLP_PointCloud.h"
#include "cisstAlgorithmICP_StdICP_PointCloud.h"
#include "cisstAlgorithmICP_RobustICP_PointCloud.h"
#include "cisstAlgorithmICP_IMLP_MahalDist_PointCloud.h"
#include "cisstAlgorithmICP_IMLP_ClosestPoint_PointCloud.h"

#include "ParameterizedTest.h"
#include "utility.h"
#include <limits.h>

// declerations: test routines
void Run_ParameterizedTest_PointCloud_PartialOverlap(TestParameters params);

// declerations: other
void SetOutputDir(TestParameters &params, std::string outputCommonDir, std::string outputSubFolder);
void GenerateSamplePointSet_PointCloud_PartialOverlap(
  cisstMesh &mesh,
  TestParameters &params,
  unsigned int randSeed, unsigned int &randSeqPos,
  std::ifstream &randnStream,
  vctDynamicVector<vct3> &samples,
  vctDynamicVector<vct3> &sampleNorms,
  vctDynamicVector<vct3> &noisySamples,
  unsigned int trialNum);

struct CustomTestParams_PointCloud_PartialOverlap : public CustomTestParams
{
  enum NoiseType{ NoiseInvalid = 0, NoiseSurface, NoiseGlobal, NoiseRandom, NoiseNone };

  //unsigned int nGoodSamples;
  //double percentOutliers;
  //double minPosOffsetOutlier, maxPosOffsetOutlier;

  //double        sampleNoise_InPlaneSD, sampleNoise_PerpPlaneSD;
  vct3          sampleNoise_EigValSqrt;
  NoiseType     sampleNoise_Type;
  //unsigned int  sampleNoise_CovRandSeed;

  // IMLP Params
  double outlier_ChiSquareThreshold;
  double sigma2Max;
  double surfaceModel_InPlaneSD, surfaceModel_PerpPlaneSD;

  bool bSampleCov_ApplyNoiseModel;
  bool bSampleCov_ApplySurfaceModel;
  bool bTargetCov_ApplySurfaceModel;
  bool bTargetCov_DynamicSurfaceModel;

  // Robust ICP Params
  //double D, D0max;

  // default constructor
  CustomTestParams_PointCloud_PartialOverlap() :
    //nGoodSamples(0),
    //percentOutliers(0.0),
    //minPosOffsetOutlier(0.0),
    //maxPosOffsetOutlier(0.0),
    //sampleNoise_InPlaneSD(0.0),
    //sampleNoise_PerpPlaneSD(0.0),
    sampleNoise_EigValSqrt(vct3(0.0)),
    sampleNoise_Type(NoiseInvalid),
    //sampleNoise_CovRandSeed(0),
    // IMLP
    outlier_ChiSquareThreshold(7.81),
    sigma2Max(std::numeric_limits<double>::max()),
    surfaceModel_InPlaneSD(0.0),
    surfaceModel_PerpPlaneSD(0.0),
    bSampleCov_ApplyNoiseModel(true),
    bSampleCov_ApplySurfaceModel(true),
    bTargetCov_ApplySurfaceModel(true),
    bTargetCov_DynamicSurfaceModel(false)
  {};

  std::string toString() {
    std::stringstream ss;
    ss
      << std::endl
      << "Custom Params: (General)" << std::endl
      //<< " nGoodSamples: " << nGoodSamples << std::endl
      //<< " percentOutliers: " << percentOutliers << std::endl
      //<< " minPosOffsetOutlier: " << minPosOffsetOutlier << std::endl
      //<< " maxPosOffsetOutlier: " << maxPosOffsetOutlier << std::endl
      //<< " sampleNoise_InPlaneSD:    " << sampleNoise_InPlaneSD << std::endl
      //<< " sampleNoise_PerpPlaneSD:  " << sampleNoise_PerpPlaneSD << std::endl
      << " sampleNoise_EigValSqrt:  " << sampleNoise_EigValSqrt << std::endl
      << " sampleNoise_Type:        " << sampleNoise_Type << std::endl
      //<< " sampleNoise_CovRandSeed: " << sampleNoise_CovRandSeed << std::endl
      << std::endl
      << "Custom Params: (IMLP)" << std::endl
      << " outlier_ChiSquareThreshold: " << outlier_ChiSquareThreshold << std::endl
      << " sigma2Max: " << sigma2Max << std::endl
      << " surfaceModel_InPlaneSD:   " << surfaceModel_InPlaneSD << std::endl
      << " surfaceModel_PerpPlaneSD: " << surfaceModel_PerpPlaneSD << std::endl
      << " bSampleCov_ApplyNoiseModel:     " << bSampleCov_ApplyNoiseModel << std::endl
      << " bSampleCov_ApplySurfaceModel:   " << bSampleCov_ApplySurfaceModel << std::endl
      << " bTargetCov_ApplySurfaceModel:   " << bTargetCov_ApplySurfaceModel << std::endl
      << " bTargetCov_DynamicSurfaceModel: " << bTargetCov_DynamicSurfaceModel << std::endl
      << std::endl
      //<< "Custom Params: (Robust ICP)" << std::endl
      //<< " D:     " << D << std::endl
      //<< " D0max: " << D0max << std::endl
      ;
    return ss.str();
  }

};



#endif