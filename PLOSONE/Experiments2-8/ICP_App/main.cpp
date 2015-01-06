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
#include <stdio.h>
#include <iostream>
#include <vector>

//#include <cisstCommon.h>
#include <cisstVector.h>


#include "ParameterizedTest_Mesh_SurfaceNoise.h"
#include "ParameterizedTest_PointCloud_SurfaceNoise.h"
#include "ParameterizedTest_PointCloud_SubShape.h"
#include "ParameterizedTest_PointCloud_PartialOverlap.h"

#include "utility.h"

int main(void)
{	


#if 1
{
  // Common Settings
  std::string workingDir = "../ICP_TestData";
  TestParameters params;
  params.inputDir = workingDir;
  params.nValidationSamples = 100;
  params.noiseMeshSD = 0.0;
  params.randSeed_MeshNoise = 0;
  params.bEnableIterationCallbacks = true;

  // Experiment 2: Mesh Target Without Outliers
  {
    std::string baseOutputFolder = "RandomTestData_Mesh_SurfaceNoise";
    std::string baseOutputDir = workingDir + "/" + baseOutputFolder;
    std::string commonOutputDir = baseOutputDir + "/CommonFiles";
    CreateDir(baseOutputDir);
    CreateDir(commonOutputDir);

    // standard test parameters
    params.TestDescription = baseOutputFolder;
    params.outputBaseDir = baseOutputDir;
    params.outputCommonDir = commonOutputDir;
    params.targetMeshFile = "RIGHTHEMIPELVIS_centered.mesh";
    params.sourceMeshFile = params.targetMeshFile;
    params.numTrials = 300;
    //params.nSamples = 100;
    // custom test parameters (general)
    CustomTestParams_PointCloud_SurfaceNoise customParams;
    params.pCustomTestParams = &customParams;
    customParams.nGoodSamples = 100;
    customParams.minPosOffsetOutlier = 10.0;
    customParams.maxPosOffsetOutlier = 20.0;
    customParams.percentOutliers = 0.0;
    // custom test parameters (IMLP)
    customParams.outlier_ChiSquareThreshold = 10000.0;
    customParams.surfaceModel_InPlaneSD = 5.0;
    customParams.surfaceModel_PerpPlaneSD = 0.5;
    customParams.bSampleCov_ApplyNoiseModel = true;
    customParams.bSampleCov_ApplySurfaceModel = false;
    customParams.bTargetCov_ApplySurfaceModel = false;
    customParams.bTargetCov_DynamicSurfaceModel = false;
    // custom test parameters (RobustICP)
    std::cout << std::endl << "Computing Avg Neighbor Distance for RobustICP... (may take a while!)" << std::endl;
    customParams.D = ComputeAvgNeighborDistance(params.inputDir + "/" + params.targetMeshFile);
    customParams.D = 0.0;
    customParams.D0max = 0.0;

    Run_ParameterizedTest_Mesh_SurfaceNoise(params);
  }

  // Experiment 3: Mesh Target With Outliers
  {
    std::string baseOutputFolder = "RandomTestData_Mesh_SurfaceNoise_Outliers";
    std::string baseOutputDir = workingDir + "/" + baseOutputFolder;
    std::string commonOutputDir = baseOutputDir + "/CommonFiles";
    CreateDir(baseOutputDir);
    CreateDir(commonOutputDir);

    // standard test parameters
    params.TestDescription = baseOutputFolder;
    params.outputBaseDir = baseOutputDir;
    params.outputCommonDir = commonOutputDir;
    params.targetMeshFile = "RIGHTHEMIPELVIS_centered.mesh";
    params.sourceMeshFile = params.targetMeshFile;
    params.numTrials = 300;
    //params.nSamples = 105;
    // custom test parameters (General)
    CustomTestParams_PointCloud_SurfaceNoise customParams;
    params.pCustomTestParams = &customParams;
    customParams.nGoodSamples = 100;
    customParams.minPosOffsetOutlier = 10.0;
    customParams.maxPosOffsetOutlier = 20.0;
    // these are set in the test routine
    //customParams.percentOutliers = 0.05;
    //customParams.outlier_ChiSquareThreshold = 7.81;
    customParams.surfaceModel_InPlaneSD = 5.0;
    customParams.surfaceModel_PerpPlaneSD = 0.5;
    customParams.bSampleCov_ApplyNoiseModel = true;
    customParams.bSampleCov_ApplySurfaceModel = false;
    customParams.bTargetCov_ApplySurfaceModel = false;
    customParams.bTargetCov_DynamicSurfaceModel = false;
    // custom test parameters (RobustICP)
    std::cout << std::endl << "Computing Avg Neighbor Distance for RobustICP... (may take a while!)" << std::endl;
    customParams.D = ComputeAvgNeighborDistance(params.inputDir + "/" + params.targetMeshFile);
    customParams.D0max = 1000.0;
    std::cout << "...Avg Neighbor Distance: " << customParams.D << std::endl << std::endl;

    Run_ParameterizedTest_Mesh_SurfaceNoise_Outliers(params);
  }
}
#endif


#if 1
{
  //-- Random Point-Cloud Tests --//

  // Common Settings
  std::string workingDir = "../ICP_TestData";
  TestParameters params;
  params.inputDir = workingDir;
  params.nValidationSamples = 100;
  params.noiseMeshSD = 0.0;
  params.randSeed_MeshNoise = 0;
  params.bEnableIterationCallbacks = true;


  // Experiment 4: Point Cloud Target Without Outliers
  {
    std::string baseOutputFolder = "RandomTestData_PointCloud_SurfaceNoise";
    std::string baseOutputDir = workingDir + "/" + baseOutputFolder;
    std::string commonOutputDir = baseOutputDir + "/CommonFiles";
    CreateDir(baseOutputDir);
    CreateDir(commonOutputDir);

    // standard test parameters
    params.TestDescription = baseOutputFolder;
    params.outputBaseDir = baseOutputDir;
    params.outputCommonDir = commonOutputDir;
    params.targetMeshFile = "RIGHTHEMIPELVIS_centered.mesh";
    params.sourceMeshFile = params.targetMeshFile;
    params.numTrials = 300;
    //params.nSamples = 100;
    // custom test parameters (general)
    CustomTestParams_PointCloud_SurfaceNoise customParams;
    params.pCustomTestParams = &customParams;
    customParams.nGoodSamples = 100;
    customParams.minPosOffsetOutlier = 10.0;
    customParams.maxPosOffsetOutlier = 20.0;
    customParams.percentOutliers = 0.0;
    // custom test parameters (IMLP)
    customParams.outlier_ChiSquareThreshold = 10000.0;
    customParams.surfaceModel_InPlaneSD = 5.0;
    customParams.surfaceModel_PerpPlaneSD = 0.5;
    customParams.bSampleCov_ApplyNoiseModel = true;
    customParams.bSampleCov_ApplySurfaceModel = true;
    customParams.bTargetCov_ApplySurfaceModel = true;
    customParams.bTargetCov_DynamicSurfaceModel = false;
    // custom test parameters (RobustICP)
    std::cout << std::endl << "Computing Avg Neighbor Distance for RobustICP... (may take a while!)" << std::endl;
    customParams.D = ComputeAvgNeighborDistance(params.inputDir + "/" + params.targetMeshFile);
    customParams.D0max = 0.0;

    Run_ParameterizedTest_PointCloud_SurfaceNoise(params);
  }


  // Experiment 5: Point Cloud Target With Outliers
  {
    std::string baseOutputFolder = "RandomTestData_PointCloud_SurfaceNoise_Outliers";
    std::string baseOutputDir = workingDir + "/" + baseOutputFolder;
    std::string commonOutputDir = baseOutputDir + "/CommonFiles";
    CreateDir(baseOutputDir);
    CreateDir(commonOutputDir);

    // standard test parameters
    params.TestDescription = baseOutputFolder;
    params.outputBaseDir = baseOutputDir;
    params.outputCommonDir = commonOutputDir;
    params.targetMeshFile = "RIGHTHEMIPELVIS_centered.mesh";
    params.sourceMeshFile = params.targetMeshFile;
    params.numTrials = 300;
    //params.nSamples = 105;
    // custom test parameters (General)
    CustomTestParams_PointCloud_SurfaceNoise customParams;
    params.pCustomTestParams = &customParams;
    customParams.nGoodSamples = 100;
    customParams.minPosOffsetOutlier = 10.0;
    customParams.maxPosOffsetOutlier = 20.0;
    //customParams.percentOutliers = 0.05;
    // custom test parameters (IMLP)
    // Note:    ChiSquare(0.6) = 2.95      
    //          ChiSquare(0.7) = 3.66      
    //          ChiSquare(0.8) = 4.64      
    //          ChiSquare(0.85) = 5.32     
    //          ChiSquare(0.9) = 6.25      
    //          ChiSquare(0.925) = 6.90    
    // default  ChiSquare(0.95) = 7.81     (1.96 Std Dev)
    //          ChiSquare(0.975) = 9.35    (2.24 Std Dev)
    //          ChiSquare(0.99) = 11.34    (2.56 Std Dev)
    //          ChiSquare(0.9973) = 14.16  (3.0 Std Dev)     MATLAB: chi2inv(0.9973,3)
    //customParams.outlier_ChiSquareThreshold = 7.81;
    //customParams.outlier_ChiSquareThreshold = 6.25;
    //customParams.outlier_ChiSquareThreshold = 4.64;
    //customParams.outlier_ChiSquareThreshold = 3.66;
    //customParams.outlier_ChiSquareThreshold = 6.90;    
    //customParams.outlier_ChiSquareThreshold = 5.32;
    customParams.surfaceModel_InPlaneSD = 5.0;
    customParams.surfaceModel_PerpPlaneSD = 0.5;
    customParams.bSampleCov_ApplyNoiseModel = true;
    customParams.bSampleCov_ApplySurfaceModel = true;
    customParams.bTargetCov_ApplySurfaceModel = true;
    customParams.bTargetCov_DynamicSurfaceModel = false;
    // custom test parameters (RobustICP)
    std::cout << std::endl << "Computing Avg Neighbor Distance for RobustICP... (may take a while!)" << std::endl;
    customParams.D = ComputeAvgNeighborDistance(params.inputDir + "/" + params.targetMeshFile);
    customParams.D0max = 1000.0;
    std::cout << "...Avg Neighbor Distance: " << customParams.D << std::endl << std::endl;

    Run_ParameterizedTest_PointCloud_SurfaceNoise_Outliers(params);
  }
}
#endif


#if 1
{
  //-- Random Sub-Shape Tests (Point Cloud) --//

  // Common Settings
  std::string workingDir = "../ICP_TestData";
  TestParameters params;
  params.inputDir = workingDir;
  params.nValidationSamples = 100;
  params.noiseMeshSD = 0.0;
  params.randSeed_MeshNoise = 0;
  params.bEnableIterationCallbacks = true;


  // Experiment 6: Sub-Shape Source to Point Cloud Target Without Outliers
  {
    std::string baseOutputFolder = "RandomTestData_PointCloud_SubShape";
    std::string baseOutputDir = workingDir + "/" + baseOutputFolder;
    std::string commonOutputDir = baseOutputDir + "/CommonFiles";
    CreateDir(baseOutputDir);
    CreateDir(commonOutputDir);

    // standard test parameters
    params.TestDescription = baseOutputFolder;
    params.outputBaseDir = baseOutputDir;
    params.outputCommonDir = commonOutputDir;
    params.targetMeshFile = "ProximalFemur.mesh";
    params.sourceMeshFile = "ProximalFemur_Crop.mesh";
    params.numTrials = 300;
    params.nSamples = 100;
    // custom test parameters (general)
    CustomTestParams_PointCloud_SubShape customParams;
    params.pCustomTestParams = &customParams;
    //customParams.nGoodSamples = 100;
    //customParams.minPosOffsetOutlier = 10.0;
    //customParams.maxPosOffsetOutlier = 20.0;
    //customParams.percentOutliers = 0.0;
    // custom test parameters (IMLP)
    customParams.outlier_ChiSquareThreshold = 10000.0;
    customParams.surfaceModel_InPlaneSD = 5.0;
    customParams.surfaceModel_PerpPlaneSD = 0.5;
    customParams.bSampleCov_ApplyNoiseModel = true;
    customParams.bSampleCov_ApplySurfaceModel = true;
    customParams.bTargetCov_ApplySurfaceModel = true;
    customParams.bTargetCov_DynamicSurfaceModel = false;
    // custom test parameters (RobustICP)
    //std::cout << std::endl << "Computing Avg Neighbor Distance for RobustICP... (may take a while!)" << std::endl;
    //customParams.D = ComputeAvgNeighborDistance(params.inputDir + "/" + params.targetMeshFile);
    //customParams.D = 0.0;
    //customParams.D0max = 0.0;

    Run_ParameterizedTest_PointCloud_SubShape(params);
  }

}
#endif

#if 1
{
  //-- Random Partial-Overlap Tests (Point Cloud) --//

  // Common Settings
  std::string workingDir = "../ICP_TestData";
  TestParameters params;
  params.inputDir = workingDir;
  params.nValidationSamples = 100;
  params.noiseMeshSD = 0.0;
  params.randSeed_MeshNoise = 0;
  params.bEnableIterationCallbacks = true;


  // Experiment 6: Sub-Shape Source to Point Cloud Target Without Outliers
  {
    std::string baseOutputFolder = "RandomTestData_PointCloud_PartialOverlap";
    std::string baseOutputDir = workingDir + "/" + baseOutputFolder;
    std::string commonOutputDir = baseOutputDir + "/CommonFiles";
    CreateDir(baseOutputDir);
    CreateDir(commonOutputDir);

    // standard test parameters
    params.TestDescription = baseOutputFolder;
    params.outputBaseDir = baseOutputDir;
    params.outputCommonDir = commonOutputDir;
    params.targetMeshFile = "sforza_bust50k_recentered_front.mesh";
    params.sourceMeshFile = "sforza_bust50k_recentered_right.mesh";
    params.numTrials = 10;
    //params.nSamples = 100;
    // custom test parameters (general)
    CustomTestParams_PointCloud_PartialOverlap customParams;
    params.pCustomTestParams = &customParams;
    //customParams.nGoodSamples = 100;
    //customParams.minPosOffsetOutlier = 10.0;
    //customParams.maxPosOffsetOutlier = 20.0;
    //customParams.percentOutliers = 0.0;
    // custom test parameters (IMLP)
    customParams.outlier_ChiSquareThreshold = 7.81;
    customParams.sigma2Max = 0.1;
    customParams.surfaceModel_InPlaneSD = 5.0;
    customParams.surfaceModel_PerpPlaneSD = 0.5;
    customParams.bSampleCov_ApplyNoiseModel = false;
    customParams.bSampleCov_ApplySurfaceModel = true;
    customParams.bTargetCov_ApplySurfaceModel = true;
    customParams.bTargetCov_DynamicSurfaceModel = false;
    // custom test parameters (RobustICP)
    //std::cout << std::endl << "Computing Avg Neighbor Distance for RobustICP... (may take a while!)" << std::endl;
    //customParams.D = ComputeAvgNeighborDistance(params.inputDir + "/" + params.targetMeshFile);
    //customParams.D = 0.0;
    //customParams.D0max = 0.0;

    Run_ParameterizedTest_PointCloud_PartialOverlap(params);
  }
}
#endif

  return 0;
}

