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
#include "ParameterizedTest_PointCloud_SurfaceNoise.h"


// Radomized trials for registering 2 point clouds
//  sample noise is dependent on surface normal orientation
//  (assumes mesh itself has no measurement noise)
void Run_ParameterizedTest_PointCloud_SurfaceNoise_Outliers(TestParameters params)
{
  std::stringstream ss;

  //std::string workingDir = "../ICP_TestData";
  //std::string baseOutputFolder = "RandomTestData_PointCloud_SurfaceNoise";  // "RandomTests_Plane2Plane_Term0.001G_5.0,0.5_ChiSqu0.95_A0_Sigma2/";
  //std::string sourceMeshFile = "RIGHTHEMIPELVIS_centered.mesh";
  //std::string targetMeshFile = sourceMeshFile;

  //// Create folders for output files
  //std::string baseOutputDir = workingDir + "/" + baseOutputFolder;
  //std::string commonOutputDir = baseOutputDir + "/CommonFiles";
  //CreateDir(baseOutputDir);
  //CreateDir(commonOutputDir);

  std::string loadSourceMeshPath = params.inputDir + "/" + params.sourceMeshFile;
  std::string loadTargetMeshPath = params.inputDir + "/" + params.targetMeshFile;
  std::string saveSourceMeshPath = params.outputCommonDir + "/SaveMeshSource";
  std::string saveTargetMeshPath = params.outputCommonDir + "/SaveMeshTarget";
  std::string saveNoisyTargetMeshPath = params.outputCommonDir + "/SaveMeshTargetNoisy";
  std::string saveTargetPointCloud = params.outputCommonDir + "/SaveTargetPointCloud.pts";
  std::string saveTargetCov = params.outputCommonDir + "/SaveTargetCov.txt";

  //// test parameters
  //TestParameters params;
  //params.TestDescription = "Run_ParameterizedTest_PointCloud_SurfaceNoise";
  //params.inputDir = workingDir;
  //params.sourceMeshFile = sourceMeshFile;
  //params.targetMeshFile = targetMeshFile;
  //params.numTrials = 2; // 300;
  //params.nSamples = 100;
  //params.nValidationSamples = 100;
  //params.noiseMeshSD = 0.0;
  //params.randSeed_MeshNoise = 0;
  params.GenerateSamplesFunc = GenerateSamplePointSet_PointCloud_SurfaceNoise;
  CustomTestParams_PointCloud_SurfaceNoise &customParams
    = dynamic_cast<CustomTestParams_PointCloud_SurfaceNoise&> (*params.pCustomTestParams);

  //// custom test parameters
  //// NOTE: target point cloud noise model
  ////        perpendicular noise must be specified
  ////        in-plane noise is set to avg distance from triangle center to each vertex
  //CustomTestParams_PointCloud_SurfaceNoise customParams;
  //params.pCustomTestParams = &customParams;
  //customParams.percentOutliers = 0.0;
  //customParams.bSampleCov_ApplyNoiseModel = true;
  //customParams.bSampleCov_ApplySurfaceModel = true;
  //customParams.bTargetCov_ApplySurfaceModel = true;
  //customParams.bTargetCov_DynamicSurfaceModel = false;
  //customParams.surfaceModel_InPlaneSD = 5.0;
  //customParams.surfaceModel_PerpPlaneSD = 0.5;
  //customParams.minPosOffsetOutlier = 10.0;
  //customParams.maxPosOffsetOutlier = 20.0;

  // load target mesh
  cisstMesh mesh_target;
  CreateMesh(mesh_target, loadTargetMeshPath, &saveTargetMeshPath);

  // load source mesh
  //  NOTE: the sample point set is randomly generated from this mesh
  cisstMesh mesh_source;
  CreateMesh(mesh_source, loadSourceMeshPath, &saveSourceMeshPath);

  // generate noisy target mesh
  //  NOTE: the sample point set is registered to this mesh
  cisstMesh noisyTargetMesh;
  if (params.noiseMeshSD != 0.0)
  {
    std::cout << "ERROR: noisy mesh currently disabled" << std::endl;
    return;
    //std::string normRVFile = params.inputDir + "/GaussianValues.txt";
    //std::ifstream randnStream_MeshNoise(normRVFile.c_str());
    //GenerateNoisyMesh(mesh_target, noisyTargetMesh, randnStream_MeshNoise, params.noiseMeshSD, &saveNoisyTargetMeshPath);
  }
  else
  {
    noisyTargetMesh = mesh_target;
  }

  // build point cloud covariance tree from noisy mesh
  //  and define point cloud noise model (but don't add more noise)
  cisstCovTree_PointCloud   *pTree;
  int    nThresh = 5;       // Cov Tree Params
  double diagThresh = 5.0;  //  ''
  std::cout << "Building covariance tree .... " << std::endl;
  double targetModel_InPlaneSD = params.noiseMeshSD;
  double targetModel_PerpPlaneSD = params.noiseMeshSD;
  if (customParams.bTargetCov_ApplySurfaceModel)
  {
    targetModel_InPlaneSD += customParams.surfaceModel_InPlaneSD;
    targetModel_PerpPlaneSD += customParams.surfaceModel_PerpPlaneSD;
  }
  if (customParams.bTargetCov_DynamicSurfaceModel && customParams.bTargetCov_ApplySurfaceModel)
  {
    // set target surface noise model dynamically from mesh
    pTree = new cisstCovTree_PointCloud(noisyTargetMesh, nThresh, diagThresh, targetModel_PerpPlaneSD);
    std::cout << "Tree built: NNodes=" << pTree->NumNodes() << " NData=" << pTree->NumData()
      << " TreeDepth=" << pTree->TreeDepth() << std::endl;
    std::cout << " Point Cloud Noise Model:" << std::endl
      << "  perp-plane variance = " << pow(customParams.surfaceModel_PerpPlaneSD, 2) << std::endl
      << "  in-plane variance = " << pTree->avgVarInPlane << " (avg)" << std::endl;
    assert(pTree->NumNodes());
  }
  else
  {
    // set static target surface noise model
    pTree = new cisstCovTree_PointCloud(noisyTargetMesh, nThresh, diagThresh, targetModel_PerpPlaneSD, targetModel_InPlaneSD);
    std::cout << "Tree built: NNodes=" << pTree->NumNodes() << " NData=" << pTree->NumData()
      << " TreeDepth=" << pTree->TreeDepth() << std::endl;
    std::cout << " Point Cloud Noise Model:" << std::endl
      << "  perp-plane variance = " << pow(customParams.surfaceModel_PerpPlaneSD, 2) << std::endl
      << "  in-plane variance = " << pow(customParams.surfaceModel_InPlaneSD, 2) << std::endl;
    assert(pTree->NumNodes());
  }
  // save target point cloud and covariance model
  pTree->SavePointCloud(saveTargetPointCloud);
  pTree->SavePointCloudCov(saveTargetCov);

  // ICP Options
  cisstICP::Options opt_ICP;
  opt_ICP.auxOutputDir = params.outputBaseDir + "/AuxiliaryFiles/";
  opt_ICP.maxIter = 100;
  opt_ICP.termHoldIter = 2;
  opt_ICP.minE = -std::numeric_limits<double>::max();
  opt_ICP.tolE = 0.0;
  opt_ICP.dPosThresh = 0.1;
  opt_ICP.dAngThresh = 0.1*(cmnPI / 180);
  opt_ICP.dPosTerm = 0.001;
  opt_ICP.dAngTerm = 0.001*(cmnPI / 180);
  params.opt_ICP = opt_ICP;

  //// CPD Options
  //CPD::Options opt_CPD;
  //opt_CPD.maxIter = opt_ICP.maxIter;
  //opt_CPD.termHoldIter = opt_ICP.termHoldIter;
  //opt_CPD.minE = opt_ICP.minE;
  //opt_CPD.tolE = opt_ICP.tolE;
  //opt_CPD.dPosThresh = opt_ICP.dPosThresh;
  //opt_CPD.dAngThresh = opt_ICP.dAngThresh;
  //opt_CPD.dPosTerm = opt_ICP.dPosTerm;
  //opt_CPD.dAngTerm = opt_ICP.dAngTerm;
  //opt_CPD.w = customParams.w;
  //params.opt_CPD = opt_CPD;


#define nNoises     9
#define nOffsets    2
#define nOutliers   6
#define nTests      (nNoises*nOffsets*nOutliers)

  double percentOutliers[nOutliers] = { 0.05, 0.1, 0.2, 0.3, 0.4, 0.5 };
  //double percentOutliers[nOutliers] = { 0.05, 0.1 };

  double ChiSqrThresh[nOutliers] = { 7.81, 6.25, 4.64, 3.66, 3.66, 3.66 };

  double noiseInPlaneSDArr[nNoises]   = { 0.5, 1.0, 2.0,  0.5, 1.0, 0.5,  1.0, 2.0, 2.0 };
  double noisePerpPlaneSDArr[nNoises] = { 0.5, 1.0, 2.0,  1.0, 2.0, 2.0,  0.5, 1.0, 0.5 };
  //double noiseInPlaneSDArr[nNoises] = { 0.5, 1.0 };
  //double noisePerpPlaneSDArr[nNoises] = { 0.5, 1.0 };

  double minOffsetAngArr[nOffsets] = { 30,15 }; //{ 10, 30 }; //{0,  10, 30, 60};
  double maxOffsetAngArr[nOffsets] = { 60,30 }; //{ 30, 60 }; //{10, 30, 60, 90};
  double minOffsetPosArr[nOffsets] = { 30,15 }; //{ 10, 30 }; //{0,  10, 30, 60};
  double maxOffsetPosArr[nOffsets] = { 60,30 }; //{ 30, 60 }; //{10, 30, 60, 100};

  std::vector<unsigned int> randSeed_SamplesArr(nTests,0);   // sample generater
  std::vector<unsigned int> randSeed_XfmsArr(nTests, 0);     // initial offset xfm generater

  for (unsigned int i = 0; i < nTests; i++)
  {
    randSeed_SamplesArr[i] = i + 23;
    randSeed_XfmsArr[i] = i + 127;
  }

  std::string outputCommonDir = params.outputCommonDir;


  //--- Parameterized Trials ---//

  for (unsigned int offset = 0; offset < nOffsets; offset++)
  {
    for (unsigned int noise = 0; noise < nNoises; noise++)
    {
      for (unsigned int outlier = 0; outlier < nOutliers; outlier++)
      {
        //unsigned int testNum = noise*nOffsets*nOutliers + offset*nOutliers + outlier;
        unsigned int testNum = offset*nNoises*nOutliers + noise*nOutliers + outlier;

        ss.str("");
        ss << "Noise_Prll" << noiseInPlaneSDArr[noise] << "Perp" << noisePerpPlaneSDArr[noise]
          << "_Offset_R" << minOffsetAngArr[offset] << "-" << maxOffsetAngArr[offset]
          << "T" << minOffsetPosArr[offset] << "-" << maxOffsetPosArr[offset]
          << "_Outliers" << percentOutliers[outlier];
        std::string outputSubFolder = ss.str();

        params.minOffsetAng = minOffsetAngArr[offset];
        params.maxOffsetAng = maxOffsetAngArr[offset];
        params.minOffsetPos = minOffsetPosArr[offset];
        params.maxOffsetPos = maxOffsetPosArr[offset];
        params.randSeed_Samples = randSeed_SamplesArr[testNum];
        params.randSeed_Xfms = randSeed_XfmsArr[testNum];

        customParams.sampleNoise_InPlaneSD = noiseInPlaneSDArr[noise];
        customParams.sampleNoise_PerpPlaneSD = noisePerpPlaneSDArr[noise];

        // hold number of good samples constant and increase sample size to
        //  achieve the desired percentage of outliers
        double p = percentOutliers[outlier];
        unsigned int nOutlierSamps = (unsigned int)(p/(1-p) * customParams.nGoodSamples);
        params.nSamples = customParams.nGoodSamples + nOutlierSamps;
        customParams.percentOutliers = p;

        customParams.outlier_ChiSquareThreshold = ChiSqrThresh[outlier];

        // TODO: remove change to output common directories below

        // Run StdICP
        {
          params.AlgorithmDescription = "StdICP";
          params.algType = TestParameters::StdICP;
          SetOutputDir(params, outputCommonDir, outputSubFolder);
          
          vctDynamicVector<vct3> dummySamples;      // these will be set later by the sample generator routine
          params.pAlg = new cisstAlgorithmICP_StdICP_PointCloud(pTree, dummySamples);
          Run_ParameterizedTest(mesh_source, pTree, params);
          delete params.pAlg; params.pAlg = NULL;
        }

        // Run RobustICP
        {
          params.AlgorithmDescription = "RobustICP";
          params.algType = TestParameters::RobustICP;
          SetOutputDir(params, outputCommonDir, outputSubFolder);

          vctDynamicVector<vct3> dummySamples;      // these will be set later by the sample generator routine
          params.pAlg = new cisstAlgorithmICP_RobustICP_PointCloud(pTree, dummySamples, customParams.D, customParams.D0max);
          Run_ParameterizedTest(mesh_source, pTree, params);
          delete params.pAlg; params.pAlg = NULL;
        }

        // Run IMLP
        {
          std::stringstream ss;
          ss << "IMLP_ChiSqu" << customParams.outlier_ChiSquareThreshold;
          params.AlgorithmDescription = ss.str();
          //params.AlgorithmDescription = "IMLP_OutlierB_SDExpFctr3";
          params.algType = TestParameters::IMLP;
          SetOutputDir(params, outputCommonDir, outputSubFolder);

          vctDynamicVector<vct3> dummySamples;      // these will be set later by the sample generator routine
          vctDynamicVector<vct3x3> dummySampleCov;  //  ''
          cisstAlgorithmICP_IMLP_PointCloud *pAlg =
            new cisstAlgorithmICP_IMLP_PointCloud(
            pTree, dummySamples, dummySampleCov, dummySampleCov, customParams.outlier_ChiSquareThreshold);
          params.pAlg = pAlg;
          Run_ParameterizedTest(mesh_source, pTree, params);
          delete params.pAlg; params.pAlg = NULL;
        }

      }
    }
  }
}
