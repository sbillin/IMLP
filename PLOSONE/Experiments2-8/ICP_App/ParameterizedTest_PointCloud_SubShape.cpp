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
#include "ParameterizedTest_PointCloud_SubShape.h"
#include "utility.h"


void GenerateSamplePointSet_PointCloud_SubShape(
  cisstMesh &mesh,
  TestParameters &params,
  unsigned int randSeed, unsigned int &randSeqPos,
  std::ifstream &randnStream,
  vctDynamicVector<vct3> &samples,
  vctDynamicVector<vct3> &sampleNorms,
  vctDynamicVector<vct3> &noisySamples,
  unsigned int trialNum)
{
  // get custom params
  CustomTestParams_PointCloud_SubShape *pCustomParams =
    dynamic_cast<CustomTestParams_PointCloud_SubShape*>(params.pCustomTestParams);

  std::stringstream saveSamplesPath;
  std::stringstream saveNoisySamplesPath;
  std::stringstream saveNoiseCovPath;

  saveSamplesPath << params.outputCommonDir << "/SaveSamples_" << trialNum << ".pts";
  saveNoisySamplesPath << params.outputCommonDir << "/SaveNoisySamples_" << trialNum << ".pts";
  saveNoiseCovPath << params.outputCommonDir << "/SaveNoiseCov_" << trialNum << ".txt";

  std::string  strSaveSamplesPath = saveSamplesPath.str();
  std::string  strSaveNoisySamplesPath = saveNoisySamplesPath.str();
  std::string  strSaveNoiseCovPath = saveNoiseCovPath.str();

  std::string *pSaveSamplesPath = &strSaveSamplesPath;
  std::string *pSaveNoisySamplesPath = &strSaveNoisySamplesPath;
  std::string *pSaveNoiseCovPath = &strSaveNoiseCovPath;
#ifdef DISABLE_OUTPUT_FILES
  pSaveSamplesPath = NULL;
  pSaveNoisySamplesPath = NULL;
  pSaveNoiseCovPath = NULL;
#endif

  vctDynamicVector<unsigned int>  sampleDatums(params.nSamples);
  vctDynamicVector<vct3x3>        sampleCov(params.nSamples, vct3x3(0.0));
  vctDynamicVector<vct3x3>        noiseCov(params.nSamples, vct3x3(0.0));
  vctDynamicVector<vct3x3>        surfaceModelCov(params.nSamples, vct3x3(0.0));
  //vctDynamicVector<vct3x3>        noiseInvCov(params.nSamples, vct3x3(0.0));

  // Generate Samples
  GenerateSamples(
    mesh,
    randSeed, randSeqPos,
    params.nSamples,
    samples, sampleNorms, sampleDatums,
    pSaveSamplesPath);

  // Generate Sample Noise
  switch (pCustomParams->sampleNoise_Type)
  {
  case CustomTestParams_PointCloud_SubShape::NoiseSurface:
  {
    // generate covariances oriented to the surface
    ComputeCovariances_SurfaceModel(
      //pCustomParams->sampleNoise_InPlaneSD, pCustomParams->sampleNoise_PerpPlaneSD,
      pCustomParams->sampleNoise_EigValSqrt[1], pCustomParams->sampleNoise_EigValSqrt[0],
      samples, sampleNorms,
      noiseCov,
      pSaveNoiseCovPath);
    break;
  }
  case CustomTestParams_PointCloud_SubShape::NoiseGlobal:
  {
    // generate one global covariance for all samples
    vctDynamicVector<vct3x3> cov;
    ComputeCovariances_Random(
      randSeed, randSeqPos, 
      pCustomParams->sampleNoise_EigValSqrt,
      cov, 1,
      pSaveNoiseCovPath);
    noiseCov.SetAll(cov(0));
    break;
  }
  case CustomTestParams_PointCloud_SubShape::NoiseRandom:
  {
    // generate different random covariance for each sample
    ComputeCovariances_Random(
      randSeed, randSeqPos,
      pCustomParams->sampleNoise_EigValSqrt,
      noiseCov, params.nSamples,
      pSaveNoiseCovPath);
    break;
  }
  default:
  {
    std::cout << "ERROR: unrecognized noise type" << std::endl;
    assert(0);
  }
  }

  // generate sample noise according to the defined covariance
  GenerateSampleErrors_Covariance(
    randnStream,
    noiseCov,
    samples, sampleNorms,
    noisySamples,    
    pSaveNoisySamplesPath);

  // Set Samples in Algorithm
  switch (params.algType)
  {
  case TestParameters::StdICP:
  {
    // configure algorithm samples
    cisstAlgorithmICP_StdICP *pAlg = dynamic_cast<cisstAlgorithmICP_StdICP*>(params.pAlg);
    pAlg->SetSamples(noisySamples);
    break;
  }

  case TestParameters::IMLP:
  {
    std::stringstream saveSampleCovPath;
    std::stringstream saveSurfaceModelCovPath;
    saveSampleCovPath << params.outputDataDir << "/SaveSampleCov_" << trialNum << ".txt";
    saveSurfaceModelCovPath << params.outputDataDir << "/SaveSurfaceModelCov_" << trialNum << ".txt";

    std::string strSaveSampleCovPath = saveSampleCovPath.str();
    std::string strSaveSurfaceModelCovPath = saveSurfaceModelCovPath.str();

    std::string *pSaveSampleCovPath = &strSaveSampleCovPath;
    std::string *pSaveSurfaceModelCovPath = &strSaveSurfaceModelCovPath;
#ifdef MINIMIZE_OUTPUT_FILES
    pSaveSampleCovPath = NULL;
    pSaveSurfaceModelCovPath = NULL;
#endif
#ifdef DISABLE_OUTPUT_FILES
    pSaveSampleCovPath = NULL;
    pSaveSurfaceModelCovPath = NULL;
#endif

    unsigned int nSamps = samples.size();

    // Set Sample Covariances
    // noise covariance
    if (pCustomParams->bSampleCov_ApplyNoiseModel)
    {
      for (unsigned int i = 0; i < nSamps; i++)
      {
        sampleCov(i) += noiseCov(i);
      }
    }
    // surface model covariance
    if (pCustomParams->bSampleCov_ApplySurfaceModel)
    {
      // compute surface model
      ComputeCovariances_SurfaceModel(
        pCustomParams->surfaceModel_InPlaneSD, pCustomParams->surfaceModel_PerpPlaneSD,
        samples, sampleNorms,
        surfaceModelCov,
        pSaveSurfaceModelCovPath);
      // add suface model to sample covariances
      for (unsigned int i = 0; i < nSamps; i++)
      {
        sampleCov(i) += surfaceModelCov(i);
      }
    }
    // save sample covariances
    if (pSaveSampleCovPath)
    {
      WriteToFile_Cov(sampleCov, *pSaveSampleCovPath);
    }

    // configure algorithm samples
    cisstAlgorithmICP_IMLP *pAlg = dynamic_cast<cisstAlgorithmICP_IMLP*>(params.pAlg);
    //cisstAlgorithmICP_IMLP_PointCloud *pAlg = dynamic_cast<cisstAlgorithmICP_IMLP_PointCloud*>(params.pAlg);
    if (!pAlg)
    {
      std::cerr << "ERROR: algorithm class unrecognized class" << std::endl;
    }
    //pAlg->SetSamples(noisySamples);
    //pAlg->SetSampleCovariances(sampleCov);
    pAlg->SetSamples(noisySamples);
    pAlg->SetSampleCovariances(sampleCov, noiseCov);

    break;
  }

  case TestParameters::RobustICP:
  {
    // configure algorithm samples
    cisstAlgorithmICP_RobustICP *pAlg = dynamic_cast<cisstAlgorithmICP_RobustICP*>(params.pAlg);
    pAlg->SetSamples(noisySamples);
    break;
  }

  default:
    std::cout << std::endl << "=====> ERROR: Algorithm Type not Recognized by Sample Generator" << std::endl << std::endl;
    assert(0);
    break;
  } // switch AlgType
}


void SetOutputDir_SubShape(TestParameters &params, std::string outputCommonDir, std::string outputSubFolder)
{
  std::string outputSubDataDir = params.outputBaseDir + "/" + params.AlgorithmDescription;
  params.outputDataDir = outputSubDataDir + "/" + outputSubFolder;
  CreateDir(outputSubDataDir);

  std::string outputSubCommonDir = outputCommonDir + "_" + params.AlgorithmDescription;
  params.outputCommonDir = outputSubCommonDir + "/" + outputSubFolder;
  CreateDir(outputSubCommonDir);
  //params.outputCommonDir = outputCommonDir + "/" + outputSubFolder;
}



void Run_ParameterizedTest_PointCloud_SubShape(TestParameters params)
{
  std::stringstream ss;

  std::string loadSourceMeshPath = params.inputDir + "/" + params.sourceMeshFile;
  std::string loadTargetMeshPath = params.inputDir + "/" + params.targetMeshFile;
  std::string saveSourceMeshPath = params.outputCommonDir + "/SaveMeshSource";
  std::string saveTargetMeshPath = params.outputCommonDir + "/SaveMeshTarget";
  std::string saveNoisyTargetMeshPath = params.outputCommonDir + "/SaveMeshTargetNoisy";

  params.GenerateSamplesFunc = GenerateSamplePointSet_PointCloud_SubShape;
  CustomTestParams_PointCloud_SubShape &customParams
    = dynamic_cast<CustomTestParams_PointCloud_SubShape&> (*params.pCustomTestParams);

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
  if (params.noiseMeshSD > 0.0)
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
  std::string saveTargetPointCloud = params.outputCommonDir + "/SaveTargetPointCloud.pts";
  std::string saveTargetCov = params.outputCommonDir + "/SaveTargetCov.txt";
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

#define nNoises     7
#define nOffsets    1
#define nTests      (nNoises*nOffsets)
  //#define nOutliers   6
  //#define nTests      (nNoises*nOffsets*nOutliers)
  //double percentOutliers[nOutliers] = { 0.05, 0.1, 0.2, 0.3, 0.4, 0.5 };

  // used for surface model noise
  // used for random and global noise
  //double noiseInPlaneSDArr[nNoises]   = { 0.0, 0.5, 1.0, 0.5, 1.0, 0.0, 0.0 };
  //double noisePerpPlaneSDArr[nNoises] = { 0.0, 0.5, 1.0, 1.0, 0.5, 0.0, 0.0 };
  double noiseEigValSqrt1[nNoises] =  { 0.0, 0.5, 1.0, 0.5, 1.0, 0.5, 0.5 };  // perp-plane std dev for surf noise
  double noiseEigValSqrt23[nNoises] = { 0.0, 0.5, 1.0, 1.0, 0.5, 1.0, 1.0 };  // in-plane std dev for surf noise
  CustomTestParams_PointCloud_SubShape::NoiseType noiseType[nNoises] = {
    CustomTestParams_PointCloud_SubShape::NoiseSurface,
    CustomTestParams_PointCloud_SubShape::NoiseSurface,
    CustomTestParams_PointCloud_SubShape::NoiseSurface,
    CustomTestParams_PointCloud_SubShape::NoiseSurface,
    CustomTestParams_PointCloud_SubShape::NoiseSurface,
    CustomTestParams_PointCloud_SubShape::NoiseGlobal,
    CustomTestParams_PointCloud_SubShape::NoiseRandom
  };
  std::string noiseTypeStr[nNoises];
  for (unsigned int i = 0; i < nNoises; i++)
  {
    switch (noiseType[i])
    {
    case CustomTestParams_PointCloud_SubShape::NoiseSurface:
      noiseTypeStr[i] = "Surf";
      break;
    case CustomTestParams_PointCloud_SubShape::NoiseGlobal:
      noiseTypeStr[i] = "Glob";
      break;
    case CustomTestParams_PointCloud_SubShape::NoiseRandom:
      noiseTypeStr[i] = "Rand";
      break;
    default:
      std::cout << "ERROR: unrecognized noise type" << std::endl;
      noiseTypeStr[i] = "ERROR";
      assert(0);
    }
  }

  double minOffsetAngArr[nOffsets] = { 10 };
  double maxOffsetAngArr[nOffsets] = { 20 };
  double minOffsetPosArr[nOffsets] = { 10 };
  double maxOffsetPosArr[nOffsets] = { 20 };

  std::vector<unsigned int> randSeed_SamplesArr(nTests, 0);   // sample generater
  std::vector<unsigned int> randSeed_XfmsArr(nTests, 0);     // initial offset xfm generater

  for (unsigned int i = 0; i < nTests; i++)
  {
    randSeed_SamplesArr[i] = i + 23;
    randSeed_XfmsArr[i] = i + 127;
  }

  std::string outputCommonDir = params.outputCommonDir;
  //params.nSamples = customParams.nGoodSamples;

  //--- Parameterized Trials ---//

  for (unsigned int offset = 0; offset < nOffsets; offset++)
  {
    for (unsigned int noise = 0; noise < nNoises; noise++)
    {
      unsigned int testNum = offset*nNoises + noise;

      ss.str("");
      //ss << "Noise_Prll" << noiseInPlaneSDArr[noise] << "Perp" << noisePerpPlaneSDArr[noise]
      ss << "Noise_" << noiseTypeStr[noise] << "EigVal" << noiseEigValSqrt1[noise] << "-" << noiseEigValSqrt23[noise]
        << "_Offset_R" << minOffsetAngArr[offset] << "-" << maxOffsetAngArr[offset]
        << "T" << minOffsetPosArr[offset] << "-" << maxOffsetPosArr[offset];        
      //<< "_Outliers" << percentOutliers[outlier];
      std::string outputSubFolder = ss.str();

      params.minOffsetAng = minOffsetAngArr[offset];
      params.maxOffsetAng = maxOffsetAngArr[offset];
      params.minOffsetPos = minOffsetPosArr[offset];
      params.maxOffsetPos = maxOffsetPosArr[offset];
      params.randSeed_Samples = randSeed_SamplesArr[testNum];
      params.randSeed_Xfms = randSeed_XfmsArr[testNum];

      //customParams.sampleNoise_InPlaneSD = noiseInPlaneSDArr[noise];
      //customParams.sampleNoise_PerpPlaneSD = noisePerpPlaneSDArr[noise];
      customParams.sampleNoise_EigValSqrt = vct3(noiseEigValSqrt1[noise], noiseEigValSqrt23[noise], noiseEigValSqrt23[noise]);
      customParams.sampleNoise_Type = noiseType[noise];

      // hold number of good samples constant and increase sample size to
      //  achieve the desired percentage of outliers
      //double p = percentOutliers[outlier];
      //unsigned int nOutlierSamps = (unsigned int)(p/(1-p) * customParams.nGoodSamples);
      //params.nSamples = customParams.nGoodSamples + nOutlierSamps;      


      // TODO: remove change to output common files to separate directories below

      // Run StdICP
      {
        params.AlgorithmDescription = "StdICP";
        //params.AlgorithmDescription = "StdICP_Verbose";
        params.algType = TestParameters::StdICP;
        SetOutputDir_SubShape(params, outputCommonDir, outputSubFolder);

        vctDynamicVector<vct3> dummySamples;      // these will be set later by the sample generator routine
        params.pAlg = new cisstAlgorithmICP_StdICP_PointCloud(pTree, dummySamples);
        Run_ParameterizedTest(mesh_source, pTree, params);
        delete params.pAlg; params.pAlg = NULL;
      }

      // Run IMLP
      {
        params.AlgorithmDescription = "IMLP";
        //params.AlgorithmDescription = "IMLP_Verbose";
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

      // Run IMLP_MD (Mahalanobis Distance Match Criteria)
      {
        params.AlgorithmDescription = "IMLP-MD";
        //params.AlgorithmDescription = "IMLP-MD_Verbose";
        params.algType = TestParameters::IMLP;
        SetOutputDir(params, outputCommonDir, outputSubFolder);

        vctDynamicVector<vct3> dummySamples;      // these will be set later by the sample generator routine
        vctDynamicVector<vct3x3> dummySampleCov;  //  ''
        cisstAlgorithmICP_IMLP_MahalDist_PointCloud *pAlg =
          new cisstAlgorithmICP_IMLP_MahalDist_PointCloud(
          pTree, dummySamples, dummySampleCov, dummySampleCov, customParams.outlier_ChiSquareThreshold);
        params.pAlg = pAlg;
        Run_ParameterizedTest(mesh_source, pTree, params);
        delete params.pAlg; params.pAlg = NULL;
      }

      // Run IMLP_CP (Closest Point Match Criteria)
      {
        params.AlgorithmDescription = "IMLP-CP";
        //params.AlgorithmDescription = "IMLP-CP_Verbose";
        params.algType = TestParameters::IMLP;
        SetOutputDir(params, outputCommonDir, outputSubFolder);

        vctDynamicVector<vct3> dummySamples;      // these will be set later by the sample generator routine
        vctDynamicVector<vct3x3> dummySampleCov;  //  ''
        cisstAlgorithmICP_IMLP_ClosestPoint_PointCloud *pAlg =
          new cisstAlgorithmICP_IMLP_ClosestPoint_PointCloud(
          pTree, dummySamples, dummySampleCov, dummySampleCov, customParams.outlier_ChiSquareThreshold);
        params.pAlg = pAlg;
        Run_ParameterizedTest(mesh_source, pTree, params);
        delete params.pAlg; params.pAlg = NULL;
      }

    }
  }
}

