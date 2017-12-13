#ifndef _testICPNormals_H
#define _testICPNormals_H

#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <limits.h>

#include "utility.h"
#include "cisstICP.h"
#include "cisstMesh.h"
#include "DirPDTree_Mesh.h"
#include "DirPDTree_PointCloud.h"

#include "algDirICP_StdICP_PointCloud.h"
#include "algDirICP_StdICP_Mesh.h"
#include "algDirICP_IMLOP_Mesh.h"
//#include "algDirICP_GIMLOP_Mesh.h"
#include "algDirICP_PIMLOP_Mesh.h"

// disable for run-time tests
#define ENABLE_CALLBACKS

enum ICPDirAlgType { DirAlgType_StdICP, DirAlgType_IMLOP, DirAlgType_PIMLOP };  // DirAlgType_GIMLOP, 


void Callback_TrackRegPath_testICPNormals(cisstICP::CallbackArg &arg, void *userData)
{
  // cast to norm callback arg type
  cisstICP::CallbackArg *argp = &arg;
  //argp = dynamic_cast<cisstICP::CallbackArgNormals*>(&arg);
  //if (argp==0)
  //{
  //  std::cerr << "ERROR: cannot cast callback argument to cisstICP type" << std::endl;
  //  assert(argp);   // terminate program
  //}

  // Save to file:
  //  - error function (-loglik)
  //  - incremental transform
  //  - IMLOP params
  //  - residual match errors
  // output format:
  //  err r00 r01 r02 r10 r11 r12 r20 r21 r22 tx ty tz normWeight posWeight avgNormError avgPosError 
  std::ofstream *fs = (std::ofstream *)(userData);
  (*fs) << argp->E << " " << argp->dF.Rotation().Row(0) << " " << argp->dF.Rotation().Row(1) << " "
    << argp->dF.Rotation().Row(2) << " " << argp->dF.Translation() << " "
    //<< argp->normWeight << " " << argp->posWeight << " "
    //<< argp->MatchPosErrorAvg << " " << argp->MatchPosErrorSD << " "
    //<< argp->MatchNormErrorAvg << " " << argp->MatchNormErrorSD 
    << std::endl;

}
void Callback_SaveIterationsToFile_testICPNormals(cisstICP::CallbackArg &arg, void *userData)
{
  std::ofstream *fs = (std::ofstream *)(userData);

  // cast to norm callback arg type
  cisstICP::CallbackArg *argp;
  //argp = dynamic_cast<cisstICP::CallbackArgNormals*>(&arg);
  //if (argp==0)
  //{
  //  std::cerr << "ERROR: cannot cast callback argument to cisstICP type" << std::endl;
  //  assert(argp);   // terminate program
  //}

  vctRodRot3 dR(arg.dF.Rotation());
  std::stringstream ss;
  ss << cmnPrintf("iter=%u  E=%.3f  tolE=%.4f  nW/pW=%.4f/%.4f  (cSD/dSD)=%.2f/%.2f  dTheta=%.2f/%.2f/%.2f  (dAng/dPos)= %.2f/%.2f  t=%.4f NNodes=%u/%u/%u NOut=%u/%u")
    // (RMS/Res)=%.4f/%.4f
    //  dTheta=%.2f/%.2f/%.2f  
    << argp->iter
    << argp->E
    << argp->tolE
    //<< argp->normWeight
    //<< argp->posWeight
    //<< argp->circSD*180.0 / cmnPI << sqrt(argp->posVar)
    //<< argp->dThetaMin*180.0 / cmnPI << argp->dThetaMax*180.0 / cmnPI << argp->dThetaAvg*180.0 / cmnPI
    << dR.Norm() * 180 / cmnPI << arg.dF.Translation().Norm()
    << argp->time;
  //<< argp->maxNodesSearched << argp->avgNodesSearched << argp->minNodesSearched
  //<< argp->nOutliersPos
  //<< argp->nOutliersNorm;    

  (*fs) << ss.str() << std::endl;
}

// TargetShapeAsMesh    true - uses mesh to represent target shape
//                      false - uses point cloud (taken from mesh) to represent target shape
void testICPNormals(bool TargetShapeAsMesh, ICPDirAlgType algType)
{
  //char k;
  //std::cout << "press key then enter to start:" << std::endl;
  //std::cin >> k;

  int    nThresh = 15;       // Cov Tree Params
  double diagThresh = 15.0;  //  ''
  //int    nThresh = 5;       // Cov Tree Params
  //double diagThresh = 5.0;  //  ''

  std::string workingDir = "..//test_data//";
  std::string outputDir = "LastRun//";

  std::string saveSourceMeshPath = workingDir + outputDir + "SaveMeshSource";
  std::string saveTargetMeshPath = workingDir + outputDir + "SaveMeshTarget";
  std::string saveSamplesPath = workingDir + outputDir + "SaveSamples";
  std::string saveNoisySamplesPath = workingDir + outputDir + "SaveNoisySamples";
  std::string savePath_L = workingDir + outputDir + "SaveL.pts";
  std::string savePath_Cov = workingDir + outputDir + "SaveCov.pts";

  cisstMesh         mesh_source;
  cisstMesh         mesh_target;
  DirPDTreeBase*      pTree;
  vctDynamicVector<vct3>    samples;
  vctDynamicVector<vct3>    sampleNorms;
  vctDynamicVector<vct3>    noisySamples;
  vctDynamicVector<vct3>    noisySampleNorms;
  vctDynamicVector<unsigned int>  sampleDatums;
  vctDynamicVector<vct3x3>  sampleNoiseCov;
  vctDynamicVector<vct3x3>  sampleNoiseInvCov;
  vctDynamicVector<vct3x2>  sampleNoiseL;

#if 1
  std::string loadSourceMeshPath = workingDir + "ProximalFemur.ply";
  std::string loadTargetMeshPath = workingDir + "ProximalFemur.ply";
  //std::string loadSourceMeshPath = workingDir + "RIGHTHEMIPELVIS_centered.mesh";
  //std::string loadTargetMeshPath = workingDir + "RIGHTHEMIPELVIS_centered.mesh";
  //std::string loadMeshPath = workingDir + "RIGHTHEMIPELVIS_centered.mesh";
  //std::string loadMeshPath = workingDir + "RIGHTHEMIPELVIS.mesh";
  //std::string loadMeshPath = workingDir + "CTBreastImage_Dec20000_Shell.mesh";  

  const int nSamples = 100;

  // Samples Noise Model
  //  NOTE: this is a generative noise model (i.e. noise is generated according
  //        to the noise properties defined here)
  //double noiseSampsSD[3] = {1.0, 1.0, 1.0};   // noise model for samples (std dev along each axis)
  double sampleNoiseInPlane = 1.0;      // standard deviation of noise in and out of plane
  double sampleNoisePerpPlane = 1.0;    //   ''
  double sampleNoiseCircSDDeg = 2.0;   // noise to apply to sample orientations
  double sampleNoiseCircSD = sampleNoiseCircSDDeg*cmnPI / 180.0;
  double sampleNoiseEccentricity = 0.5; // eccentricity of orientation noise

  double minOffsetPos = 5.0;
  double maxOffsetPos = 10.0;
  double minOffsetAng = 10.0;
  double maxOffsetAng = 20.0;

  unsigned int randSeed1 = 0;     // generates samples
  unsigned int randSeqPos1 = 0;
  unsigned int randSeed2 = 1;     // generates offsets
  unsigned int randSeqPos2 = 0;

  double percentOutliers = 0.0;
  double minPosOffsetOutlier = 15.0;
  double maxPosOffsetOutlier = 20.0;
  double minAngOffsetOutlier = 15.0;
  double maxAngOffsetOutlier = 20.0;


  // load source mesh
  CreateMesh(mesh_source, loadSourceMeshPath, &saveSourceMeshPath);
  // load target mesh
  CreateMesh(mesh_target, loadTargetMeshPath, &saveTargetMeshPath);

  if (TargetShapeAsMesh)
  {
    // build PD tree on the mesh directly
    //  Note: defines measurement noise to be zero
    printf("Building mesh PD tree .... ");
    pTree = new DirPDTree_Mesh(mesh_target, nThresh, diagThresh);
    //tree.RecomputeBoundingBoxesUsingExistingCovFrames();      //*** is this ever needed?
    printf("Tree built: NNodes=%d  NData=%d  TreeDepth=%d\n", pTree->NumNodes(), pTree->NumData(), pTree->TreeDepth());
  }
  else
  {
    // build Point Cloud PD tree from mesh
    // (uses triangle center points as the point cloud)
    //  Note: the mesh constructor for the point cloud PD tree
    //        assumes the specified noise in direction perpendicular to surface
    //        and sets the in-plane noise based on the triangle size in order
    //        to allow for greater freedom of match anywhere along the triangle surface
    //        even though each triangle surface is represented by only one point.
    printf("Building point cloud PD tree .... ");
    cisstPointCloud pointCloud(mesh_target);
    DirPDTree_PointCloud *pPointCloudTree;
    pPointCloudTree = new DirPDTree_PointCloud(pointCloud, nThresh, diagThresh);
    pTree = pPointCloudTree;
    //tree.RecomputeBoundingBoxesUsingExistingCovFrames();      //*** is this ever needed?
    printf("Tree built: NNodes=%d  NData=%d  TreeDepth=%d\n", pTree->NumNodes(), pTree->NumData(), pTree->TreeDepth());
  }

  // Random Numbers: Normal RV's
  std::string normRVFile = workingDir + "GaussianValues.txt";
  std::ifstream randnStream(normRVFile.c_str());  // streams N(0,1) RV's

  // Generate random samples from mesh
  GenerateSamples(mesh_source, randSeed1, randSeqPos1, nSamples,
    samples, sampleNorms, sampleDatums,
    &saveSamplesPath);

  // Add noise to samples
  GenerateSampleSurfaceNoise(randSeed1, randSeqPos1, randnStream,
    sampleNoiseInPlane, sampleNoisePerpPlane,
    sampleNoiseCircSDDeg*cmnPI / 180.0, sampleNoiseEccentricity,
    samples, sampleNorms,
    noisySamples, noisySampleNorms,
    sampleNoiseCov, sampleNoiseInvCov,
    sampleNoiseL,
    percentOutliers,
    minPosOffsetOutlier, maxPosOffsetOutlier,
    minAngOffsetOutlier, maxAngOffsetOutlier,
    &saveNoisySamplesPath,
    &savePath_Cov, &savePath_L);

  // Generate random initial offset
  vctFrm3 Fi;
  GenerateRandomTransform(randSeed2, randSeqPos2,
    minOffsetPos, maxOffsetPos,
    minOffsetAng, maxOffsetAng,
    Fi);
#else
  // Replay Randomized Trial
  vctFrm3 Fi;
  //std::string baseFolder = "..\\ICP_TestData\\RandomTests\\";
  std::string baseFolder = "//jhdfs/data/lcsr$/CIIS/sbillin3/Research/ICOP/Data/RandomTests_FemurPatch/TestD_BivariateNormalTheta_OffsetR10-20P10-20_NoOutlierDetection/";  
  std::string testDir = "Noise_Plane1Norm1_Offset_R10-20T10-20_NSamps50/";
  std::string trialName = "IMLOP";     // Note: Finit and samples are same for all trials
  std::string trialNum = "19";
  std::string commonDir = "CommonFiles/";  
  std::string loadTargetMeshFile = baseFolder + commonDir + "SaveMeshTarget.mesh";
  std::string loadSourceMeshFile = baseFolder + commonDir + "SaveMeshSource.mesh";
  std::string loadNoisySamplesFile = baseFolder + testDir + commonDir + "SaveNoisySamples_" + trialNum + ".pts";  
  // need this to get FGuess
  std::string loadTrackPathFile = baseFolder + testDir + trialName + "/SaveTrackRegPath_" + trialName + "_" + trialNum + ".txt";

  //std::string workingFolder = "..\\ICP_TestData\\LastRun_CovEst\\";
  //std::string loadMeshFile = workingFolder + "SaveMesh.mesh";
  //std::string loadSamplesFile = workingFolder + "SaveSamples" + ".pts";
  //std::string loadNoisySamplesFile = workingFolder + "SaveNoisySamples" + ".pts";
  //std::string loadTrackPathFile = workingFolder + "SaveTrackRegPath" + ".txt";

  mesh_target.LoadMeshFile( loadTargetMeshFile );
  cisstICP::MeshSave(mesh_target, saveTargetMeshPath);
  if (TargetShapeAsMesh)
  { // target shape is a mesh
    pTree = new DirPDTree_Mesh(mesh_target, nThresh, diagThresh);
  }
  else
  { // target shape is a point cloud
    pTree = new DirPDTree_PointCloud(mesh_target, nThresh, diagThresh);
  }
  printf("Tree built: NNodes=%d  NData=%d  TreeDepth=%d\n", pTree->NumNodes(), pTree->NumData(), pTree->TreeDepth());
  cisstICP::SamplesLoad( noisySamples, noisySampleNorms, loadNoisySamplesFile );
  cisstICP::SamplesSave( noisySamples, noisySampleNorms, saveNoisySamplesPath );
  std::ifstream fs(loadTrackPathFile.c_str());
  double temp;
  double r00,r01,r02,r10,r11,r12,r20,r21,r22,p1,p2,p3;
  fs >> temp >> r00 >> r01 >> r02 >> r10 >> r11 >> r12 >> r20 >> r21 >> r22 >> p1 >> p2 >> p3;
  vctRot3 Rm( 
    r00, r01, r02,
    r10, r11, r12,
    r20, r21, r22,
    VCT_NORMALIZE );
  vct3 pm( p1,p2,p3 );
  Fi.Assign(Rm,pm);

  //double meshstdDevPerpPlanePlane = 0.0;         // noise model for mesh
  //double meshStdDevInPlane = 0.0;
  //double noisePosSD[3] = {0.5, 0.5, 0.5};   // noise model for samples (std dev along each axis)
#endif

  // creating ICP solver
  cisstICP ICP;

#ifdef ENABLE_CALLBACKS
  // adding ICP callbacks
  std::vector<cisstICP::Callback> callbacks;
  //  callback: iteration file
  cisstICP::Callback iterCallback;
  iterCallback.cbFunc = Callback_SaveIterationsToFile_testICPNormals;
  std::stringstream iterFile;
  iterFile << workingDir << outputDir << "SaveIterations.txt";
  std::ofstream iterFileStream(iterFile.str().c_str());
  iterCallback.userData = (void*)(&iterFileStream);
  callbacks.push_back(iterCallback);
  //  callback: track path file
  cisstICP::Callback xfmCallback;
  xfmCallback.cbFunc = Callback_TrackRegPath_testICPNormals;
  std::stringstream trackPathFile;
  trackPathFile << workingDir << outputDir << "SaveTrackRegPath.txt";
  std::ofstream xfmFileStream(trackPathFile.str().c_str());
  xfmCallback.userData = (void*)(&xfmFileStream);
  callbacks.push_back(xfmCallback);
#endif

  std::cout << std::endl << "Applying Sample Offset Fi: " << std::endl << Fi << std::endl;
  vctFrm3 FGuess = Fi;


  // ICP Algorithm
  algDirICP *pICPAlg = NULL;
  switch (algType)
  {
  case DirAlgType_StdICP:
  {
    if (TargetShapeAsMesh)
    { // target shape is a mesh
      DirPDTree_Mesh *pTreeMesh = dynamic_cast<DirPDTree_Mesh*>(pTree);
      pICPAlg = new algDirICP_StdICP_Mesh(pTreeMesh, noisySamples, noisySampleNorms);
    }
    else
    { // target shape is a point cloud
      DirPDTree_PointCloud *pTreePointCloud = dynamic_cast<DirPDTree_PointCloud*>(pTree);
      pICPAlg = new algDirICP_StdICP_PointCloud(pTreePointCloud, noisySamples, noisySampleNorms);
    }
    break;
  }
  case DirAlgType_IMLOP:
  {
    if (!TargetShapeAsMesh)
    {
      std::cout << "ERROR: Currently only mesh target supported for IMLOP" << std::endl;
      assert(0);
    }
    DirPDTree_Mesh *pTreeMesh = dynamic_cast<DirPDTree_Mesh*>(pTree);
    //algDirICP_IMLOP_Mesh *pAlg = new algDirICP_IMLOP_Mesh( pTreeMesh, &ICP );
    //pAlg->k_init = 0.0;
    //pAlg->sigma2_init = 1.0;
    //pAlg->wRpos = 0.5;
    double k = 1.0 / (sampleNoiseCircSD*sampleNoiseCircSD);
    double sigma2 = sampleNoiseInPlane*sampleNoiseInPlane;
    double wRpos = 0.5;
    double kfactor = 1.0;
    bool dynamicParamEst = true;
    algDirICP_IMLOP_Mesh *pAlg = new algDirICP_IMLOP_Mesh(
      pTreeMesh, noisySamples, noisySampleNorms,
      k, sigma2, wRpos, kfactor, dynamicParamEst);
    pICPAlg = pAlg;
    break;
  }
  //case DirAlgType_GIMLOP:
  //{
  //  if (!TargetShapeAsMesh)
  //  {
  //    std::cout << "ERROR: Currently only mesh target supported for GIMLOP" << std::endl;
  //    assert(0);
  //  }
  //  DirPDTree_Mesh *pTreeMesh = dynamic_cast<DirPDTree_Mesh*>(pTree);
  //  // define GIMLOP parameters
  //  double k = 1.0 / (sampleNoiseCircSD*sampleNoiseCircSD);
  //  double B = sampleNoiseEccentricity*k / 2.0;
  //  std::cout << "k: " << k << " B: " << B << std::endl;
  //  vctDynamicVector<double> argK(nSamples, k);
  //  vctDynamicVector<double> argB(nSamples, B);
  //  //vctDynamicVector<double> argB( nSamples,0.0 );

  //  //double sigma2 = sampleNoiseInPlane*sampleNoiseInPlane;
  //  //vctDynamicVector<double> argSigma2( nSamples,sigma2 );
  //  //vctDynamicVector<vctFixedSizeMatrix<double,3,2>> argL( nSamples );
  //  //vct3 xProd, L1,L2;
  //  //vct3 xAxis( 1.0,0.0,0.0 );
  //  //vct3 yAxis( 0.0,1.0,0.0 );
  //  //for (unsigned int i=0; i<nSamples; i++)
  //  //{ // set argL as isotropic
  //  //  xProd = vctCrossProduct( noisySampleNorms(i),xAxis );
  //  //  if (xProd.Norm() < 0.01)
  //  //  {
  //  //    xProd = vctCrossProduct( noisySampleNorms(i),yAxis );
  //  //  }
  //  //  L1 = xProd.Normalized();
  //  //  L2 = vctCrossProduct(noisySampleNorms(i),L1).Normalized();
  //  //  argL(i).Column(0) = L1;
  //  //  argL(i).Column(1) = L2;
  //  //}
  //  // create algorithm
  //  algDirICP_GIMLOP_Mesh *pAlg = new algDirICP_GIMLOP_Mesh(
  //    pTreeMesh, noisySamples, noisySampleNorms,
  //    argK, argB, sampleNoiseL, sampleNoiseInvCov);
  //  //pAlg->SetNoiseModel(argK, argB, sampleNoiseL, sampleNoiseInvCov);
  //  pICPAlg = pAlg;
  //  break;
  //}
  case DirAlgType_PIMLOP:
  {
    if (!TargetShapeAsMesh)
    {
      std::cout << "ERROR: Currently only mesh target supported for PIMLOP" << std::endl;
      assert(0);
    }
    DirPDTree_Mesh *pTreeMesh = dynamic_cast<DirPDTree_Mesh*>(pTree);
    // define noise model parameters
    vctDynamicVector<vct2> Xpln(nSamples, vct2(0.0));
    double k = 0.0;
    vct3x3 M = vct3x3::Eye();
    vctRot3 Rx_pln;
    vctDynamicVector<double> argK(nSamples, k);
    vctDynamicVector<vct3x3> argM(nSamples, M);
    vctDynamicVector<vctRot3> argRx_pln(nSamples, Rx_pln);

    // create algorithm
    algDirICP_PIMLOP_Mesh *pAlg = new algDirICP_PIMLOP_Mesh(
      pTreeMesh, 
      noisySamples, Xpln, argRx_pln, argK, argM );
    pICPAlg = pAlg;
    break;
  }
  default:
  {
    std::cout << "ERROR: unknown algorithm type" << std::endl;
    assert(0);
  }
  }

  // ICP Options
  //cisstICP::OptionsNormals opt;
  cisstICP::Options opt;
  opt.auxOutputDir = workingDir + outputDir;
  opt.maxIter = 100;
  opt.termHoldIter = 2;
  opt.minE = -std::numeric_limits<double>::max();
  opt.tolE = 0.0;
  opt.dPosThresh = 0.25;
  opt.dAngThresh = 0.25*(cmnPI / 180);
  opt.dPosTerm = 0.001;
  opt.dAngTerm = 0.001*(cmnPI / 180);
  //opt.alg = pICPAlg;

  // Run ICP
  int numRuns = 1;
  vctFrm3 Freg;
  double runtime = 0.0;
  cisstICP::ReturnType rv;
  for (int i = 0; i < numRuns; i++)
  {
    rv = ICP.RunICP(pICPAlg, opt, FGuess);
    //rv = ICP.RunICP_Rigid( samples, sampleNorms, *pTree, opt, FGuess, Freg );
    std::cout << rv.termMsg;
#ifdef ENABLE_CALLBACKS
    iterFileStream << rv.termMsg;
#endif
    runtime += rv.runTime;
    Freg = rv.Freg;
  }
  std::cout << std::endl << " ===> Avg RunTime: " << runtime / numRuns << std::endl;

  // Freg now includes Fi as FGuess => Freg should be identity for perfect registration
  vctFrm3 Ferr = Freg;
  vctRodRot3 Rerr(Ferr.Rotation());
  double terr = Ferr.Translation().Norm();
  double rerr = Rerr.Norm();
  vctRodRot3 Rinit(Fi.Rotation());
  double tinit = Fi.Translation().Norm();
  double rinit = Rinit.Norm();

  std::stringstream resultStream;
  resultStream << std::endl;
  resultStream << "Starting Offset:   \tdAng: " << rinit * 180 / cmnPI << "\tdPos: " << tinit << std::endl;
  resultStream << "Registration Error:\tdAng: " << rerr * 180 / cmnPI << "\tdPos: " << terr << std::endl << std::endl;
  std::cout << resultStream.str();
#ifdef ENABLE_CALLBACKS
  iterFileStream << resultStream.str();
#endif

  if (pICPAlg) delete pICPAlg;

  //std::cout << "press key then enter to quit:" << std::endl;
  //std::cin >> k;
}

#endif // _testICPNormals_H
