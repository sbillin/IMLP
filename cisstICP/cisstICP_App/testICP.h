#ifndef _testICP_H
#define _testICP_H

#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <limits.h>

#include "utility.h"
#include "cisstICP.h"
#include "cisstMesh.h"
#include "cisstPointCloud.h"
#include "PDTree_Mesh.h"
#include "PDTree_PointCloud.h"

#include "algICP_StdICP_Mesh.h"
#include "algICP_IMLP_Mesh.h"

#include "algICP_StdICP_PointCloud.h"
#include "algICP_IMLP_PointCloud.h"

enum ICPAlgType { AlgType_StdICP, AlgType_IMLP };

void Callback_TrackRegPath_testICP( cisstICP::CallbackArg &arg, void *userData )
{
  // Save to file:
  //  - error function
  //  - incremental transform
  // output format:
  //  error r00 r01 r02 r10 r11 r12 r20 r21 r22 tx ty tz
  std::ofstream *fs = (std::ofstream *)(userData);
  (*fs) << arg.E << " " << arg.dF.Rotation().Row(0) << " " << arg.dF.Rotation().Row(1) << " " 
    << " " << arg.dF.Rotation().Row(2) << " " << arg.dF.Translation() << std::endl;
}
void Callback_SaveIterationsToFile_testICP( cisstICP::CallbackArg &arg, void *userData )
{
  std::ofstream *fs = (std::ofstream *)(userData);

  vctRodRot3 dR(arg.dF.Rotation());
  std::stringstream ss;
  ss << cmnPrintf("iter=%u  E=%.3f  tolE=%.4f (dAng/dPos)= %.2f/%.2f  t=%.3f NNodes=%u/%u/%u NOut=%u") 
    << arg.iter 
    << arg.E 
    << arg.tolE
    << dR.Norm()*180/cmnPI << arg.dF.Translation().Norm()
    << arg.time 
    //<< arg.maxNodesSearched << arg.avgNodesSearched << arg.minNodesSearched
    << arg.nOutliers;

  (*fs) << ss.str() << std::endl;
}

// TargetShapeAsMesh    true - uses mesh to represent target shape
//                      false - uses point cloud (taken from mesh) to represent target shape
void testICP(bool TargetShapeAsMesh, ICPAlgType algType)
{
  //char k;
  //std::cout << "press key then enter to start:" << std::endl;
  //std::cin >> k;

  int    nThresh = 5;       // Cov Tree Params
  double diagThresh = 5.0;  //  ''

  std::string workingDir = "../test_data/";
  std::string outputDir =  "LastRun/";

  std::string saveMeshPath =          workingDir + outputDir + "SaveMesh";
  std::string saveSamplesPath =       workingDir + outputDir + "SaveSamples.pts";
  std::string saveNoisySamplesPath =  workingDir + outputDir + "SaveNoisySamples.pts";
  std::string saveCovPath =           workingDir + outputDir + "SaveSampleCov.txt";
  std::string saveOffsetXfmPath =     workingDir + outputDir + "SaveOffsetXfm.txt";
  std::string saveRegXfmPath =        workingDir + outputDir + "SaveRegXfm.txt";
  //std::string saveLPath =             workingDir + outputDir + "SaveSampleL";

  cisstMesh         mesh;
  PDTreeBase*         pTree;
  vctDynamicVector<vct3>    samples;
  vctDynamicVector<vct3>    sampleNorms;
  vctDynamicVector<vct3>    noisySamples;
  vctDynamicVector<vct3>    noisySampleNorms;
  vctDynamicVector<unsigned int>  sampleDatums;
  vctDynamicVector<vct3x3>  sampleNoiseCov;
  vctDynamicVector<vct3x3>  sampleNoiseInvCov;
  vctDynamicVector<vct3x2>  sampleNoiseL;

#if 1

  std::string loadMeshPath = workingDir + "ProximalFemur.ply";
  //std::string loadMeshPath = workingDir + "RIGHTHEMIPELVIS_centered.mesh";
  //std::string loadMeshPath = workingDir + "RIGHTHEMIPELVIS.mesh";
  //std::string loadMeshPath = workingDir + "CTBreastImage_Dec20000_Shell.mesh";  

  const int nSamples = 100;

  // Samples Noise Model
  //  NOTE: this is a generative noise model (i.e. noise is generated according
  //        to the noise properties defined here)
  //double noiseSampsSD[3] = {1.0, 1.0, 1.0};   // noise model for samples (std dev along each axis)
  double sampleNoiseInPlane = 1.0;      // standard deviation of noise in and out of plane
  double sampleNoisePerpPlane = 2.0;    //   ''

  // Target Noise Model (for point cloud target only)
  //  NOTE: this is a descriptive model, not a generative one
  //        i.e. no noise is added to the point cloud, the noise model is merely
  //        allow for errors at intermediate locations between the points and penalize
  //        errors offset from the surface
  double PointCloudNoisePerpPlane = 0.5;  // noise model for point cloud using mesh constructor
                                          //  Note: in-plane noise set automatically relative to triangle size

  double minOffsetPos = 50.0;
  double maxOffsetPos = 100.0;
  double minOffsetAng = 30.0;
  double maxOffsetAng = 60.0;
  
  double percentOutliers = 0.05;
  double minPosOffsetOutlier = 5.0;
  double maxPosOffsetOutlier = 10.0;
  double minAngOffsetOutlier = 0.0;
  double maxAngOffsetOutlier = 0.0;

  unsigned int randSeed1 = 0;       // generates samples
  unsigned int randSeqPos1 = 0;
  unsigned int randSeed2 = 17;      // generates offsets
  unsigned int randSeqPos2 = 28;

  // load mesh
  CreateMesh( mesh, loadMeshPath, &saveMeshPath );

  // Create target shape from mesh (as a PD tree)
  if (TargetShapeAsMesh)
  {
    // build PD tree on the mesh directly
    //  Note: defines measurement noise to be zero
    printf("Building mesh PD tree .... ");
    pTree = new PDTree_Mesh(mesh, nThresh, diagThresh);
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
    cisstPointCloud pointCloud(mesh, PointCloudNoisePerpPlane);
    PDTree_PointCloud *pPointCloudTree;
    pPointCloudTree = new PDTree_PointCloud(pointCloud, nThresh, diagThresh);
    pTree = pPointCloudTree;
    //tree.RecomputeBoundingBoxesUsingExistingCovFrames();      //*** is this ever needed?
    printf("Tree built: NNodes=%d  NData=%d  TreeDepth=%d\n", pTree->NumNodes(), pTree->NumData(), pTree->TreeDepth());
    //printf(" Point Cloud Noise Model:\n  perp-plane variance = %f\n  in-plane variance = %f (avg)\n\n", 
    //  PointCloudNoisePerpPlane, pPointCloudTree->avgVarInPlane);
  }

  // Random Numbers: Normal RV's
  std::string normRVFile = workingDir + "GaussianValues.txt";
  std::ifstream randnStream(normRVFile.c_str());  // streams N(0,1) RV's

  //// initialize random numbers
  //cmnRandomSequence &cisstRandomSeq = cmnRandomSequence::GetInstance();
  //cisstRandomSeq.SetSeed( randSeed1 );
  //cisstRandomSeq.SetSequencePosition( randSeqPos1 );

  // Generate random samples from mesh
  GenerateSamples(mesh, randSeed1, randSeqPos1, nSamples,
                   samples, sampleNorms, sampleDatums,
                   &saveSamplesPath );

  // Add noise to samples
  GenerateSampleSurfaceNoise(randSeed1, randSeqPos1, randnStream,
      sampleNoiseInPlane, sampleNoisePerpPlane, 0.0, 0.0,
      samples, sampleNorms,
      noisySamples, noisySampleNorms,
      sampleNoiseCov, sampleNoiseInvCov, sampleNoiseL,
      percentOutliers,
      minPosOffsetOutlier, maxPosOffsetOutlier,
      minAngOffsetOutlier, maxAngOffsetOutlier,
      &saveNoisySamplesPath,
      &saveCovPath);

  //// initialize random numbers
  //cisstRandomSeq.SetSeed( randSeed2 );
  //cisstRandomSeq.SetSequencePosition( randSeqPos2 );

  // Generate random initial offset
  vctFrm3 Fi;
  GenerateRandomTransform( randSeed2, randSeqPos2,
                           minOffsetPos, maxOffsetPos,
                           minOffsetAng, maxOffsetAng,
                           Fi );
  // save initial offset
  transform_write(Fi, saveOffsetXfmPath);

#else
  // Replay Randomized Trial
  vctFrm3 Fi;
  std::string baseFolder = "..\\ICP_TestData\\RandomTests\\";
  std::string noiseDir = "NoiseIP0.5PP0.5_OffsetR0-10P0-10\\";
  std::string trialName = "IMLP";     // this can be any of the trials (since Finit and samples are same for all trials)
  std::string trialNum = "1";
  std::string commonDir = "CommonFiles\\";
  std::string loadMeshFile = baseFolder + commonDir + "SaveMesh.mesh";
  //std::string loadSamplesFile = baseFolder + noiseLevel + commonDir + "SaveSamples_" + trialNum + ".pts";
  std::string loadNoisySamplesFile = baseFolder + noiseDir + commonDir + "SaveNoisySamples_" + trialNum + ".pts";
  // need this to get FGuess
  std::string loadTrackPathFile = baseFolder + noiseDir + trialName + "\\SaveTrackRegPath_" + trialName + "_" + trialNum + ".txt";

  //std::string workingFolder = "..\\ICP_TestData\\LastRun_CovEst\\";
  //std::string loadMeshFile = workingFolder + "SaveMesh.mesh";
  //std::string loadSamplesFile = workingFolder + "SaveSamples" + ".pts";
  //std::string loadNoisySamplesFile = workingFolder + "SaveNoisySamples" + ".pts";
  //std::string loadTrackPathFile = workingFolder + "SaveTrackRegPath" + ".txt";

  mesh.LoadMeshFile( loadMeshFile );
  cisstICP::MeshSave(mesh, saveMeshFile);
  pTree = new PDTree_PointCloud(mesh, nThresh, diagThresh);
  printf("Tree built: NNodes=%d  NData=%d  TreeDepth=%d\n", pTree->NumNodes(), pTree->NumData(), pTree->TreeDepth());
  cisstICP::SamplesLoad( noisySamples, loadNoisySamplesFile );
  cisstICP::SamplesSave( noisySamples, saveNoisySamplesFile );
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

  std::cout << Fi << std::endl;

  double meshstdDevPerpPlanePlane = 0.0;         // noise model for mesh
  double meshStdDevInPlane = 0.0;
  double noisePosSD[3] = {0.5, 0.5, 0.5};   // noise model for samples (std dev along each axis)
#endif


  // creating ICP solver
  cisstICP ICP;

  // user callbacks for ICP
  std::vector<cisstICP::Callback> userCallbacks;
  //  callback: iteration file
  cisstICP::Callback iterCallback;
  iterCallback.cbFunc = Callback_SaveIterationsToFile_testICP;
  std::stringstream iterFile;
  iterFile << workingDir << outputDir << "SaveIterations.txt";
  std::ofstream iterFileStream(iterFile.str().c_str());
  iterCallback.userData = (void*)(&iterFileStream);
  userCallbacks.push_back(iterCallback);
  //  callback: track path file
  cisstICP::Callback xfmCallback;
  xfmCallback.cbFunc = Callback_TrackRegPath_testICP;
  std::stringstream trackPathFile;
  trackPathFile << workingDir << outputDir << "SaveTrackRegPath.txt";
  std::ofstream xfmFileStream(trackPathFile.str().c_str());
  xfmCallback.userData = (void*)(&xfmFileStream);
  userCallbacks.push_back( xfmCallback );

  std::cout << std::endl << "Applying Sample Offset Fi: " << std::endl << Fi << std::endl;
  vctFrm3 FGuess = Fi;


  // ICP Algorithm
  algICP *pICPAlg = NULL;
  switch (algType)
  {
  case AlgType_StdICP:
    {
      if (TargetShapeAsMesh)
      { // target shape is a mesh
        PDTree_Mesh *pTreeMesh = dynamic_cast<PDTree_Mesh*>(pTree);
        pICPAlg = new algICP_StdICP_Mesh(pTreeMesh, noisySamples);
      }
      else
      { // target shape is a point cloud
        PDTree_PointCloud *pTreePointCloud = dynamic_cast<PDTree_PointCloud*>(pTree);
        pICPAlg = new algICP_StdICP_PointCloud(pTreePointCloud, noisySamples);
      }
      break;
    }
  case AlgType_IMLP:
    {
      if (TargetShapeAsMesh)
      { // target shape is a mesh
        PDTree_Mesh *pTreeMesh = dynamic_cast<PDTree_Mesh*>(pTree);
        algICP_IMLP_Mesh *pAlg;
        pAlg = new algICP_IMLP_Mesh(pTreeMesh, noisySamples, sampleNoiseCov, sampleNoiseCov);
        // set sample noise model
        //pAlg->SetSampleCovariances( sampleNoiseCov );
        // set mesh noise model to zero noise
        mesh.TriangleCov.SetSize( mesh.NumTriangles() );
        mesh.TriangleCovEig.SetSize( mesh.NumTriangles() );
        mesh.TriangleCov.SetAll( vct3x3(0.0) );
        mesh.TriangleCovEig.SetAll( vct3(0.0) );
        pTreeMesh->ComputeNodeNoiseModels();
        //double noiseSDInPlane = 0.5;
        //double noiseSDPerpPlane = 1.0;
        //SetMeshTriangleCovariances( mesh, noiseSDInPlane, noiseSDPerpPlane );
        pICPAlg = pAlg;
      }
      else
      { // target shape is a point cloud
        PDTree_PointCloud *pTreePointCloud = dynamic_cast<PDTree_PointCloud*>(pTree);
        algICP_IMLP_PointCloud *pAlg;
        pAlg = new algICP_IMLP_PointCloud(pTreePointCloud, noisySamples, sampleNoiseCov, sampleNoiseCov);
        // set sample noise model
        //pAlg->SetSampleCovariances( sampleNoiseCov );
        pICPAlg = pAlg;
        // NOTE: (PD tree noise model was already defined
        //       by the point cloud cov tree constructor)
      }
      break;
    }
  default:
    {
      std::cout << "ERROR: unknown algorithm type" << std::endl;
      assert(0);
    }
  }

  // ICP Options
  cisstICP::Options opt;
  opt.auxOutputDir = workingDir + outputDir;
  opt.maxIter = 100;
  opt.termHoldIter = 2;
  opt.minE = -std::numeric_limits<double>::max();
  opt.tolE = 0.0;
  opt.dPosThresh = 0.1;
  opt.dAngThresh = 0.1*(cmnPI/180);
  opt.dPosTerm = 0.01;
  opt.dAngTerm = 0.01*(cmnPI/180);

  // Run ICP
  int numRuns = 1;
  vctFrm3 Freg;
  double runtime = 0.0;
  cisstICP::ReturnType rv;
  for (int i=0; i<numRuns; i++)
  {
    rv = ICP.RunICP(pICPAlg, opt, FGuess, &userCallbacks);
    //rv = ICP.RunICP_Rigid( samples, *pTree, opt, FGuess, Freg );
    std::cout << rv.termMsg;
    iterFileStream << rv.termMsg;
    runtime += rv.runTime;
    Freg = rv.Freg;
  }
  std::cout << std::endl << " ===> Avg RunTime: " << runtime / numRuns << std::endl;

  // save registration result
  transform_write(Freg, saveRegXfmPath);

  // Freg now includes Fi as FGuess => Freg should be identity for perfect registration
  vctFrm3 Ferr = Freg;    
  //vctFrm3 Ferr = Fi * Freg;
  vctRodRot3 Rerr(Ferr.Rotation());
  double terr = Ferr.Translation().Norm();
  double rerr = Rerr.Norm();
  vctRodRot3 Rinit(Fi.Rotation());
  double tinit = Fi.Translation().Norm();
  double rinit = Rinit.Norm();

  std::stringstream resultStream;
  resultStream << std::endl;
  resultStream << "Starting Offset:   \tdAng: " << rinit*180/cmnPI << "\tdPos: " << tinit << std::endl;
  resultStream << "Registration Error:\tdAng: " << rerr*180/cmnPI << "\tdPos: " << terr << std::endl << std::endl;
  std::cout << resultStream.str();
  iterFileStream << resultStream.str();

  if (pICPAlg) delete pICPAlg;

  //std::cout << "press key then enter to quit:" << std::endl;
  //std::cin >> k;
}

#endif // _testICP_H
