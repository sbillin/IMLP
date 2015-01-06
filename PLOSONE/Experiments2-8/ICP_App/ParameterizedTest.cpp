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
#include "ParameterizedTest.h"
#include "cisstICP.h"
#include "cisstPointCloud.h"
#include "utility.h"

//#include "cpd.h"
#include "cisstCovTree_PointCloud.h"

#include <windows.h>

// place anonymous namespace around global
//  functions to avoid collision with other cpp files
namespace
{  // anonymous namspace

  //--- cisstICP Callbacks ---//

  void Callback_TrackRegPath(cisstICP::CallbackArg &arg, void *userData)
  {
    // Save to file:
    //  - error function
    //  - incremental transform
    // output format:
    //  error r00 r01 r02 r10 r11 r12 r20 r21 r22 tx ty tz nOutliers
    std::ofstream *fs = (std::ofstream *)(userData);
    (*fs) << arg.E << "    " << arg.dF.Rotation().Row(0) << " " << arg.dF.Rotation().Row(1) << " "
      << " " << arg.dF.Rotation().Row(2) << " " << arg.dF.Translation() << "   " << arg.nOutliers << std::endl;
  }

  void Callback_SaveIterationsToFile(cisstICP::CallbackArg &arg, void *userData)
  {
    std::ofstream *fs = (std::ofstream *)(userData);
    vctRodRot3 dR(arg.dF.Rotation());
    std::stringstream ss;
    ss << cmnPrintf("iter=%u  E=%.3f  tolE=%.4f (dAng/dPos)= %.2f/%.2f  t=%.3f NOut=%u")  //NNodes = %u / %u / %u NOut = %u")
      << arg.iter
      << arg.E
      << arg.tolE
      << dR.Norm() * 180 / cmnPI << arg.dF.Translation().Norm()
      << arg.time
      //<< arg.maxNodesSearched << arg.avgNodesSearched << arg.minNodesSearched
      << arg.nOutliers;
    (*fs) << ss.str();
    (*fs) << std::endl;
  }

  ////--- CPD Callbacks ---//

  //void Callback_TrackRegPath_CPD(CPD::IterationData &arg, void *userData)
  //{
  //  // Save to file:
  //  //  - error function
  //  //  - incremental transform
  //  // output format:
  //  //  error r00 r01 r02 r10 r11 r12 r20 r21 r22 tx ty tz
  //  std::ofstream *fs = (std::ofstream *)(userData);
  //  (*fs) << arg.E << " " << arg.dF.Rotation().Row(0) << " " << arg.dF.Rotation().Row(1) << " "
  //    << " " << arg.dF.Rotation().Row(2) << " " << arg.dF.Translation() << std::endl;
  //}

  //void Callback_SaveIterationsToFile_CPD(CPD::IterationData &arg, void *userData)
  //{
  //  std::ofstream *fs = (std::ofstream *)(userData);
  //  vctRodRot3 dR(arg.dF.Rotation());
  //  std::stringstream ss;
  //  ss << cmnPrintf("iter=%u  E=%.3f  tolE=%.4f (dAng/dPos)= %.2f/%.2f  sigma2=%f")  //t=%.3f  NOut=%u  NNodes = %u / %u / %u NOut = %u")
  //    << arg.iter
  //    << arg.E
  //    << arg.dL
  //    << dR.Norm() * 180 / cmnPI << arg.dF.Translation().Norm()
  //    << arg.sigma2
  //    ;
  //    //<< arg.time
  //    //<< arg.maxNodesSearched << arg.avgNodesSearched << arg.minNodesSearched
  //    //<< arg.nOutliers;
  //  (*fs) << ss.str();
  //  (*fs) << std::endl;
  //}

  // compute ground truth error statistics
  // Note: this routine assumes that true registration is identity
  //       => designed for trials where FGuess is used to offset
  //       the two shapes to be registered
  void ComputeGroundTruthStats(
    const vctDynamicVector<vct3> &samples,
    const vctFrm3 &Freg,
    double &posErr_mean, double &posErr_SD)

  {
    unsigned int nSamps = samples.size();

    vctDynamicVector<double> posErr(nSamps);
    posErr_mean = 0.0;
    vct3 t1;
    double t2;
    for (unsigned int i = 0; i < nSamps; i++)
    { // compute mean distances
      t1 = Freg*samples[i] - samples[i];
      posErr[i] = t1.Norm();
      posErr_mean += posErr[i];
    }
    posErr_mean /= nSamps;

    double posErr_Var = 0.0;
    for (unsigned int i = 0; i < nSamps; i++)
    { // compute standard deviation of distance error
      t2 = (posErr[i] - posErr_mean);
      posErr_Var += t2*t2;
    }
    posErr_Var /= nSamps;
    posErr_SD = sqrt(posErr_Var);
  }

} // namspace anonymous



// Parameterized trial for position based registration
//   NOTE: noise model over mesh must be initialized before calling this function
void Run_ParameterizedTest(
  cisstMesh &mesh, cisstCovTreeBase *pTree,
  TestParameters params)
{
  int rv;
  cisstICP ICP;
  //CPD cpd;

  std::cout << params.toString() << std::endl;

  // Create folders for output files
  std::string inputDir = params.inputDir;
  std::string dataOutputDir = params.outputDataDir;
  std::string commonOutputDir = params.outputCommonDir;

#ifndef DISABLE_OUTPUT_FILES
  CreateDir(dataOutputDir);
  CreateDir(commonOutputDir);

  // Registration Settings File
  //  the line below was creating an error due to file names becoming too long
  //std::string settingsFile = dataOutputDir + "/RegistrationSettings_" + params.AlgorithmDescription + ".txt";
  std::string settingsFile = dataOutputDir + "/RegistrationSettings.txt";
  std::ofstream settingsStream(settingsFile.c_str());
  assert(!settingsStream.fail());
  settingsStream << params.toString() << std::endl;
  settingsStream.close();

  // stats files (output)
  std::string statsFile = dataOutputDir + "/stats_" + params.AlgorithmDescription + ".txt";
  std::ofstream statsStream(statsFile.c_str());
  assert(!statsStream.fail());
#endif


  // random variable generators
  // normal distributed random variables input stream
  //  supplies standard normal random variables N(0,1)
  std::string normRVFile = inputDir + "/GaussianValues.txt";
  std::ifstream randnStream_SampleNoise(normRVFile.c_str());
  // uniform random number generator is a singleton class
  //  => need to keep track of sequence positions
  //  if using it for multiple purposes and if wanting to keep same samples
  //  across multiple tests
  cmnRandomSequence &cisstRandomSeq = cmnRandomSequence::GetInstance();
  unsigned int randSeqPos_Samples = 0;
  unsigned int randSeqPos_Xfms = 0;
  unsigned int randSeqPos_Validation = 0;

  vctDynamicVector<vct3>  samples;
  vctDynamicVector<vct3>  sampleNorms;
  vctDynamicVector<vct3>  noisySamples;
  vctDynamicVector<vct3>  noisySampleNorms;
  vctDynamicVector<vct3>  validationSamples;
  vctDynamicVector<vct3>  validationSampleNorms;
  vctDynamicVector<unsigned int>  validationSampleDatums;

  // Generate validation samples
  std::string saveValidationSamplesPath = commonOutputDir + "/SaveValidationSamples.pts";
  std::string *pSaveValidationSamplesPath = &saveValidationSamplesPath;
#ifdef DISABLE_OUTPUT_FILES
  pSaveValidationSamplesPath = NULL;
#endif
  GenerateSamples(
    mesh,
    params.randSeed_Validation, randSeqPos_Validation,
    params.nValidationSamples,
    validationSamples, validationSampleNorms, validationSampleDatums,
    pSaveValidationSamplesPath);


  //=== Run Trials ===//

  for (unsigned int i = 0; i < params.numTrials; i++)
  {
    unsigned int trialNum = i + 1;

    // Generate samples
    params.GenerateSamplesFunc(
      mesh, params,
      params.randSeed_Samples, randSeqPos_Samples, randnStream_SampleNoise,
      samples, sampleNorms, noisySamples,
      trialNum);

    // Generate random misalignment
    vctFrm3 Fi;
    GenerateRandomTransform(
      params.randSeed_Xfms, randSeqPos_Xfms,
      params.minOffsetPos, params.maxOffsetPos,
      params.minOffsetAng, params.maxOffsetAng,
      Fi);

#ifndef DISABLE_OUTPUT_FILES
    // save random transform
    std::stringstream saveXfm;
    saveXfm << commonOutputDir << "/SaveXfm_" << trialNum << ".tfm";
    transform_write(Fi, saveXfm.str());
#endif

    vctFrm3 FGuess = Fi;
    vctFrm3 Freg;
    double runTime, runTimeFirstMatch;
    unsigned int numIter;
    unsigned int nOutliers;
    double matchErrorAvg, matchErrorStdDev;

#ifndef DISABLE_OUTPUT_FILES
    std::stringstream iterFile;
    //iterFile << dataOutputDir << "/SaveIterations_" << params.AlgorithmDescription << "_" << trialNum << ".txt";
    iterFile << dataOutputDir << "/SaveIterations_" << trialNum << ".txt";
    std::ofstream iterFileStream(iterFile.str().c_str());
    assert(!iterFileStream.fail());

    std::stringstream trackPathFile;
    //trackPathFile << dataOutputDir << "/SaveTrackRegPath_" << params.AlgorithmDescription << "_" << trialNum << ".txt";
    trackPathFile << dataOutputDir << "/SaveTrackRegPath_" << trialNum << ".txt";
    std::ofstream xfmFileStream(trackPathFile.str().c_str());
    assert(!xfmFileStream.fail());
#endif

    std::stringstream initStream;
    initStream << std::endl << std::endl << "Starting Randomized Trial: " << params.AlgorithmDescription
      << " (Trial #" << trialNum << ")" << std::endl;
    initStream << "Applying Random Offset: " << std::endl << Fi << std::endl;
    std::cout << initStream.str();

#ifndef DISABLE_OUTPUT_FILES
    iterFileStream << initStream.str();
#endif

    switch (params.algType)
    {
    case TestParameters::StdICP:
    case TestParameters::IMLP:
    case TestParameters::RobustICP:
    {
      // user callbacks for ICP
      std::vector<cisstICP::Callback> userCallbacks;
#ifndef DISABLE_OUTPUT_FILES
      // callback: iteration file
      userCallbacks.push_back(cisstICP::Callback(Callback_SaveIterationsToFile, (void*)(&iterFileStream)));
      // callback: track path file
      userCallbacks.push_back(cisstICP::Callback(Callback_TrackRegPath, (void*)(&xfmFileStream)));
#endif

      // user callback enable/disable
      std::vector<cisstICP::Callback> *pUserCallbacks;
      if (params.bEnableIterationCallbacks)
      {
        pUserCallbacks = &userCallbacks;
      }
      else
      {
        pUserCallbacks = NULL;
      }

      // Run Registration    
      cisstICP::ReturnType rvICP;
      rvICP = ICP.RunICP(params.pAlg, params.opt_ICP, FGuess, pUserCallbacks, params.bEnableIterationCallbacks);
      Freg = rvICP.Freg;
      runTime = rvICP.runTime;
      runTimeFirstMatch = rvICP.runTimeFirstMatch;
      numIter = rvICP.numIter;
      nOutliers = rvICP.nOutliers;
      matchErrorAvg = rvICP.MatchPosErrAvg;
      matchErrorStdDev = rvICP.MatchPosErrSD;

      std::stringstream termStream;
      termStream << rvICP.termMsg;
#ifndef DISABLE_OUTPUT_FILES
      iterFileStream << termStream.str();
#endif
      std::cout << termStream.str();

      break;
    }

      //case TestParameters::CPD:
      //{
      //  // user callbacks for CPD
      //  std::vector<CPD::Callback> userCallbacks;
      //  // callback: iteration file
      //  userCallbacks.push_back(CPD::Callback(Callback_SaveIterationsToFile_CPD, (void*)(&iterFileStream)));
      //  // callback: track path file
      //  userCallbacks.push_back(CPD::Callback(Callback_TrackRegPath_CPD, (void*)(&xfmFileStream)));

      //  // user callback enable/disable
      //  std::vector<CPD::Callback> *pUserCallbacks;
      //  if (params.bEnableIterationCallbacks)
      //  {
      //    pUserCallbacks = &userCallbacks;
      //  }
      //  else
      //  {
      //    pUserCallbacks = NULL;
      //  }
      //  cpd.ClearIterationCallbacks();
      //  cpd.AddIterationCallbacks(userCallbacks);

      //  // Run Registration    
      //  CPD::ReturnType rvCPD;
      //  cisstCovTree_PointCloud *pTreePointCloud;
      //  pTreePointCloud = dynamic_cast<cisstCovTree_PointCloud*>(pTree);  // CPD is a point cloud only method
      //  rvCPD = cpd.cpd_register(noisySamples, pTreePointCloud->points, params.opt_CPD, FGuess);
      //  Freg = rvCPD.Freg;
      //  runTime = rvCPD.runTime;
      //  runTimeFirstMatch = 0.0;  // N/A
      //  numIter = rvCPD.numIter;

      //  std::stringstream termStream;
      //  termStream << rvCPD.termMsg;
      //  std::cout << termStream.str();
      //  iterFileStream << termStream.str();
      //  break;
      //}

    default:
    {
      std::cout << std::endl
        << "========> ERROR: Unrecognized Algorithm Type for Parameterized Test"
        << std::endl << std::endl;
      assert(0);
      return;
    }
    }

    // Registration Analysis 
    vctRodRot3 Rinit(Fi.Rotation());
    double tinit = Fi.Translation().Norm();
    double rinit = Rinit.Norm();

    vctFrm3 Ferr = Freg;    // Freg includes Finit as FGuess => Freg should be identity for perfect registration
    vctRodRot3 Rerr(Ferr.Rotation());
    double terr = Ferr.Translation().Norm();
    double rerr = Rerr.Norm();

    std::stringstream resultStream;
    resultStream << std::endl;
    resultStream << "Starting Offset:   dAng: " << rinit * 180 / cmnPI << "\tdPos: " << tinit << std::endl;
    resultStream << "Reg Error:         dAng: " << rerr * 180 / cmnPI << "\tdPos: " << terr << std::endl;
    resultStream << "Runtime: " << runTime << std::endl;
    std::cout << resultStream.str() << std::endl;
#ifndef DISABLE_OUTPUT_FILES
    iterFileStream << resultStream.str() << std::endl;
#endif

    // Summary Statistics
    //   TRE_mean TRE_SD ...
    //   F_AngErr F_PosErr runTime runTimeFirstMatch numIter nOutliers ...
    //   MatchError_mean MatchError_SD
    double TRE_mean, TRE_SD;
    // This routine assumes true registration is identity
    ComputeGroundTruthStats(
      validationSamples,
      Freg,
      TRE_mean, TRE_SD);

#ifndef DISABLE_OUTPUT_FILES
    statsStream
      << TRE_mean << " " << TRE_SD << "\t"
      << rerr*180.0 / cmnPI << " " << terr << "\t"
      << runTime << "\t" << runTimeFirstMatch << "\t"
      << numIter << "\t" << nOutliers << "\t"
      << matchErrorAvg << "\t" << matchErrorStdDev << std::endl;
#endif

  }
}

