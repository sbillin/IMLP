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
#include <sstream>
#include <limits>

#include <cisstOSAbstraction.h>

#include "cisstICP.h"
#include "algICP.h"

// debug
//#define ENABLE_CODE_TRACE
//#define ENABLE_CODE_PROFILER

cisstICP::ReturnType cisstICP::RunICP(
  algICP *pAlg,
  const Options &opt,
  const vctFrm3 &FGuess,
  std::vector<Callback> *pUserCallbacks,
  bool bEnableAlgorithmCallbacks)
{
  if (pAlg == NULL)
  {
    std::cout << "ERROR: no registration algorithm specified" << std::endl;
    assert(0);
    return ReturnType();
  }
  this->pAlgorithm = pAlg;
  this->opt = opt;
  this->FGuess = FGuess;

  // setup iteration callbacks
  ClearIterationCallbacks();
  if (pUserCallbacks)
  {
    // add user callbacks
    AddIterationCallbacks(*pUserCallbacks);
  }
  if (bEnableAlgorithmCallbacks)
  {
    // add algorithm callbacks
    // NOTE: for some reason, linux requires callbacks to be stored
    //       to a local variable before use as a function argument
    std::vector<Callback> callbacks = pAlg->ICP_GetIterationCallbacks();
    AddIterationCallbacks(callbacks);
  }

  // begin registration
  return IterateICP();
}


cisstICP::ReturnType cisstICP::IterateICP()
{
  //bool JustDidAccelStep = false;
  std::stringstream termMsg;
  double dAng, dPos;
  double dAng01, dAng12 = 0.0;
  double dPos01, dPos12 = 0.0;
  vctFrm3 dF;
  vctFrm3 F01, F12, F02;
  vctFrm3 Freg0, Freg1, Freg2;
  vctRodRot3 dR;
  double E0, E1, E2;
  double tolE;
  ReturnType rt;
  unsigned int terminateIter = 0;  // consecutive iterations satisfying termination

#ifdef ENABLE_CODE_PROFILER
  osaStopwatch codeProfiler;
  double time_Callbacks = 0.0;
  double time_Extras = 0.0;
  double time_Match = 0.0;
  double time_UpdateParams_PostMatch = 0.0;
  double time_FilterMatches = 0.0;
  double time_EvalErrorFunc = 0.0;
  double time_Register = 0.0;
  double time_UpdateParams_PostRegister = 0.0;
  codeProfiler.Reset();
  codeProfiler.Start();
#endif 

  if (opt.printOutput)
  {
    std::cout << "\n===================== Beginning Registration ==================\n";
  }

  osaStopwatch totalTimer;
  osaStopwatch iterTimer;
  totalTimer.Reset();
  totalTimer.Start();
  iterTimer.Reset();
  iterTimer.Start();

  //--- ICP Initialize ---//

  // initialize algorithm
  Freg0 = Freg1 = vctFrm3::Identity();
  dF = Freg2 = Freg = FGuess;
  pAlgorithm->ICP_InitializeParameters(FGuess);

#ifdef ENABLE_CODE_PROFILER
  time_Extras = codeProfiler.GetElapsedTime();
  codeProfiler.Reset();
  codeProfiler.Start();
#endif

  //------------ ICP Iterate ----------------//

  unsigned int iter;
  for (iter = 1; iter <= opt.maxIter; iter++)
  {

#ifdef ENABLE_CODE_TRACE
    std::cout << "ComputeMatches()" << std::endl;
#endif

    // compute matches
    pAlgorithm->ICP_ComputeMatches();

#ifdef ENABLE_CODE_PROFILER
    time_Match = codeProfiler.GetElapsedTime();
    codeProfiler.Reset();
    codeProfiler.Start();
#endif

#ifdef ENABLE_CODE_TRACE
    std::cout << "UpdateParameters_PostMatch()" << std::endl;
#endif

    // update algorithm's post-match parameters
    pAlgorithm->ICP_UpdateParameters_PostMatch();

#ifdef ENABLE_CODE_PROFILER
    time_UpdateParams_PostMatch = codeProfiler.GetElapsedTime();
    codeProfiler.Reset();
    codeProfiler.Start();
#endif

#ifdef ENABLE_CODE_TRACE
    std::cout << "FilterMatches()" << std::endl;
#endif

    // filter matches for outliers
    nOutliers = pAlgorithm->ICP_FilterMatches();

#ifdef ENABLE_CODE_PROFILER
    time_FilterMatches = codeProfiler.GetElapsedTime();
    codeProfiler.Reset();
    codeProfiler.Start();
#endif

    //--- First Iteration: report initial match statistics ---//
    if (iter == 1)
    {
#ifdef ENABLE_CODE_TRACE
      std::cout << "EvaluateErrorFunction()" << std::endl;
#endif

      // Compute initial error function value
      E = pAlgorithm->ICP_EvaluateErrorFunction();
      E0 = E1 = std::numeric_limits<double>::max();
      E2 = E;
      tolE = 0.0;
      Ebest = E;
      iterBest = 0;
      Fbest = FGuess;

#ifdef ENABLE_CODE_PROFILER
      time_EvalErrorFunc = codeProfiler.GetElapsedTime();
      codeProfiler.Reset();
      codeProfiler.Start();
#endif

#ifdef ENABLE_CODE_TRACE
      std::cout << "Callbacks" << std::endl;
#endif

      // initial callback
      iterData.iter = 0;
      iterData.E = E;
      iterData.tolE = tolE;
      iterData.Freg.Assign(FGuess);
      iterData.dF.Assign(FGuess);
      iterData.time = iterTimer.GetElapsedTime();
      iterData.nOutliers = nOutliers;
      //iterData.isAccelStep = false;
      std::vector<Callback>::iterator cbIter;
      for (cbIter = this->iterationCallbacks.begin(); cbIter != this->iterationCallbacks.end(); cbIter++)
      {
        cbIter->cbFunc(iterData, cbIter->userData);
      }      
      iterTimer.Reset();
      iterTimer.Start();

      rt.runTimeFirstMatch = iterData.time;

#ifdef ENABLE_CODE_PROFILER
      time_Callbacks = codeProfiler.GetElapsedTime();
      codeProfiler.Reset();
      codeProfiler.Start();
#endif

#ifdef ENABLE_CODE_PROFILER
      std::cout
        << " time_Match:              " << time_Match << std::endl
        << " time_FilterMatches:      " << time_FilterMatches << std::endl
        << " time_UpdateParams_PostMatch: " << time_UpdateParams_PostMatch << std::endl
        << " time_EvalErrorFunc:      " << time_EvalErrorFunc << std::endl
        << " time_Callbacks:          " << time_Callbacks << std::endl
        << " time_Extras:             " << time_Extras << std::endl;
      codeProfiler.Reset();
      codeProfiler.Start();
#endif
    }

#ifdef ENABLE_CODE_TRACE
    std::cout << "RegisterMatches()" << std::endl;
#endif

    // compute registration
    Freg = pAlgorithm->ICP_RegisterMatches();
    Freg0 = Freg1;
    Freg1 = Freg2;
    Freg2 = Freg;

#ifdef ENABLE_CODE_TRACE
    std::cout << "F:" << std::endl << Freg << std::endl;
#endif

    // dF = xfm from Freg1 to Freg2
    //  first go back along Freg1 then go forward along Freg2
    dF = Freg2 * Freg1.Inverse();

#ifdef ENABLE_CODE_PROFILER
    time_Register = codeProfiler.GetElapsedTime();
    codeProfiler.Reset();
    codeProfiler.Start();
#endif

#ifdef ENABLE_CODE_TRACE
    std::cout << "UpdateParameters_PostRegister()" << std::endl;
#endif

    // update algorithm's post-registration step parameters
    pAlgorithm->ICP_UpdateParameters_PostRegister(Freg);

#ifdef ENABLE_CODE_PROFILER
    time_UpdateParams_PostRegister = codeProfiler.GetElapsedTime();
    codeProfiler.Reset();
    codeProfiler.Start();
#endif

#ifdef ENABLE_CODE_TRACE
    std::cout << "EvaluateErrorFunction()" << std::endl;
#endif

    // On a rare iteration (typically one that involves adding a new outlier) the cost function
    //  has been seen to slightly increase, despite using the outlier threshold error correction.
    // Note: One idea is this may occur when the avg error over the entire set increases,
    //       thereby increasing the outlier threshold (and hence the outlier error contribution)
    //       for all outliers in the set. This has not been confirmed, and the true source is
    //       yet unknown.
    // Note: On the rare iterations when the rms error increases, the weighted point match
    //       registration function still reduces the rms error for the points being matched.
    //       So the increased error has something to do with how the outliers are handled.
    //if (E2 > E1)
    //{ std::cout << "  ---> Runtime Warning: cost function increased!" << std::endl; }

    // compute error function value
    E = pAlgorithm->ICP_EvaluateErrorFunction();
    E0 = E1;
    E1 = E2;
    E2 = E;
    tolE = fabs((E2 - E1) / E1);
    if (E <= Ebest)
    {
      Ebest = E;
      Fbest = Freg;
      iterBest = iter;
    }

#ifdef ENABLE_CODE_PROFILER
    time_EvalErrorFunc = codeProfiler.GetElapsedTime();
    codeProfiler.Reset();
    codeProfiler.Start();
#endif

#ifdef ENABLE_CODE_TRACE
    std::cout << "Callbacks" << std::endl;
#endif

    //-- Callbacks --//
    iterData.iter = iter;
    iterData.E = E;
    iterData.tolE = tolE;
    iterData.Freg.Assign(Freg);
    iterData.dF.Assign(dF);
    iterData.time = iterTimer.GetElapsedTime();
    iterData.nOutliers = nOutliers;
    //iterData.isAccelStep = JustDidAccelStep;
    std::vector<Callback>::iterator cbIter;
    for (cbIter = this->iterationCallbacks.begin(); cbIter != this->iterationCallbacks.end(); cbIter++)
    {
      cbIter->cbFunc(iterData, cbIter->userData);
    }
    iterTimer.Reset();
    iterTimer.Start();

#ifdef ENABLE_CODE_PROFILER
    time_Callbacks = codeProfiler.GetElapsedTime();
    codeProfiler.Reset();
    codeProfiler.Start();
#endif

#ifdef ENABLE_CODE_PROFILER
    std::cout
      << " time_Register:                  " << time_Register << std::endl
      << " time_UpdateParams_PostRegister: " << time_UpdateParams_PostRegister << std::endl
      << " time_Match:                     " << time_Match << std::endl
      << " time_UpdateParams_PostMatch:    " << time_UpdateParams_PostMatch << std::endl
      << " time_FilterMatches:             " << time_FilterMatches << std::endl
      << " time_EvalErrorFunc:             " << time_EvalErrorFunc << std::endl
      << " time_Callbacks:                 " << time_Callbacks << std::endl
      << " time_Extras:                    " << time_Extras << std::endl;
    codeProfiler.Reset();
    codeProfiler.Start();
#endif


#ifdef ENABLE_CODE_TRACE
    std::cout << "Termination Test" << std::endl;
#endif

    //-- Termination Test --//

    dR.From(dF.Rotation());   // convert rotation to Rodrigues form
    dAng = dR.Norm();
    dPos = dF.Translation().Norm();
    dAng01 = dAng12;
    dAng12 = dAng;
    dPos01 = dPos12;
    dPos12 = dPos;

    // Algorithm specific termination
    //  also enables algorithm to update the registration
    //  to a different iteration if desired
    if (pAlgorithm->ICP_Terminate(Freg))
    {
      totalTimer.Stop();
      termMsg << std::endl << "Termination Condition:  Termination Requested by Algorithm" << std::endl;
      break;  // exit iteration loop
    }

    // Consider termination
    if (dAng < opt.dAngThresh && dPos < opt.dPosThresh)
    //if (JustDidAccelStep == false
    //  && (dAng < opt.dAngThresh && dPos < opt.dPosThresh)
    //  )
    {
      // Termination Test
      //  Note: max iterations is enforced by for loop
      if ((dAng < opt.dAngTerm && dPos < opt.dPosTerm)
        || E < opt.minE
        || tolE < opt.tolE)
      {
        // termination condition must be satisfied for min number of consecutive iterations
        terminateIter++;
        if (terminateIter >= opt.termHoldIter)
        {
          // prepare termination message
          totalTimer.Stop();
          termMsg << std::endl << "Termination Condition satisfied for " << opt.termHoldIter << " iterations: " << std::endl;
          if (E < opt.minE) termMsg << "reached minE (" << opt.minE << ")" << std::endl;
          else if (tolE < opt.tolE) termMsg << "reached min dE/E (" << opt.tolE << ")" << std::endl;
          else termMsg << "reached min dAngle & dTrans (" << opt.dAngTerm * 180 / cmnPI << "/" << opt.dPosTerm << ")" << std::endl;
          break;  // exit iteration loop
        }
      }
      else
      {
        terminateIter = 0;
      }
    }
    else
    {
      terminateIter = 0;
    }

    if (iter == opt.maxIter)
    {
      // prepare termination message
      totalTimer.Stop();
      termMsg << std::endl << "Termination Condition: reached max iteration (" << opt.maxIter << ")" << std::endl;
    }

#ifdef ENABLE_CODE_PROFILER
    time_Extras = codeProfiler.GetElapsedTime();
    codeProfiler.Reset();
    codeProfiler.Start();
#endif

  }
  iterTimer.Stop();

  // complete termination message
  termMsg << " E: " << E << std::endl;
  termMsg << " dE/E: " << tolE << std::endl;
  termMsg << " dAng: " << dAng12 * 180 / cmnPI << " " << dAng01*180.0 / cmnPI << " (deg)" << std::endl;
  termMsg << " dPos: " << dPos12 << " " << dPos01 << std::endl;
  termMsg << " iter: " << iter << std::endl;
  termMsg << " runtime: " << totalTimer.GetElapsedTime() << std::endl;
  termMsg << std::endl << Freg << std::endl;
  if (iterData.iter != iterBest)
  {
    termMsg << std::endl << "WARNING: best iteration (" << iterBest << ") is not the final iteration" << std::endl;
  }
  //std::cout << termMsg.str().c_str();

  rt.termMsg = termMsg.str();
  rt.Freg = Freg;    
  rt.runTime = totalTimer.GetElapsedTime();
  rt.numIter = iter;  
  rt.nOutliers = nOutliers;

  // compute final match distance
  pAlgorithm->ComputeMatchStatistics(rt.MatchPosErrAvg, rt.MatchPosErrSD);

  return rt;
}


void cisstICP::AddIterationCallback(Callback &callback)
{
  this->iterationCallbacks.push_back(callback);
}

void cisstICP::AddIterationCallbacks(std::vector<Callback> &callbacks)
{
  if (callbacks.size() > 0)
  {
    Callback *callbackArray = callbacks.data();
    iterationCallbacks.insert(iterationCallbacks.end(), callbackArray, callbackArray + callbacks.size());
  }
}

void cisstICP::ClearIterationCallbacks()
{
  this->iterationCallbacks.clear();
}


//// Code for Accelerated ICP
////  NOTE: there are bugs in this code in that the accelerated ICP step seems to
////        not actually be applied once computed
//template <class T> inline T __cisststMIN(const T& a, const T& b) { return a < b ? a : b; };
//// Accelerated ICP
//int DecreaseCount = 0;
//int IncreaseCount = 0;
//int TinyChangeCount=0;
//double res0Save;
//double resPreAccel;
//int DecreaseCountPreAccel;
//vctFrm3 FregPreAccel, Freg0Save;
//double e01,e12,e2,s01,s12,s2;
//double res0, res1, res2;
//double accelThreshConst;
////// number of sequential increases / decreases in residual error
////DecreaseCount = (res2 < res1) ? (1 + DecreaseCount) : 0;
////IncreaseCount = (res2 > res1) ? (1 + IncreaseCount) : 0;
////
//// check for tiny change in transformation parameters
////double dRmin = 0.00001;
////double dPmin = 0.00001;
////if ((e12 < dPmin) && (s12 < dRmin)) ++TinyChangeCount;
////else TinyChangeCount = 0;
////
////// Consider termination
//////  *** why if decrease count > 1 ? 
////if ( JustDidAccelStep == 0 &&   // don't terminate immediately after accelation step
////     // don't terminate when E is oscillating, unless dF is repeatedly very small
////     (TinyChangeCount > 2 || DecreaseCount > 1 || res2 == res1 || IncreaseCount > 5) &&
////     iter > 2)                  // run at least 2 iterations
////{
////  if ( s01 < dAngThresh && s12 < dAngThresh && e01 < dPosThresh && e12 < dPosThresh )
////  {
////    // vctFrm3 is not changing very much
////    double ratio1 = res1 / res0;
////    double ratio2 = res2 / res1;
//
////    if ( TinyChangeCount > 2 ||
////         ((ratio1 > resRatioThresh && ratio1 < (1.0 / resRatioThresh)) &&
////         (ratio2 > resRatioThresh && ratio2 < (1.0 / resRatioThresh))))
////    { // residual ratio has also stopped changing
////      // *** matches not filtered in calculation of Fbest?
////      Freg = Fbest;
////      // SDB: commented these out since now calculating resRMS after
////      //      rather than before each solved registration
////      //ICP_Match();
////      //resRMS = UpdateSampleXfmPositions();
////      break;  
////    }
////  }
////}
//
////--- Accelerated ICP ---//
////
//// NOTE: this adds about 4ms to the loop runtime
////
//
//// TODO: F12 should be Freg2 * Freg1.Inverse()???
////F01 = F12;
//F01 = Freg0.Inverse() * Freg1;
//F12 = Freg1.Inverse() * Freg2;
//F02 = Freg0.Inverse() * Freg2;
//
//vctQuatRot3 q01(F01.Rotation()); vct3 q01v(q01.X(), q01.Y(), q01.Z());
//vctQuatRot3 q12(F12.Rotation()); vct3 q12v(q12.X(), q12.Y(), q12.Z());
//vctQuatRot3 q02(F02.Rotation()); vct3 q02v(q02.X(), q02.Y(), q02.Z());
//
//e01 = F01.Translation().Norm();
//e12 = F12.Translation().Norm();
//e2 = vctDotProduct(F01.Translation(), F12.Translation());  // *** forgot sqrt?
//s01 = q01v.Norm();
//s12 = q12v.Norm();
//s2 = q01v * q12v;   // *** forgot sqrt?
//
//if (JustDidAccelStep)
//{ // did it help???
//  if (res2 > res1)
//  { // oops, better undo (error got worse)
//    Freg1 = Freg0; res1 = res0;
//    Freg0 = Freg0Save; res0 = res0Save;
//    Freg = Freg2 = FregPreAccel; resRMS = res2 = resPreAccel;
//    DecreaseCount = DecreaseCountPreAccel;
//    JustDidAccelStep = false;
//
//    std::cout << "     undoing accel step ... " << std::endl;
//
//    continue; // skip accel step
//  }
//}
//
//// Consider Accelerating
//// *** Linear Acceleration not fully implemented?
////   i.e. parabolic acceleration is considered, but not linear acceleration
//// *** accelThreshConst not initialized
//JustDidAccelStep = false;
//if (accelThreshConst > 0.0 &&
//  s2 > 0.9 * s01 * s12 &&
//  e2 > 0.9 * e01 * e12 &&
//  iter > 2)   // need dF values from past two iterations before accelerating
//{
//  // Linear Acceleration
//  //  find linear fit for recent xfm parameters ("thetas" in dF.Translation() and q.v)
//  //  i.e., find a value lam1 that best approximates:
//  //    theta1[j] = lam1*theta2[j]
//  double lam1 = 2 * (F01.Translation() * F02.Translation() + q01v * q02v) /
//    (F02.Translation() * F02.Translation() + q02v * q02v);
//  // sanity check
//  if (lam1 <= 0 || lam1 >= 1)
//  { // not sensible
//    continue;
//  }
//
//  // Parabolic Acceleration
//  //  find parabolic fit for recent xfm parameters
//  //  i.e. fit the model res = a+b*lam+c*lam**2
//  //      a = res0
//  //      b+c = res2-a
//  //      b*lam1+c*lam1**2 = res1-a
//  //  i.e.
//  //      b*lam1+(res2-a-b)*lam1**2 = res1-a
//  //  i.e.
//  //      b*(lam1-lam1**2)=res1-a-(res2-a)*lam1**2
//  double a = res0;
//  double b = ((res1 - res0) - (res2 - res0) * lam1 * lam1) / (lam1 - lam1 * lam1);
//  double c = res2 - res0 - b;
//  // find a value of lambda that minimizes |a+b*lambda+c*lambda**2|
//  double lambdaMin;
//  double d = b * b - 4 * a * c;
//  if (d <= 0)
//  {
//    lambdaMin = -b / (2 * c);
//  }
//  else
//  {
//    double dd = sqrt(d);
//    double root1 = (-b + dd) / (2 * c);
//    double root2 = (-b - dd) / (2 * c);
//    lambdaMin = (root1 < 0) ? root2 : ((root2 < 0) ? root1 : __cisststMIN(root1, root2));
//  }
//  if (lambdaMin <= 0) continue;
//
//  // Now cut this off at some maximum
//  double lambda = __cisststMIN(lambdaMin, accelThreshConst);
//  double s02squared = q02v * q02v;
//  if (lambda * lambda * s02squared > 0.01)
//  {
//    lambda = sqrt(0.01 / s02squared);
//  }
//
//  // Now guess an accelerated dF
//  vct3 dqAccel;
//  vct3 dPaccel;
//  for (int j = 0; j<3; j++) {
//    // compute accelerated angles
//    dqAccel(j) = q02v(j)*lambda;
//    dPaccel(j) = F02.Translation()[j] * lambda;
//  }
//
//  FregPreAccel = Freg;
//  resPreAccel = res2;
//  DecreaseCountPreAccel = DecreaseCount;
//  vctFrm3 dFaccel(vctRot3(vctRodRot3(dqAccel[0], dqAccel[1], dqAccel[2])), dPaccel);
//  Freg = Freg2 = Freg0 * dFaccel;
//  JustDidAccelStep = true;
//
//  std::cout << std::endl << "=== Doing Accel Step ===" << std::endl << std::endl;
//}
