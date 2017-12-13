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
#include "algDirICP.h"


#ifdef ValidatePDTreeSearch
//#define ValidateByEuclideanDist
//std::ofstream validFS("../ICP_TestData/LastRun/debugDirPDTreeSearchValidation.txt");
std::ofstream validFS("debugDirPDTreeSearchValidation.txt");
vct3 validPoint;
vct3 validNorm;
int  validDatum;
double validDist;
double validAng;
double validError;
double searchError;
unsigned int numValidDatums;
unsigned int numInvalidDatums;
double validPercent;
double doubleEps = 1e-16;
int validIter = 0;
#endif

#ifdef SaveMatchesToFile
std::string saveMatchesDir("../ICP_TestData/LastRun");
int saveMatchesIter = 0;
#endif

// declerations
void DirICPCallback_PrintIteration(cisstICP::CallbackArg &arg, void *userData);

// default callback
void DirICPCallback_PrintIteration(cisstICP::CallbackArg &arg, void *userData)
{
  std::stringstream ss;
  algDirICP *pThis = (algDirICP*)userData;
  vctAxAnRot3 dRAxAn(arg.dF.Rotation());

  //// for all iterations we display post-registration values, except for
  ////  iteration 0 which occurs before the first registration step
  //double circSD;
  //double sumSqrDist;
  //if (arg.iter == 0)
  //{
  //  circSD = pThis->circSD;
  //  sumSqrDist = pThis->sumSqrDist_PostRegister;
  //}
  //else
  //{
  //  circSD = pThis->circSD_PostReg;
  //  sumSqrDist = pThis->sumSqrDist_PostReg;
  //}

  //std::string fstring("i=%u E=%.1f tolE=%.3f  t=%.3f NNodes=%u/%u/%u NOut=%u");
  //std::string fstring("i=%u E=%.1f tolE=%.3f  nW/pW=%.1f/%.2f (cSD/pSD)=%.1f/%.1f NNodes=%u/%u/%u");
  std::string fstring("i=%u E=%.1f nW/pW=%.1f/%.2f  (dR/dT)= %.3f/%.3f  NN=%u/%u/%u NOut=%u");
#ifdef ValidatePDTreeSearch
  fstring.append("  vld=%.2f");
#endif

  ss << cmnPrintf(fstring.c_str())
    // (RMS/Res)=%.4f/%.4f
    // t=%.3f
    // dTheta=%.2f/%.2f/%.2f
    // (dAng/dPos)= %.2f/%.2f
    // (cSD/dSD)=%.1f/%.1f 
    // (dA/dP)=%.1f/%.1f
    // NNodes=%u/%u/%u
    // nOut=%u
    // nOut=%u/%u    
    << arg.iter
    << arg.E
    //<< arg.tolE
    << pThis->errFuncNormWeight << pThis->errFuncPosWeight
    << dRAxAn.Angle()*180.0 / cmnPI << arg.dF.Translation().Norm()
    //<< circSD * 180.0/cmnPI << sqrt(sumSqrDist / pThis->nSamples)
    << pThis->maxNodesSearched << pThis->avgNodesSearched << pThis->minNodesSearched
    << arg.nOutliers    
    //<< argp->dAngAvg*180.0/cmnPI << argp->dPosAvg
    //<< argp->normCircSD*180.0/cmnPI << argp->posCircSD*180.0/cmnPI
    //<< argp->dThetaMin*180.0/cmnPI << argp->dThetaMax*180.0/cmnPI << argp->dThetaAvg*180.0/cmnPI
    //<< arg.time
    //<< argp->nOutliers
    //<< argp->nOutliersPos << argp->nOutliersNorm
#ifdef ValidatePDTreeSearch
    << validPercent
#endif
    ;

  std::cout << ss.str() << std::endl;
}

// constructor
algDirICP::algDirICP(
  DirPDTreeBase *pDirTree,
  const vctDynamicVector<vct3> &samplePts,
  const vctDynamicVector<vct3> &sampleNorms)
  : algICP(NULL, samplePts),
  pDirTree(pDirTree)
{
  SetSamples(samplePts, sampleNorms);
}

std::vector<cisstICP::Callback> algDirICP::ICP_GetIterationCallbacks()
{
  std::vector<cisstICP::Callback> callbacks;
  cisstICP::Callback defaultICPCallback(DirICPCallback_PrintIteration, this);
  callbacks.push_back(defaultICPCallback);
  return callbacks;
}

void algDirICP::SetSamples(
  const vctDynamicVector<vct3> &argSamplePts,
  const vctDynamicVector<vct3> &argSampleNorms)
{
  if (argSamplePts.size() != argSampleNorms.size())
  {
    std::cout << "ERROR: sample points and normals are different sizes" << std::endl;
    assert(0);
    return;
  }

  // base class
  algICP::SetSamples(argSamplePts);

  // copy normals
  sampleNorms = argSampleNorms;

  // initialize normals variables dependent on sample size
  sampleNormsXfmd.SetSize(nSamples);
  matchNorms.SetSize(nSamples);

  //normProducts_PostMatch.SetSize(nSamples);
  //normProducts_PostRegister.SetSize(nSamples);
}

void algDirICP::ICP_InitializeParameters(vctFrm3 &FGuess)
{
  // do not call base class parameter initialization since
  //  the match initialization requires a different routine currently
  // TODO: make the match initialization routine a method of the base class that does not
  //       depend on Dir vs non-Dir PD tree, since normals
  //       are not used for this anyway

  // set starting sample positions
  UpdateSampleXfmPositions(FGuess);

  // initialize matches with accelerated approximate search
  if (pDirTree)
  {
    for (unsigned int i = 0; i < nSamples; i++)
    {
      matchDatums[i] = pDirTree->FastInitializeProximalDatum(
        samplePtsXfmd[i], sampleNormsXfmd[i], matchPts[i], matchNorms[i]);
    }
  }

  //// initialize matches to any model point
  ////  i.e. we don't know the closest match => set it to anything valid
  //for (unsigned int s = 0; s < nSamples; s++)
  //{
  //  matchDatums.Element(s) = 0;
  //  matchPts.Element(s) = pDirTree->DatumSortPoint(0);
  //  matchNorms.Element(s) = pDirTree->DatumNorm(0);
  //}

  nOutliers = 0;
  Freg = FGuess;
}

void  algDirICP::ICP_UpdateParameters_PostMatch()
{
  //ComputeErrors_PostMatch();
}

void  algDirICP::ICP_UpdateParameters_PostRegister(vctFrm3 &Freg)
{
  UpdateSampleXfmPositions(Freg);

  //ComputeErrors_PostRegister();
}

void algDirICP::ICP_ComputeMatches()
{
  // Find the point on the model having lowest match error
  //  for each sample point

#ifdef ValidatePDTreeSearch  
  numInvalidDatums = 0;
  numValidDatums = 0;
#endif

  unsigned int nodesSearched = 0;
  minNodesSearched = std::numeric_limits<unsigned int>::max();
  maxNodesSearched = std::numeric_limits<unsigned int>::min();
  avgNodesSearched = 0;

  for (unsigned int s = 0; s < nSamples; s++)
  {
    // inform algorithm beginning new match
    SamplePreMatch(s);

    // Find best match for this sample
    matchDatums.Element(s) = pDirTree->FindClosestDatum(
      samplePtsXfmd.Element(s), sampleNormsXfmd.Element(s),
      matchPts.Element(s), matchNorms.Element(s),
      matchDatums.Element(s),
      matchErrors.Element(s),
      nodesSearched);

    avgNodesSearched += nodesSearched;
    minNodesSearched = (nodesSearched < minNodesSearched) ? nodesSearched : minNodesSearched;
    maxNodesSearched = (nodesSearched > maxNodesSearched) ? nodesSearched : maxNodesSearched;

#ifdef ValidatePDTreeSearch        
    validDatum = pDirTree->ValidateClosestDatum(
      samplePtsXfmd.Element(s), sampleNormsXfmd.Element(s),
      validPoint, validNorm);
    if (validDatum != matchDatums.Element(s))
    {
      // It is possible to have different datums for same point if the match
      //  lies on a datum edge; if this is the case, then the search didn't
      //  actually fail
      // Note: cannot compare two double values for exact equality due to
      //       inexact representation of decimal values in binary arithmetic
      searchError = (validPoint - matchPts.Element(s)).NormSquare();
      if (searchError > doubleEps)
      {
        numInvalidDatums++;
        //matchDist = (matchPts.Element(s)-samplePtsXfmd.Element(s)).Norm();
        //matchAng  = acos( vctDotProduct(matchNorms.Element(s),sampleNormsXfmd.Element(s)) );
        validDist = (validPoint - samplePtsXfmd.Element(s)).Norm();
        validAng = acos(vctDotProduct(validNorm, sampleNormsXfmd.Element(s)));
        vct3 tmp1, tmp2;
        //searchError = algorithm->FindClosestPointOnDatum(
        //                    samplePtsXfmd.Element(s), sampleNormsXfmd.Element(s),
        //                    tmp1, tmp2, matchDatums.Element(s));
        validError = pDirTree->pAlgorithm->FindClosestPointOnDatum(
          samplePtsXfmd.Element(s), sampleNormsXfmd.Element(s),
          tmp1, tmp2, validDatum);
        validFS << "Correspondence Search Data:  (foundMatch / validMatch)" << std::endl
          << " MatchError = " << matchErrors.Element(s) << "/" << validError << std::endl
          << " dPos = " << dist_PostMatch.Element(s) << "/" << validDist << std::endl
          //<< " dAng = " << dTheta*180.0/cmnPI << "/" << validAng*180.0/cmnPI << std::endl
          << " XfmSamplePoint = [" << samplePtsXfmd.Element(s) << "]" << std::endl
          << " XfmSampleNorm  = [" << sampleNormsXfmd.Element(s) << "]" << std::endl
          << " MatchPoint =     [" << matchPts.Element(s) << "]" << std::endl
          << " MatchNorm =      [" << matchNorms.Element(s) << "]" << std::endl
          << " MatchDatum = " << matchDatums.Element(s) << std::endl
          << " ValidPoint =     [" << validPoint << "]" << std::endl
          << " ValidNorm =      [" << validNorm << "]" << std::endl
          << " ValidDatum = " << validDatum << std::endl
          << " SampleIndex = " << s << std::endl;

        //validFS << "Fact = [" << std::endl << Fact << "]" << std::endl;

        DirPDTreeNode *termNode = 0;
        pDirTree->FindTerminalNode(validDatum, &termNode);
        if (!termNode)
        {
          std::cout << "ERROR: did not find terminal node for datum: " << validDatum << std::endl;
          assert(0);
        }
        validFS << " Valid Terminal Node:" << std::endl;
        validFS << "   MinCorner: " << termNode->Bounds.MinCorner << std::endl;
        validFS << "   MaxCorner: " << termNode->Bounds.MaxCorner << std::endl;
        validFS << "   NData: " << termNode->NData << std::endl;
        validFS << "Fnode_valid = [" << std::endl << termNode->F << "]" << std::endl;
        validFS << std::endl;
      }
      else
      {
        numValidDatums++;
      }
    }
    else
    {
      numValidDatums++;
    }
#endif
  }

  avgNodesSearched /= nSamples;

#ifdef ValidatePDTreeSearch  
  validPercent = (double)numValidDatums / (double)nSamples;
  validFS << "iter " << validIter << ":  NumMatches(valid/invalid): "
    << numValidDatums << "/" << numInvalidDatums << "  valid% = "
    << validPercent << std::endl;
  validIter++;
#endif

#ifdef SaveMatchesToFile
  {
    std::stringstream ss;
    ss << saveMatchesDir + "/matches-" << saveMatchesIter << ".txt";
    std::ofstream fsCP(ss.str().c_str());

    ss.str("");
    ss << saveMatchesDir + "/samplesXfmd-" << saveMatchesIter << ".txt";
    std::ofstream fsSP(ss.str().c_str());

    for (unsigned int i = 0; i < nSamples; i++)
    {
      fsCP << matchPts.at(i) << " " << matchNorms.at(i) << std::endl;
      fsSP << samplePtsXfmd.at(i) << " " << sampleNormsXfmd.at(i) << std::endl;
    }
    fsCP.close();
    fsSP.close();
    saveMatchesIter++;
  }
#endif

}

unsigned int algDirICP::ICP_FilterMatches()
{
  return 0;
}

void algDirICP::UpdateSampleXfmPositions(const vctFrm3 &F)
{
  vctRot3 R = F.Rotation();
  for (unsigned int s = 0; s < nSamples; s++)
  {
    samplePtsXfmd.Element(s) = F * samplePts.Element(s);
    sampleNormsXfmd.Element(s) = R * sampleNorms.Element(s);
  }
}

//void algDirICP::ComputeErrors_PostMatch()
//{
//  // base class
//  algICP::ComputeErrors_PostMatch();
//
//  //double dTheta;
//  //dThetaMin = std::numeric_limits<double>::max();
//  //dThetaMax = std::numeric_limits<double>::min();
//  //dThetaAvg = 0.0;
//
//  sumNormProducts_PostMatch = 0.0;
//
//  for (unsigned int s = 0; s < nSamples; s++)
//  {
//    normProducts_PostMatch.Element(s) = vctDotProduct(sampleNormsXfmd.Element(s), matchNorms.Element(s));
//    sumNormProducts_PostMatch += normProducts_PostMatch.Element(s);
//
//    //dTheta = acos(normProducts_PostMatch.Element(s));
//    //dTheta_PostMatch.Element(s) = dTheta;
//    //dThetaMin = dTheta < dThetaMin ? dTheta : dThetaMin;
//    //dThetaMax = dTheta > dThetaMax ? dTheta : dThetaMax;
//    //dThetaAvg += dTheta;
//  }
//
//  //dThetaAvg /= nSamples;
//}

//void algDirICP::ComputeErrors_PostRegister()
//{
//  // base class
//  algICP::ComputeErrors_PostRegister();
//
//  sumNormProducts_PostRegister = 0.0;
//
//  for (unsigned int s = 0; s < nSamples; s++)
//  {
//    normProducts_PostRegister.Element(s) = vctDotProduct(sampleNormsXfmd.Element(s), matchNorms.Element(s));
//    sumNormProducts_PostRegister += normProducts_PostRegister.Element(s);
//    //errDeg = acos(normProd) * (180.0 / cmnPI);
//  }
//}

// *** May be required post-match for filtering matches
void algDirICP::ComputeCircErrorStatistics(double sumNormProducts, double &R, double &circSD)
{
  // Angular Error of Normal Vectors:
  //   Circular Standard Deviation: (radians)
  //     R = Sum_i(dot(Ny,Nx))/N
  //     circSD = stdDev(theta) = sqrt(-2*ln(R))
  //     where theta is the angle between matched vectors; i.e. theta = acos(dot(Ny,Nx))
  R = sumNormProducts / nSamples;

  // R must be <= 1
  if (R > 1.0)
  {
    std::cout << "WARNING: R > 1, thresholding R at 1" << std::endl;
    R = 1.0;
  }

  // R must be >= 0 for a dist'n that is biased towards its mean direction
  //  (i.e. for a Fisher-like distn')
  // It is possible that all matches point in opposite directions
  //   => R could be negative from the calculation above.
  if (R <= 0.0)
  {
    R = 0.0;
    circSD = cmnPI;
    //std::cout << "WARNING: R <= 0.0; thresholding at 0.0" << std::endl;
  }
  else
  {
    circSD = sqrt(-2.0 * log(R));
    // theta can be no larger than 180 degrees
    //  => circSD must not be larger than 180 degrees either
    if (circSD > cmnPI)
    {
      circSD = cmnPI;
      //std::cout << "WARNING: circSD > 180.0 deg; thresholding at 180.0 deg" << std::endl;
    }
  }
}

void algDirICP::ComputeMatchStatistics(
  double &PosAvg, double &PosStdDev,
  double &AngAvg, double &AngStdDev)
{
  double sumSqrMatchDist = 0.0;
  double sumMatchDist = 0.0;
  double sqrMatchDist;
  double sumSqrMatchAngle = 0.0;
  double sumMatchAngle = 0.0;
  double matchAngle, sqrMatchAngle;

  // return the average match distance of the inliers
  for (unsigned int i = 0; i < nSamples; i++)
  {
    sqrMatchDist = (matchPts[i] - Freg * samplePts[i]).NormSquare();

    sumSqrMatchDist += sqrMatchDist;
    sumMatchDist += sqrt(sqrMatchDist);

    matchAngle = acos(matchNorms[i] * (Freg.Rotation() * sampleNorms[i]));

    sumMatchAngle += matchAngle;
    sumSqrMatchAngle += matchAngle * matchAngle;
  }

  PosAvg = sumMatchDist / nSamples;
  PosStdDev = (sumSqrMatchDist / nSamples) + PosAvg*PosAvg;
  AngAvg = sumMatchAngle / nSamples;
  AngStdDev = (sumSqrMatchAngle / nSamples) + AngAvg*AngAvg;
}