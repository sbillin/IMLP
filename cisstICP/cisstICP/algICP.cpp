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

#include "algICP.h"


#ifdef ValidatePDTreeSearch
std::ofstream validFS("../ICP_TestData/LastRun/debugPDTreeSearchValidation.txt");
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

namespace {

  // declerations
  void ICPCallback_PrintIteration(cisstICP::CallbackArg &arg, void *userData);

  // default callback
  void ICPCallback_PrintIteration(cisstICP::CallbackArg &arg, void *userData)
  {
    vctRodRot3 dR(arg.dF.Rotation());
    std::stringstream ss;
    algICP *pThis = (algICP*)userData;

    std::string fstring("i=%u E=%.1f tolE=%.3f  t=%.3f NNodes=%u/%u/%u NOut=%u");
    ss << cmnPrintf(fstring.c_str())
      // \tt=%.3f\tNNodes=%u/%u/%u
      // (RMS/Res)=%.4f/%.4f
      // (dAng/dPos)= %.2f/%.2f )
      << arg.iter
      << arg.E
      << arg.tolE
      << arg.time
      << pThis->maxNodesSearched << pThis->avgNodesSearched << pThis->minNodesSearched
      << arg.nOutliers
      //<< dR.Norm()*180/cmnPI << arg.dF.Translation().Norm()
      ;
#ifdef ValidatePDTreeSearch
      fstring = "  vld=%.2f";
      ss << cmnPrintf(fstring.c_str())
        << validPercent;
#endif

    std::cout << ss.str() << std::endl;
  }

} // namespac anonymous

std::vector<cisstICP::Callback> algICP::ICP_GetIterationCallbacks()
{
  std::vector<cisstICP::Callback> callbacks;
  cisstICP::Callback defaultICPCallback(ICPCallback_PrintIteration, this);
  callbacks.push_back(defaultICPCallback);
  return callbacks;
}

// constructor
algICP::algICP(PDTreeBase *pTree, const vctDynamicVector<vct3> &samplePts)
  : pTree(pTree)
{
  SetSamples(samplePts);
}

void algICP::SetSamples(const vctDynamicVector<vct3> &argSamplePts)
{
  // copy sample points
  samplePts = argSamplePts;

  // initialize variables dependent on sample size
  nSamples = samplePts.size();

  samplePtsXfmd.SetSize(nSamples);
  matchPts.SetSize(nSamples);
  matchDatums.SetSize(nSamples);
  matchErrors.SetSize(nSamples);

  //residuals_PostMatch.SetSize(nSamples);
  //sqrDist_PostMatch.SetSize(nSamples);
  //dist_PostMatch.SetSize(nSamples);

  //residuals_PostRegister.SetSize(nSamples);
  //sqrDist_PostRegister.SetSize(nSamples);
  //dist_PostRegister.SetSize(nSamples);
}

void algICP::ICP_InitializeParameters(vctFrm3 &FGuess)
{
  // set starting sample positions
  UpdateSampleXfmPositions(FGuess);

  // initialize matches with accelerated approximate search
  for (unsigned int i = 0; i < nSamples; i++)
  {
    matchDatums[i] = pTree->FastInitializeProximalDatum(
      samplePtsXfmd[i], matchPts[i]);
  }

  //// initialize matches to any model point
  ////  i.e. we don't know the closest match => set it to anything valid
  //for (unsigned int s = 0; s < nSamples; s++)
  //{
  //  matchDatums.Element(s) = 0;
  //  matchPts.Element(s) = pTree->DatumSortPoint(0);
  //}

  nOutliers = 0;
  Freg = FGuess;
}

void algICP::ICP_UpdateParameters_PostMatch()
{
  //ComputeErrors_PostMatch();
}

void algICP::ICP_UpdateParameters_PostRegister(vctFrm3 &Freg)
{
  UpdateSampleXfmPositions(Freg);

  //ComputeErrors_PostRegister();
}

void algICP::ICP_ComputeMatches()
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
    matchDatums.Element(s) = pTree->FindClosestDatum(
      samplePtsXfmd.Element(s),
      matchPts.Element(s),
      matchDatums.Element(s),
      matchErrors.Element(s),
      nodesSearched);

    avgNodesSearched += nodesSearched;
    minNodesSearched = (nodesSearched < minNodesSearched) ? nodesSearched : minNodesSearched;
    maxNodesSearched = (nodesSearched > maxNodesSearched) ? nodesSearched : maxNodesSearched;

#ifdef ValidatePDTreeSearch
#ifdef ValidateByEuclideanDist
    validDatum = pTree->ValidateClosestDatum_ByEuclideanDist(samplePtsXfmd.Element(s), validPoint);
#else
    validDatum = pTree->ValidateClosestDatum(samplePtsXfmd.Element(s), validPoint);
#endif
    if (validDatum != matchDatums.Element(s))
    {
      // It is possible to have different datums for same point if the match
      //  lies on a datum edge; if this is the case, then the search did not
      //  actually fail
      searchError = (validPoint - matchPts.Element(s)).NormSquare();
      // Note: cannot compare two double values for exact equality due to
      //       inexact representation of decimal values in binary arithmetic
      //       => use an epsilon value for comparison
      if (searchError > doubleEps)
      {
        double ResidualDistance = (matchPts.Element(s) - samplePtsXfmd.Element(s)).Norm();
        validDist = (validPoint - samplePtsXfmd.Element(s)).Norm();
        numInvalidDatums++;
        vct3 tmp1;

        searchError = pTree->pAlgorithm->FindClosestPointOnDatum(
          samplePtsXfmd.Element(s), tmp1, matchDatums.Element(s));
        validError = pTree->pAlgorithm->FindClosestPointOnDatum(
          samplePtsXfmd.Element(s), tmp1, validDatum);
        validFS << "Match Errors = " << searchError << "/" << validError
          << "\t\tdPos = " << ResidualDistance << "/" << validDist << std::endl;
        validFS << " XfmSamplePoint = [" << samplePtsXfmd.Element(s) << "]" << std::endl;
        validFS << " SearchPoint = [" << matchPts.Element(s) << "]" << std::endl;
        validFS << " SearchDatum = " << matchDatums.Element(s) << std::endl;
        validFS << " ValidPoint =  [" << validPoint << "]" << std::endl;
        validFS << " ValidDatum = " << validDatum << std::endl;
        validFS << " SampleIndex = " << s << std::endl;

        //algICP_IMLP_PointCloud *alg;
        //alg = static_cast<algICP_IMLP_PointCloud*>(algorithm);
        //algICP_MahalDist_PointCloud *alg;
        //alg = static_cast<algICP_MahalDist_PointCloud*>(algorithm);
        //PDTree_PointCloud *tree;
        //tree = static_cast<PDTree_PointCloud*>(Trees[0]);

        //validFS << " RMxiRt = [" << std::endl << alg->R_Mxi_Rt[s] << "]" << std::endl;
        //validFS << " RMxRt_sigma2 = [" << std::endl << alg->sample_RMxRt_sigma2 << "]" << std::endl;
        //validFS << " Mxi = [" << std::endl << alg->Mxi[s] << "]" << std::endl;
        //validFS << " Myi_search = [" << std::endl << tree->pointCov(matchDatums.Element(s)) << "]" << std::endl;
        //validFS << " Myi_valid = [" << std::endl << tree->pointCov(validDatum) << "]" << std::endl;
        //validFS << "Fact = [" << std::endl << pICP->Fact << "]" << std::endl;

        PDTreeNode *termNode = 0;
        pTree->FindTerminalNode(validDatum, &termNode);
        if (!termNode)
        {
          std::cout << "ERROR: did not find terminal node for datum: " << validDatum << std::endl;
          assert(0);
        }
        validFS << " Valid Terminal Node:" << std::endl;
        validFS << "   MinCorner: " << termNode->Bounds.MinCorner << std::endl;
        validFS << "   MaxCorner: " << termNode->Bounds.MaxCorner << std::endl;
        validFS << "   NData: " << termNode->NData << std::endl;
        validFS << "Fnode = [" << std::endl << termNode->F << "]" << std::endl;

        // crash program
        std::cout << "Invalid Node Search; Terminating Program" << std::endl;
        termNode = NULL;
        termNode->NData = 0;
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

    // inform algorithm that match completed
    SamplePostMatch(s);
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
      fsCP << matchPts.at(i) << std::endl;
      fsSP << samplePtsXfmd.at(i) << std::endl;
    }
    fsCP.close();
    fsSP.close();
    saveMatchesIter++;
  }
#endif

}

unsigned int algICP::ICP_FilterMatches()
{
  return 0;
}

void algICP::UpdateSampleXfmPositions(const vctFrm3 &F)
{
  for (unsigned int s = 0; s < nSamples; s++)
  {
    samplePtsXfmd.Element(s) = F * samplePts.Element(s);
  }
}

//// this is a convenience routine that is not called by the base algorithm
//void algICP::ComputeErrors_PostMatch()
//{
//  //matchErrorAvg_PostMatch = 0.0;
//  matchDistAvg_PostMatch = 0.0;
//  sumSqrDist_PostMatch = 0.0;
//
//  for (unsigned int s = 0; s < nSamples; s++)
//  {
//    residuals_PostMatch.Element(s) = samplePtsXfmd.Element(s) - matchPts.Element(s);
//    sqrDist_PostMatch.Element(s) = residuals_PostMatch.Element(s).NormSquare();
//    dist_PostMatch.Element(s) = sqrt(sqrDist_PostMatch.Element(s));
//
//    sumSqrDist_PostMatch += sqrDist_PostMatch.Element(s);
//    matchDistAvg_PostMatch += dist_PostMatch.Element(s);
//    //matchErrorAvg_PostMatch += matchErrors.Element(s);
//  }
//
//  //matchErrorAvg_PostMatch /= nSamples;
//  matchDistAvg_PostMatch /= nSamples;
//}

//// this is a convenience routine that is not called by the base algorithm
//void algICP::ComputeErrors_PostRegister()
//{
//  //matchErrorAvg_PostRegister = 0.0;
//  matchDistAvg_PostRegister = 0.0;
//  sumSqrDist_PostRegister = 0.0;
//
//  for (unsigned int s = 0; s < nSamples; s++)
//  {
//    residuals_PostRegister.Element(s) = samplePtsXfmd.Element(s) - matchPts.Element(s);
//    sqrDist_PostRegister.Element(s) = residuals_PostRegister.Element(s).NormSquare();
//    dist_PostRegister.Element(s) = sqrt(sqrDist_PostRegister.Element(s));
//
//    sumSqrDist_PostRegister += sqrDist_PostRegister.Element(s);
//    matchDistAvg_PostRegister += dist_PostRegister.Element(s);
//    //matchErrorAvg_PostRegister += matchErrors.Element(s);
//  }
//
//  //matchErrorAvg_PostRegister /= nSamples;
//  matchDistAvg_PostRegister /= nSamples;
//}

void algICP::ComputeMatchStatistics(double &Avg, double &StdDev)
{
  double sumSqrMatchDist = 0.0;
  double sumMatchDist = 0.0;
  double sqrMatchDist;

  // NOTE: if using a method with outlier rejection, it may be desirable to
  //       compute statistics on only the inliers
  for (unsigned int i = 0; i < nSamples; i++)
  {
    sqrMatchDist = (matchPts[i] - Freg * samplePts[i]).NormSquare();

    sumSqrMatchDist += sqrMatchDist;
    sumMatchDist += sqrt(sqrMatchDist);
  }

  Avg = sumMatchDist / nSamples;
  StdDev = (sumSqrMatchDist / nSamples) + Avg*Avg;
}