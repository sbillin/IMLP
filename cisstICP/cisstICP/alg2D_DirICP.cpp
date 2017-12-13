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

#include "alg2D_DirICP.h"

#include <limits>


#ifdef ValidatePDTreeSearch
//#define ValidateByEuclideanDist
std::ofstream validFS("../ICP_TestData/LastRun/debugDirPDTreeSearchValidation.txt");
vct2 validPoint;
vct2 validNorm;
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


// constructor
alg2D_DirICP::alg2D_DirICP(
  DirPDTree2DBase *pDirTree,
  vctDynamicVector<vct2> &samplePts,
  vctDynamicVector<vct2> &sampleNorms)
  : pDirTree(pDirTree)
{
  SetSamples(samplePts, sampleNorms);
}

//std::vector<cisstICP::Callback> alg2D_DirICP::ICP_GetIterationCallbacks()
//{
//  CISST_THROW("Not Implemented");
//}

void alg2D_DirICP::SetSamples(
  vctDynamicVector<vct2> &argSamplePts,
  vctDynamicVector<vct2> &argSampleNorms)
{
  if (argSamplePts.size() != argSampleNorms.size())
  {
    CISST_THROW("ERROR: sample points and normals are different sizes");
  }

  // copy sample points
  samplePts = argSamplePts;
  sampleNorms = argSampleNorms;

  nSamples = samplePts.size();

  // initialize variables dependent on sample size

  samplePtsXfmd.SetSize(nSamples);
  matchPts.SetSize(nSamples);
  matchDatums.SetSize(nSamples);
  matchErrors.SetSize(nSamples);

  sampleNormsXfmd.SetSize(nSamples);
  matchNorms.SetSize(nSamples);

  //residuals_PostMatch.SetSize(nSamples);
  //sqrDist_PostMatch.SetSize(nSamples);
  //dist_PostMatch.SetSize(nSamples);

  //residuals_PostRegister.SetSize(nSamples);
  //sqrDist_PostRegister.SetSize(nSamples);
  //dist_PostRegister.SetSize(nSamples);

  //normProducts_PostMatch.SetSize(nSamples);
  //normProducts_PostRegister.SetSize(nSamples);
}

void alg2D_DirICP::ICP_InitializeParameters(vctFrm2 &FGuess)
{
  // set starting sample positions
  UpdateSampleXfmPositions(FGuess);

  // initialize matches with accelerated approximate search
  for (unsigned int i = 0; i < nSamples; i++)
  {
    matchDatums[i] = pDirTree->FastInitializeProximalDatum(
      samplePtsXfmd[i], sampleNormsXfmd[i],
      matchPts[i], matchNorms[i]);
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
}

void  alg2D_DirICP::ICP_UpdateParameters_PostMatch()
{
  //ComputeErrors_PostMatch();;
}

void  alg2D_DirICP::ICP_UpdateParameters_PostRegister(vctFrm2 &Freg)
{
  UpdateSampleXfmPositions(Freg);

  //ComputeErrors_PostRegister();
}

void alg2D_DirICP::ICP_ComputeMatches()
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
      validPoint, validNorm );
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
        validAng  = acos( vctDotProduct(validNorm,sampleNormsXfmd.Element(s)) );            
        vct2 tmp1,tmp2;           
        //searchError = algorithm->FindClosestPointOnDatum(
        //                    samplePtsXfmd.Element(s), sampleNormsXfmd.Element(s),
        //                    tmp1, tmp2, matchDatums.Element(s));
        validError = pDirTree->pAlgorithm->FindClosestPointOnDatum(
          samplePtsXfmd.Element(s), sampleNormsXfmd.Element(s),
          tmp1, tmp2, validDatum);
        validFS << "Correspondence Search Data:  (foundMatch / validMatch)" << std::endl
          << " MatchError = " << matchErrors.Element(s) << "/" << validError << std::endl
          << " dPos = " << dist_PostMatch.Element(s) << "/" << validDist << std::endl
          << " dAng = " << dTheta*180.0/cmnPI << "/" << validAng*180.0/cmnPI << std::endl
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
        pDirTree->FindTerminalNode( validDatum, &termNode );
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
      fsCP << matchPts.at(i) << " " << matchNorms.at(i) << std::endl;
      fsSP << samplePtsXfmd.at(i) << " " << sampleNormsXfmd.at(i) << std::endl;
    }
    fsCP.close();
    fsSP.close();
    saveMatchesIter++;
  }
#endif

}

unsigned int alg2D_DirICP::ICP_FilterMatches()
{
  return 0;
}

void alg2D_DirICP::UpdateSampleXfmPositions(const vctFrm2 &F)
{
  vctRot2 R = F.Rotation();
  for (unsigned int s = 0; s < nSamples; s++)
  {
    samplePtsXfmd.Element(s) = F * samplePts.Element(s);
    sampleNormsXfmd.Element(s) = R * sampleNorms.Element(s);
  }
}

//void alg2D_DirICP::ComputeErrors_PostMatch()
//{
//  //matchErrorAvg_PostMatch = 0.0;
//  matchDistAvg_PostMatch = 0.0;
//  sumSqrDist_PostMatch = 0.0;
//  sumNormProducts_PostMatch = 0.0;
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
//
//    normProducts_PostMatch.Element(s) = vctDotProduct(sampleNormsXfmd.Element(s), matchNorms.Element(s));
//    sumNormProducts_PostMatch += normProducts_PostMatch.Element(s);
//  }
//
//  //matchErrorAvg_PostMatch /= nSamples;
//  matchDistAvg_PostMatch /= nSamples;
//}
//
//void alg2D_DirICP::ComputeErrors_PostRegister()
//{
//  //matchErrorAvg_PostRegister = 0.0;
//  matchDistAvg_PostRegister = 0.0;
//  sumSqrDist_PostRegister = 0.0;
//  sumNormProducts_PostRegister = 0.0;
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
//
//    normProducts_PostRegister.Element(s) += vctDotProduct(sampleNormsXfmd.Element(s), matchNorms.Element(s));
//    sumNormProducts_PostRegister += normProducts_PostRegister.Element(s);
//  }
//
//  //matchErrorAvg_PostRegister /= nSamples;
//  matchDistAvg_PostRegister /= nSamples;
//}


// *** May be required post-match for filtering matches
void alg2D_DirICP::ComputeCircErrorStatistics(double sumNormProducts, double &R, double &circSD)
{
  // Angular Error of Normal Vectors:
  //   Circular Standard Deviation: (radians)
  //     R = Sum_i(dot(Ny,Nx))/N
  //     circSD = stdDev(theta) = sqrt(-2*ln(R))
  //     where theta is the angle between matched vectors; i.e. theta = acos(dot(Ny,Nx))
  R = sumNormProducts / nSamples;

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
