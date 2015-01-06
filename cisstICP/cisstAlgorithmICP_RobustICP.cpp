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
#include "cisstAlgorithmICP_RobustICP.h"
#include "cisstCovTreeNode.h"
#include "cisstICP.h"
#include "RegisterP2P.h"

#include <cmath>

namespace {

  // declerations
  void ICPCallback_PrintIteration(cisstICP::CallbackArg &arg, void *userData);

  // default callback
  void ICPCallback_PrintIteration(cisstICP::CallbackArg &arg, void *userData)
  {
    vctRodRot3 dR(arg.dF.Rotation());
    std::stringstream ss;
    cisstAlgorithmICP_RobustICP *pThis = (cisstAlgorithmICP_RobustICP*)userData;

    std::string fstring("i=%u E=%.1f tolE=%.3f  t=%.3f NNodes=%u/%u/%u NOut=%u DImax=%.1f u=%.1f sd=%.1f");
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
      << pThis->DImax
      << pThis->distAvg
      << pThis->distSD
      ;
      //<< dR.Norm()*180/cmnPI << arg.dF.Translation().Norm()
    if (arg.iter == 0)
    {
      std::string fstring(" D=%.1f eps=%f");
      ss << cmnPrintf(fstring.c_str())
        << pThis->D
        << pThis->epsilon
        ;
    }
#ifdef ValidateCovTreeSearch
    fstring = "  vld=%.2f";
    ss << cmnPrintf(fstring.c_str())
      << validPercent
      ;
#endif

    std::cout << ss.str() << std::endl;
  }

} // namespac anonymous

std::vector<cisstICP::Callback> cisstAlgorithmICP_RobustICP::ICP_GetIterationCallbacks()
{
  std::vector<cisstICP::Callback> callbacks;
  cisstICP::Callback defaultICPCallback(ICPCallback_PrintIteration, this);
  callbacks.push_back(defaultICPCallback);
  return callbacks;
}


void cisstAlgorithmICP_RobustICP::ICP_InitializeParameters(vctFrm3 &FGuess)
{
  // initialize base class params
  cisstAlgorithmICP_StdICP::ICP_InitializeParameters(FGuess);

  // size buffers
  matchDist.SetSize(nSamples);
  filterSamplePtsBuf.SetSize(nSamples);
  filterMatchPtsBuf.SetSize(nSamples);
  filterMatchDistBuf.SetSize(nSamples);

  DImax = D0max;
  bFirstIter_Matches = true;
}

void cisstAlgorithmICP_RobustICP::ICP_UpdateParameters_PostMatch()
{
  // base class
  cisstAlgorithmICP_StdICP::ICP_UpdateParameters_PostMatch();

  // compute match distances
  for (unsigned int i = 0; i < nSamples; i++)
  {
    matchDist[i] = (samplePtsXfmd[i] - matchPts[i]).Norm();
  }

  // compute epsilon from the initial matches
  if (bFirstIter_Matches)
  {    
    epsilon = ComputeEpsilon(matchDist);
  }

  bFirstIter_Matches = false;
}

unsigned int cisstAlgorithmICP_RobustICP::ICP_FilterMatches()
{
  //
  // Outlier detection method:
  //  Zhang, "Ierative Point Matching for Registration of Free-Form Curves and Surfaces", IJCV 1994
  //

  //--- Round 1: remove points with distance > DI-1max ---//

  // filter match distances > DI-1max
  nFilteredSamples = 0;
  for (unsigned int i = 0; i < nSamples; i++)
  {
    if (matchDist[i] <= DImax)  // DImax is currently = DI-1max
    {
      filterSamplePtsBuf[nFilteredSamples] = samplePts[i];
      filterMatchPtsBuf[nFilteredSamples] = matchPts[i];
      filterMatchDistBuf[nFilteredSamples] = matchDist[i];
      nFilteredSamples++;
    }
  }
  filterSamplePts.SetRef(filterSamplePtsBuf, 0, nFilteredSamples);
  filterMatchPts.SetRef(filterMatchPtsBuf, 0, nFilteredSamples);
  filterMatchDist.SetRef(filterMatchDistBuf, 0, nFilteredSamples);

  //--- Round 2: remove points with distance > DImax ---//

  // compute mean and standard deviation of filtered match distances
  double sumDist = 0.0;
  double sumSqrDist = 0.0;
  for (unsigned int i = 0; i < nFilteredSamples; i++)
  {
    sumDist += filterMatchDist[i];
    sumSqrDist += filterMatchDist[i] * filterMatchDist[i];
  }
  distAvg = sumDist / (double)nFilteredSamples;
  distSD = sqrt(sumSqrDist / (double)nFilteredSamples + distAvg*distAvg);

  // update DImax
  if (distAvg < D)
  { // registration is quite good
    DImax = distAvg + 3.0*distSD;
  }
  else if (distAvg < 3.0*D)
  { // registration is still good
    DImax = distAvg + 2.0*distSD;
  }
  else if (distAvg < 6.0*D)
  { // registration is not too bad
    DImax = distAvg + distSD;
  }
  else
  { // registration is really bad
    DImax = epsilon;
  }

  // filter match distances > DImax
  nGoodSamples = 0;
  for (unsigned int i = 0; i < nFilteredSamples; i++)
  {
    if (filterMatchDist[i] <= DImax)
    {
      goodSamplePtsBuf[nGoodSamples] = filterSamplePts[i];
      goodMatchPtsBuf[nGoodSamples] = filterMatchPts[i];
      nGoodSamples++;
    }
  }
  goodSamplePts.SetRef(goodSamplePtsBuf, 0, nGoodSamples);
  goodMatchPts.SetRef(goodMatchPtsBuf, 0, nGoodSamples);

  nOutliers = nSamples - nGoodSamples;

  return nOutliers;
}

double cisstAlgorithmICP_RobustICP::ComputeEpsilon(vctDynamicVector<double> &sampleDist)
{
  unsigned int numSamps = sampleDist.size();

  double minDist = sampleDist.MinElement();
  double maxDist = sampleDist.MaxElement();

  unsigned int numBins = 16;
  double binWidth = (maxDist - minDist) / (double)numBins;

  // build histogram of match distances
  vctDynamicVector<unsigned int> bins(numBins, (unsigned int)0);
  unsigned int sampleBin;
  for (unsigned int i = 0; i < numSamps; i++)
  {
    if (sampleDist[i] == maxDist)
    { // handle max case
      sampleBin = numBins - 1;
    }
    else
    {
      sampleBin = (unsigned int)floor((sampleDist[i] - minDist) / binWidth);
    }

    bins(sampleBin)++;
  }

  // find histogram peak
  unsigned int peakBin = numBins;  // initialize to invalid bin
  unsigned int peakBinSize = 0;
  for (unsigned int i = 0; i < numBins; i++)
  {
    if (bins(i) >= peakBinSize)
    {
      peakBin = i;
      peakBinSize = bins(i);
    }
  }
  // find valley following peak
  //  (valley bin must be <= 60% of peak bin size)  
  double valleyThresh = 0.6 * (double)peakBinSize;
  unsigned int valleyBin = peakBin + 1;
  for (unsigned int i = peakBin + 1; i < numBins; i++)
  {
    if ((double)bins(i) <= valleyThresh)
    {
      break;
    }
    valleyBin = i + 1;
  }

  // set epsilon to the smallest distance in the valley bin
  double epsilon = minDist + valleyBin * binWidth;

  //printHistogram(bins, peakBin, valleyBin, minDist, maxDist, binWidth);

  return epsilon;
}

void cisstAlgorithmICP_RobustICP::printHistogram(
  vctDynamicVector<unsigned int> bins, 
  unsigned int peakBin, unsigned int valleyBin,
  double minDist, double maxDist,
  double binWidth)
{
  std::stringstream ss;
  
  unsigned int numBins = bins.size();

  ss << std::endl << "Match Distance Histogram:" << std::endl;
  ss << "minDist: " << minDist << std::endl;
  ss << "maxDist: " << maxDist << std::endl;
  ss << std::endl;

  // bin header
  ss << " bins | size | dist " << std::endl;
  for (unsigned int i = 0; i < numBins; i++)
  {
    ss << " "; ss.width(4); ss << i;
    ss << " | "; ss.width(4); ss << bins(i);
    ss << " | " << minDist + i*binWidth << " ";

    if (i == peakBin)
    {
      ss << " <--- peak bin";
    }
    if (i == valleyBin)
    {
      ss << " <--- valley bin";
    }

    ss << std::endl;
  }
  ss << std::endl;

  //// bin contents
  //ss << "      ";
  //for (unsigned int i = 0; i < numBins; i++)
  //{
  //  ss << "| "; ss.width(4); ss << bins(i) << " ";
  //}
  //ss << std::endl;

  std::cout << ss.str() << std::endl;

}