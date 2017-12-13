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
#include "algDirICP_StdICP.h"
#include "DirPDTreeNode.h"
#include "cisstICP.h"
#include "RegisterP2P.h"

#define EPS  1e-12


void algDirICP_StdICP::ICP_InitializeParameters(vctFrm3 &FGuess)
{
  // base class
  algDirICP::ICP_InitializeParameters(FGuess);

  errFuncPosWeight = 1.0;
  errFuncNormWeight = 0.0;
}


double algDirICP_StdICP::ICP_EvaluateErrorFunction()
{
  // Cost Function = RMS (root mean square) error of the non-outlier matches

  vct3 residual;
  double sumSqrDist = 0.0;

  for (unsigned int s = 0; s < nGoodSamples; s++)
  {
    residual = goodSamplePts.Element(s) - goodMatchPts.Element(s);
    sumSqrDist += residual.NormSquare();
  }

  return sqrt(sumSqrDist / nSamples);

  ////
  ////  Note: add thresholded outlier error as well for
  ////        smoother cost function
  ////
  //double SqrErr = gSumSqrDist_PostMatch + oSumSqrDist_PostMatch;
  //return sqrt(SqrErr / nSamples);
}


vctFrm3 algDirICP_StdICP::ICP_RegisterMatches()
{
  // non-weighted ICP
  RegisterP2P_LSQ(goodSamplePts, goodMatchPts, Freg);

  return Freg;

  // weighted ICP
  //return RegisterP2P_WLSQ( AllGoodSamples, 
  //                         AllGoodClosestPoints,
  //                         AllGoodSampleWeights,
  //                         Fact );
}


unsigned int algDirICP_StdICP::ICP_FilterMatches()
{

#if 1
  // Method 1: No Outlier Detection
  nGoodSamples = 0;
  nOutliers = 0;

  // use all samples
  for (unsigned int s = 0; s < nSamples; s++)
  { // copy all samples to buffer
    goodSamplePtsBuf.Element(s) = samplePts.Element(s);
    goodMatchPtsBuf.Element(s) = matchPts.Element(s);

    goodSampleNormsBuf.Element(s) = sampleNorms.Element(s);
    goodMatchNormsBuf.Element(s) = matchNorms.Element(s);
  }
  nGoodSamples = nSamples;

  // Non-destructively resize good sample reference vectors
  goodSamplePts.SetRef(goodSamplePtsBuf, 0, nGoodSamples);
  goodMatchPts.SetRef(goodMatchPtsBuf, 0, nGoodSamples);
  
  goodSampleNorms.SetRef(goodSampleNormsBuf, 0, nGoodSamples);
  goodMatchNorms.SetRef(goodMatchNormsBuf, 0, nGoodSamples);

  //gSumSqrDist_PostMatch = sumSqrDist_PostMatch;
  //oSumSqrDist_PostMatch = 0.0;

  return nOutliers;
#endif 

#if 0
  // Method 2: ChiSquare test using sigma2
  // Filer Matches for Outliers
  //  Outlier filtering based on Chi-Square test statistic on Mahalanobis distance of matches.
  //  Copy all non-outlier samples to the good samples vector.
  //
  // The Square Mahalanobis Distance of the matches follow a chi-square distribution
  //  with 3 degrees of freedom (1 DOF for each spatial dimension).
  //
  //  Detect outliers as:  Square Mahalanobis Distance > ChiSquare(c)
  //
  //  Note:  ChiSquare(0.95) = 7.81     (1.96 Std Dev)
  //         ChiSquare(0.975) = 9.35    (2.24 Std Dev)
  //         ChiSquare(0.99) = 11.34    (2.56 Std Dev)
  //         ChiSquare(0.9973) = 14.16  (3.0 Std Dev)     MATLAB: chi2inv(0.9973,3)
  //
  //  Outlier If:  SqrDist / sigma2 > ChiSquare(%)
  //
  //  Note: this implementation assumes a zero-mean residual.
  //

  // set the chi-square threshold
  //double ChiSquareThresh = 11.34;
  //double ChiSquareThresh = 14.16;
  double ChiSquareThresh = 9.35;

  // estimate the variance of match error from residuals
  double sigma2 = pICP->SumSqrDist_PostMatch / pICP->nSamplesTotal;
  double sqrDistThresh = sigma2*ChiSquareThresh;

  nOutliersTotal = 0;
  nGoodSamplesTotal = 0;
  gSumSqrDist_PostMatch = 0.0;
  oSumSqrDist_PostMatch = 0.0;

  // for all model/sample sets
  for (unsigned int set = 0; set < pICP->nSets; set++)
  {
    unsigned int nGood_Set = 0;
    unsigned int nSamps_Set = pICP->nSamplesSets[set];

    vctDynamicVectorRef<vct3>   Samples( pICP->SampleSets[set] );
    vctDynamicVectorRef<vct3>   ClosestPoints( pICP->ClosestPointSets[set] );    
    vctDynamicVectorRef<double> SquareDistances( pICP->SquareDistanceSets_PostMatch[set] );

    for (unsigned int s = 0; s < nSamps_Set; s++)
    {	
      // find outliers
      if ( SquareDistances.Element(s) > sqrDistThresh )
      { // an outlier
        nOutliersTotal++;
        // apply thresholded match error to outlier error
        oSumSqrDist_PostMatch += sqrDistThresh;
      }
      else
      {	// inlier
        AllGoodSamplesBuf.Element(nGoodSamplesTotal) = Samples.Element(s);          
        AllGoodClosestPointsBuf.Element(nGoodSamplesTotal) = ClosestPoints.Element(s);
        // apply match error to good samples error
        gSumSqrDist_PostMatch += SquareDistances.Element(s);
        nGoodSamplesTotal++;
        nGood_Set++;
      }
    }
    // set sample weights for this set
    for (unsigned int k = nGoodSamplesTotal-nGood_Set; k < nGoodSamplesTotal; k++)
    { // must know # of good samples in the set before setting the sample weights
      AllGoodSampleWeightsBuf.at(k) = pICP->Weights[set]/nGood_Set;
    }
    nOutliersSets[set] = nSamps_Set - nGood_Set;
  }

  // Non-destructively resize good sample reference vectors
  AllGoodSamples.SetRef( AllGoodSamplesBuf,0,nGoodSamplesTotal );
  AllGoodClosestPoints.SetRef( AllGoodClosestPointsBuf,0,nGoodSamplesTotal );
  AllGoodSampleWeights.SetRef( AllGoodSampleWeightsBuf,0,nGoodSamplesTotal );

  return nOutliersTotal;
#endif

#if 0
  // Method 3: filter matches having distance more than n*meanDistance
  // Discard Outliers
  //
  // Detect outliers based on standard deviations of positional error.  
  // Consider a point as outlier if point distance > # standard deviations
  //  from mean.
  //  
  //   outlier if:  abs(||Py - Px|| - meanDist) > #*sigma  
  //
  //   equivalent square distance threshold:  (||Py - Px|| - meanDist)^2 > (#^2)*sigma^2
  //
  // If outlierDistFactor<=0 then
  //		copies "Samples" to "GoodSamples" directly
  // else copies only non outlier points to "GoodSamples" and 
  //    nondestructively resizes "GoodSamples" vectors to the number 
  //    of points transferred
  //
  // TODO: for multiple sample set support, the outlier parameters (StdDev) and 
  //       the avg residuals should be computed seperately for each set
  //
  nGoodSamplesTotal = 0;
  nOutliersTotal = 0;

  gSumSqrDist_PostMatch = 0.0;
  oSumSqrDist_PostMatch = 0.0;

  double outlierFactor = 3.0;
  double outlierDist2Factor = outlierFactor*outlierFactor;
  vctDynamicVector<double>  SquaredDistFromMean(nSamples);


  // for all model/sample sets
  for (unsigned int set = 0; set < pICP->nSets; set++)
  {
    unsigned int nSamps_Set = pICP->nSamplesSets[set];
    unsigned int nGood_Set = 0;

    vctDynamicVectorRef<vct3>   Samples(pICP->SampleSets[set]);
    vctDynamicVectorRef<vct3>   ClosestPoints( pICP->ClosestPointSets[set] );    
    vctDynamicVectorRef<double> SquareDistances( pICP->SquareDistanceSets_PostMatch[set] );
    vctDynamicVectorRef<double> Distances( pICP->DistanceSets_PostMatch[set] );

    if (outlierFactor <= 0.0)
    {	// outlier filtering is disabled => use all samples
      for (unsigned int s = 0; s < nSamps_Set; s++)
      { // copy all samples to buffer
        AllGoodSamplesBuf.Element(nGoodSamplesTotal) = Samples.Element(s);
        AllGoodClosestPointsBuf.Element(nGoodSamplesTotal) = ClosestPoints.Element(s);
        AllGoodSampleWeightsBuf.Element(nGoodSamplesTotal) = pICP->Weights[set] / nSamps_Set;
        nGoodSamplesTotal++;
      }
      nGood_Set = nSamps_Set;
      // apply all errors to good sample errors
      gSumSqrDist_PostMatch += pICP->SumSquareDistanceSets_PostMatch[set];
    }
    else
    { // outlier filtering is enabled

      // compute variance of errors
      double posVar = 0.0;
      double temp;
      double meanDist = pICP->matchDistAvg_PostMatch;
      for (unsigned int s = 0; s < nSamps_Set; s++)
      {
        temp = Distances.Element(s) - meanDist;
        SquaredDistFromMean.Element(s) = temp*temp;
        posVar += SquaredDistFromMean.Element(s);
      }
      posVar /= nSamps_Set;

      // compute outlier threshold
      double sqrDistThresh = outlierDist2Factor * posVar;

      // filter outliers
      for (unsigned int s = 0; s < nSamps_Set; s++)
      {	// copy only good samples to buffer

        if ( SquaredDistFromMean.Element(s) > sqrDistThresh )
        { // outlier
          // apply thresholded match error to outlier error
          oSumSqrDist_PostMatch += sqrDistThresh;
        }
        else
        {	// inlier
          AllGoodSamplesBuf.Element(nGoodSamplesTotal) = Samples.Element(s);          
          AllGoodClosestPointsBuf.Element(nGoodSamplesTotal) = ClosestPoints.Element(s);
          // apply match error to good samples error
          gSumSqrDist_PostMatch += SquareDistances.Element(s);
          nGoodSamplesTotal++;
          nGood_Set++;
        }
      }
      // set sample weights for this set
      for (unsigned int k = nGoodSamplesTotal-nGood_Set; k < nGoodSamplesTotal; k++)
      { // must know # of good samples in the set before setting the sample weights
        AllGoodSampleWeightsBuf.at(k) = pICP->Weights[set]/nGood_Set;
      }
      nOutliersSets[set] = nSamps_Set - nGood_Set;
      nOutliersTotal += nSamps_Set - nGood_Set;
    }
  }

  // Non-destructively resize good sample reference vectors
  AllGoodSamples.SetRef( AllGoodSamplesBuf,0,nGoodSamplesTotal );
  AllGoodClosestPoints.SetRef( AllGoodClosestPointsBuf,0,nGoodSamplesTotal );
  AllGoodSampleWeights.SetRef( AllGoodSampleWeightsBuf,0,nGoodSamplesTotal );

  return nOutliersTotal;
#endif

}


// PD Tree Methods

int algDirICP_StdICP::NodeMightBeCloser(
  const vct3 &v, const vct3 &n,
  DirPDTreeNode const *node,
  double ErrorBound)
{
  vct3 Fv = node->F*v;          // transform point into local coordinate system of node

  // Check if point lies w/in search range of the bounding box for this node
  // Rather than comparing only the x-axis value, check all coordinate directions
  //  of the node bounding box to further refine whether this node may be within 
  //  the search range of this point. Using the node coordinate frame is still 
  //  useful in this context, because it ensures a node bounding box of minimum size.
  // Conceptually, this check places another bounding box of size search distance
  //  around the point and then checks if this bounding box intersects the bounding
  //  box of this node.
  return node->Bounds.Includes(Fv, ErrorBound);
}
