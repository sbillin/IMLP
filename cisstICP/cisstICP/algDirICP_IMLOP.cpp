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

#include "algDirICP_IMLOP.h"
#include "DirPDTreeNode.h"
#include "RegisterP2P.h"
#include "utilities.h"

#define EPS  1e-12

void  algDirICP_IMLOP::SetNoiseModel(
  double initK, double initSigma2, double w_Rpos, bool dynamicallyEstParams)
{
  k_init = initK;
  sigma2_init = initSigma2;
  wRpos = w_Rpos;
  dynamicParamEst = dynamicallyEstParams;
}

void algDirICP_IMLOP::SetSamples(
  const vctDynamicVector<vct3> &argSamplePts,
  const vctDynamicVector<vct3> &argSampleNorms)
{
  // base class
  algDirICP::SetSamples(argSamplePts, argSampleNorms);

  std::cout << "algDirICP_IMLOP::SetSamples()" << std::endl;

  unsigned int nSamples = samplePts.size();

  // allocate buffers
  goodSamplePtsBuf.SetSize(nSamples);
  goodSampleNormsBuf.SetSize(nSamples);
  goodMatchPtsBuf.SetSize(nSamples);
  goodMatchNormsBuf.SetSize(nSamples);
}

double algDirICP_IMLOP::ICP_EvaluateErrorFunction()
{
  //// Return the negative log likelihood of the vonMises-Fisher
  ////  and Gaussian distributions under the assumption
  ////   of independence between the two distributions
  ////
  ////   Negative Log-Likelihood:
  ////    -log[ C * exp( k*dot(Ny,Nx) - B*||Y-X||^2 ) ]
  ////        where C  =  product of normalizations terms
  ////                    for Fisher and 3D Gaussian distributions
  ////              C  =  [k/(2*PI*(e^k-e^-k))]*[1/(2*PI*sigma^2)^(3/2)]
  ////              B  =  1/(2*sigma^2)
  //double logC;
  //static const double log2PI = log(2 * cmnPI); // compute this once for efficiency

#ifdef TEST_STD_ICP
  //// Test Standard ICP Condition:
  ////   compute the log of the normalization constant C for
  ////   only a Gaussian distribution
  ////   (k=0 and orientations are not being considered)
  //logC = -(3.0/2.0)*log(2*cmnPI*sigma2);

  // include match errors of good samples and of thresholded outliers
  //  to improve smoothness of cost function
  return B*(pICP->gSumSqrDist_PostMatch + pICP->oSumSqrDist_PostMatch); //- nSamples*logC;
#endif

  //// compute normalization constant
  //double logEkE_k = 0.0;
  //if (k >= 5)
  //{ // exp(k)>>exp(-k)  =>  exp(k)-exp(-k) ~= exp(k)
  //  //                  =>  log(exp(k)-exp(-k)) ~= log(exp(k)) = k
  //  logEkE_k = k;
  //}
  //else
  //{
  //  logEkE_k = log(exp(k) - exp(-k));
  //}
  //// to avoid overflow/underflow, calculate logC directly 
  ////  (rather than calculating C and then logC)
  //logC = log(k) - (5.0 / 2.0)*log2PI - logEkE_k - (3.0 / 2.0)*log(sigma2);

  vct3 residual;
  double sumSqrDist = 0.0;
  double sumNormProducts = 0.0;
  for (unsigned int s = 0; s < nGoodSamples; s++)
  {
    residual = goodSamplePts.Element(s) - goodMatchPts.Element(s);
    sumSqrDist += residual.NormSquare();

    sumNormProducts += vctDotProduct(goodSampleNorms.Element(s), goodMatchNorms.Element(s));
  }

  // include match errors of good samples and of thresholded outliers
  //  to improve smoothness of cost function
  // NOTE: adding an extra k*nSamples to the cost produces a cost function >= 0
  return B*(sumSqrDist)+k*(nSamples - sumNormProducts); // -logC*nSamples;

  //// including match errors of thresholded outliers helps improve
  ////  smoothness of cost function
  //return B*(gSumSqrDist_PostMatch + oSumSqrDist_PostMatch)
  //  - k*(gSumNormProducts_PostMatch + oSumNormProducts_PostMatch)
  //  - logC*nSamples;

  // This method for computing logC is unstable
  //   It produces C=0 => infinity for log(C) when k is large
  //   The problem comes when computing log(exp(k)-exp(-k))
  //   in that exp(k) blows up for big k
  ////  von Mises-Fisher normalizing constant
  //double Cp = k/(2*cmnPI*(exp(k)-exp(-k)));
  ////  Gaussian normalizing constant
  //double Cn = pow((2*cmnPI*sigma2),3/2);
  //C = Cp/Cn;
}

vctFrm3 algDirICP_IMLOP::ICP_RegisterMatches()
{
  RegisterP2P_Normals_vMFG(
    goodSamplePts,
    goodMatchPts,
    goodSampleNorms,
    goodMatchNorms,
    B, k, Freg);

  return Freg;
}


void algDirICP_IMLOP::ICP_InitializeParameters(vctFrm3 &FGuess)
{
  // initialize base class
  algDirICP::ICP_InitializeParameters(FGuess);

  k = k_init;
  sigma2 = sigma2_init;
  B = 1.0 / (2.0*sigma2_init);
  //B = 1.0;
  //sigma2 = 1/(2.0*B);
  //k = 1.0;

#ifdef TEST_STD_ICP
  k = 0.0;
#endif

  // monitoring variables
  errFuncNormWeight = k;
  errFuncPosWeight = B;
}

//void algDirICP_IMLOP::ICP_UpdateParameters_PostMatch()
//{
//  // base class
//  algDirICP::ICP_UpdateParameters_PostMatch();
//
//  //// update noise model
//  //UpdateNoiseModel(SumSqrDist_PostMatch, sumNormProducts_PostMatch);
//  //// update monitoring variables
//  //pICP->errFuncNormWeight = k;
//  //pICP->errFuncPosWeight = B;
//}

void algDirICP_IMLOP::ICP_UpdateParameters_PostRegister(vctFrm3 &Freg)
{
  // base class
  algDirICP::ICP_UpdateParameters_PostRegister(Freg);

  double sumSqrDist_PostRegister = 0.0;  
  double sumNormProducts_PostRegister = 0.0;
  for (unsigned int s = 0; s < nSamples; s++)
  { 
    sumSqrDist_PostRegister += 
      (samplePtsXfmd.Element(s) - matchPts.Element(s)).NormSquare();
    sumNormProducts_PostRegister += 
      vctDotProduct(sampleNormsXfmd.Element(s), matchNorms.Element(s));
  }

  //matchDistSD_PostRegister = 
  //  (sumSqrDist_PostRegister/nSamples) - matchDistAvg_PostRegister*matchDistAvg_PostRegister;
  //matchDegErrSD_PostRegister =
  //  (sumSqrDegErr_PostRegister / nSamples) - matchDegErrAvg_PostRegister*matchDegErrAvg_PostRegister;

  // update noise model
  UpdateNoiseModel(sumSqrDist_PostRegister, sumNormProducts_PostRegister);

  // update monitoring variables
  errFuncNormWeight = k;
  errFuncPosWeight = B;
}

unsigned int algDirICP_IMLOP::ICP_FilterMatches()
{

  // Method 1: No Outlier Detection
#if 1
  nGoodSamples = 0;
  nOutliers = 0;
  nPosOutliers = 0;
  nNormOutliers = 0;

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
  //gSumNormProducts_PostMatch = sumNormProducts_PostMatch;
  //oSumSqrDist_PostMatch = 0.0;
  //oSumNormProducts_PostMatch = 0.0;

  return nOutliers;
#endif 

  // Method 2: Chi-Square Outlier Test
#if 0
  // Filer Matches for Outliers
  //  
  // Positions:
  //  The Square Mahalanobis Distance of the matches follow a chi-square distribution
  //  with 3 degrees of freedom (1 DOF for each spatial dimension) assuming that the
  //  residuals of position are normally distributed about the mean. Since the
  //  positional noise model is isotropic, the square Mahalanobis distance reduces
  //  to the square Euclidean distance divided by the variance along one spatial dimension.
  //
  //  Detect position outliers as: ||Py - (R*Px+t)||^2 > ChiSquare(c) * Var(Pos)
  //
  //  Note:  ChiSquare(0.95) = 7.81     (1.96 Std Dev)
  //         ChiSquare(0.975) = 9.35    (2.24 Std Dev)
  //         ChiSquare(0.99) = 11.34    (2.56 Std Dev)
  //         ChiSquare(0.9973) = 14.16  (3.0 Std Dev)     MATLAB: chi2inv(0.9973,3)
  //
  //  Note:  Var(Pos) = (1/3*N)*Sum_i(||Pyi - (R*Pxi+t)||^2)
  //          - here we assume the residual error to have zero mean
  //          - we use 3*N as the denominator rather than N because
  //            the square distance is a sum of three independent Guassian
  //            RVs (one for each axis) and it is the variance along a single
  //            axis that appears in the Mahalanobis distance and Gaussian
  //            probability equations.
  //
  // Orientations:
  // Consider a point as outlier if orientation angle (acos(Na'Nb)) > 3 circular
  //  standard deviations
  //  
  //  Detect Orientation outlier as:  acos(dot(Nx,Ny)) > 3*circStdDev
  //
  //  Copies only non outlier points to "GoodsamplePts" and 
  //   nondestructively resizes "Good Sample" reference vectors 
  //   to the number of points transferred.
  //

  //double ChiSquareThresh = 11.34;
  double ChiSquareThresh = 14.16;
  double S2 = pICP->SumSqrDist_PostMatch / (3*nSamples); // variance estimate along a single axis
  double SquareDistThresh = S2 * ChiSquareThresh;

  double ThetaThresh = 3.0*pICP->circSD_PostMatch;
  ThetaThresh = ThetaThresh > cmnPI ? cmnPI : ThetaThresh;
  double NormProductThresh = cos(ThetaThresh);

  nGoodSamples = 0;
  nOutliers = 0;
  nPosOutliers = 0;
  nNormOutliers = 0;

  gSumSqrDist_PostMatch = 0.0;
  gSumNormProducts_PostMatch = 0.0;
  oSumSqrDist_PostMatch = 0.0;
  oSumNormProducts_PostMatch = 0.0;
  //gSumSqrNormDist_PostMatch = 0.0;
  //oSumSqrNormDist_PostMatch = 0.0;

  // for all model/sample sets
  for (unsigned int set = 0; set < pICP->nSets; set++)
  {
    unsigned int nSamps_Set = pICP->nsamplePtsSets[set];
    unsigned int nGood_Set = 0;

    vctDynamicVectorRef<vct3>   samplePts( pICP->SampleSets[set] );
    vctDynamicVectorRef<vct3>   sampleNorms( pICP->SampleNormSets[set] );
    vctDynamicVectorRef<vct3>   matchPts( pICP->ClosestPointSets[set] );
    vctDynamicVectorRef<vct3>   matchNorms( pICP->ClosestPointNormSets[set] );
    vctDynamicVectorRef<double> SquareDistances( pICP->SquareDistanceSets_PostMatch[set] );
    vctDynamicVectorRef<double> NormProducts( pICP->NormProductSets_PostMatch[set] );
    //vctDynamicVectorRef<double> Distances( pICP->DistanceSets_PostMatch[set] );
    //vctDynamicVectorRef<vct3>   TransformedsamplePts( pICP->TransformedSampleSets[set] );    
    //vctDynamicVectorRef<vct3>   TransformedsampleNorms( pICP->TransformedSampleNormSets[set] );
    //vctDynamicVectorRef<double> SquareDistances( pICP->SquareDistanceSets_PostMatch[set] );
    //vctDynamicVectorRef<vct3>   Residuals( pICP->ResidualSets_PostMatch[set] );
    //vctDynamicVectorRef<double> dThetas( pICP->dThetaSets_PostMatch[set] );

    // Filter Outliers
    for (unsigned int s = 0; s < nSamps_Set; s++)
    {
      // Positional Outlier Test
      if ( SquareDistances.Element(s) > SquareDistThresh )
      {	// a positional outlier
        nPosOutliers += 1;

        // Theshold the outlier error on position and (if applicable) orientation 
        //  so that the cost function remains smooth, but the outlier does not
        //  continue to increase the error.
        //
        // Dual error thresholds (one for position, one for orientation) are
        //  not ideal, as the error of an outlier may then
        //  change during its time as an outlier due to thresholding on
        //  both positional and orientation error. The alternative is a single
        //  threshold error value, but this would again create jumps in the error
        //  when adding/removing outliers.
        //
        // Either way the optimization runs correctly, what this effects is the
        //  error function value, which may be being monitored for the termination 
        //  condition. We have two choices:
        //   a) some change in outlier error affects the error function
        //   b) jumpy cost function in between optimizations, when adding/removing
        //      outliers, but cost of an outlier is always the same
        // Since a) is smooth and b) may be jumpy, we choose situation a).
        //
        double normProduct;
        //double sqrNormDist;
        // check if orientation is outlier as well
        //  if so, threshold orientation component of error as well
        // Note: smaller norm products mean larger error
        if ( NormProducts.Element(s) < NormProductThresh )
        { // orientation is outlier as well => threshold its error
          normProduct = NormProductThresh;
          //sqrNormDist = sqrNormDistThresh;
        }
        else
        {
          normProduct = NormProducts.Element(s);
          //sqrNormDist = SquareNormDistances.Element(s);
        }
        // apply match error to outlier error
        oSumSqrDist_PostMatch += SquareDistThresh;
        oSumNormProducts_PostMatch += normProduct;
        //oSumSqrNormDist_PostMatch += sqrNormDist; 
      }
      // Orientation Outlier Test
      else if ( NormProducts.Element(s) < NormProductThresh )
      { // an orientation outlier
        nNormOutliers += 1;

        // Threshold the outlier error on orientation so that the cost function is
        //  smooth, but the outlier does not continue to increase the error.
        //  (Don't have to threshold position error, because the position was
        //   already checked and cannot be an outlying value here)
        // apply match error to outlier error
        oSumSqrDist_PostMatch += SquareDistances.Element(s);
        oSumNormProducts_PostMatch += NormProductThresh;
        //oSumSqrNormDist_PostMatch += sqrNormDistThresh;
      }
      else
      { // inlier
        // copy this sample to the good sample buffers
        goodSamplePtsBuf.Element(nGoodSamples) = samplePts.Element(s);        
        goodSampleNormsBuf.Element(nGoodSamples) = sampleNorms.Element(s);
        goodMatchPtsBuf.Element(nGoodSamples) = matchPts.Element(s);
        goodMatchNormsBuf.Element(nGoodSamples) = matchNorms.Element(s);
        //AllGoodTransformedsamplePtsBuf.Element(pICP->nGoodSamples) = TransformedsamplePts.Element(s);
        //AllGoodTransformedsampleNormsBuf.Element(pICP->nGoodSamples) = TransformedsampleNorms.Element(s);
        //GoodSampleResiduals_PostMatch.Row(pICP->nGoodSamples).Assign(Residuals.Element(s));

        // apply match error to good sample error
        gSumSqrDist_PostMatch += SquareDistances.Element(s);
        gSumNormProducts_PostMatch += NormProducts.Element(s);
        //gSumSqrNormDist_PostMatch += SquareNormDistances.Element(s);

        nGoodSamples++;
        nGood_Set++;
      }
    }

    //// renormalize weights for good samples
    //for (unsigned int k = nGoodSamples-nGood_Set; k < nGoodSamples; k++)
    //{ // set sample weights
    //  // must know # of good samples in the set before setting the sample weights
    //  //  due to normalization factor
    //  AllGoodSampleWeightsBuf.at(k) = pICP->Weights[set]/nGood_Set;
    //}
    //nOutliersSets[set] = nSamps_Set - nGood_Set;
    nOutliers += nSamps_Set - nGood_Set;
  }

  goodSamplePts.SetRef( goodSamplePtsBuf,0,nGoodSamples );  
  goodSampleNorms.SetRef( goodSampleNormsBuf,0,nGoodSamples );  
  goodMatchPts.SetRef( goodMatchPtsBuf,0,nGoodSamples );
  goodMatchNorms.SetRef( goodMatchNormsBuf,0,nGoodSamples );
  //AllGoodSampleWeights.SetRef( AllGoodSampleWeightsBuf,0,nGoodSamples );
  //AllGoodTransformedsamplePts.SetRef( AllGoodTransformedsamplePtsBuf,0,nGoodSamples );
  //AllGoodTransformedsampleNorms.SetRef( AllGoodTransformedsampleNormsBuf,0,nGoodSamples );

  pICP->iterData->nOutliersPos = nPosOutliers;
  pICP->iterData->nOutliersNorm = nNormOutliers;

  return nOutliers;

  //if ( outlierDistFactor <= EPS && outlierNormFactor <= EPS )
  //{	// outlier filtering is disabled => use all samples
  //  for (unsigned int s = 0; s < nSamps_Set; s++)
  //  { // copy all samples to buffer
  //    goodSamplePtsBuf.Element(nGoodSamples) = samplePts.Element(s);
  //    goodSampleNormsBuf.Element(nGoodSamples) = sampleNorms.Element(s);
  //    goodMatchPtsBuf.Element(nGoodSamples) = matchPts.Element(s);        
  //    goodMatchNormsBuf.Element(nGoodSamples) = matchNorms.Element(s);
  //    AllGoodSampleWeightsBuf.Element(nGoodSamples) = pICP->Weights[set] / nSamps_Set;
  //    nGoodSamples++;
  //    //AllGoodTransformedsamplePtsBuf.Element(nGoodSamples) = TransformedsamplePts.Element(s);
  //    //AllGoodTransformedsampleNormsBuf.Element(nGoodSamples) = TransformedsampleNorms.Element(s);
  //  }
  //  nGood_Set = nSamps_Set;
  //  // apply all errors to good sample errors
  //  gSumSqrDist_PostMatch += pICP->SumSquareDistanceSets_PostMatch[set];
  //  gSumNormProducts_PostMatch += pICP->SumNormProductSets_PostMatch[set];
  //  //gSumSqrNormDist_PostMatch += SumSquareNormDistanceSets_PostMatch[set];      
  //}

  //// compute outlier thresholds for position and orientation
  //double sqrDistThresh = outlierDist2Factor * posVar;
  //double thetaThresh = outlierNormFactor * pICP->circSD_PostMatch;    
  //thetaThresh = thetaThresh > cmnPI ? cmnPI : thetaThresh;
  //double normProductThresh = cos(thetaThresh);
  ////// Use theta threshold to compute a threshold on the square geometric 
  //////  distance between normal vectors; to compute this, compare a unit
  //////  vector to itself rotated by theta threshold degrees
  //////  (used for thresholding outlier square norm distance error, not
  //////   used for detecting outlier)
  ////vct3 t1(1.0,0.0,0.0);
  ////vctAxAnRot3 AxAn(vct3(0.0,0.0,1.0),thetaThresh);
  ////vct3 t2 = vctRot3(AxAn)*t1;   // *** why must convert to vctRot3 type?
  ////double sqrNormDistThresh = (t1-t2).NormSquare();
#endif
}


// PD Tree Methods

int algDirICP_IMLOP::NodeMightBeCloser(
  const vct3 &v, const vct3 &n,
  DirPDTreeNode const *node,
  double ErrorBound)
{
  vct3 Fv = node->F*v;          // transform point into local coordinate system of node
  //vct3 Fn = F.Rotation()*n;   // don't need this since normal statistics for a node are calculated
  //  wrt the world reference frame, not local reference frame

  // Check if point lies w/in search range of the bounding box for this node
  //
  // Enlarge node boundary relative to the error bound:
  //
  //  cost:  k*(1-N'*Nclosest) + B*||v - closest||^2 
  //                                  (dist^2)
  //
  //  Improved Bound:  (Assume Best Match Normal is aligned by Max Angular Deviation
  //                    from the Mean Orientation of the Node)
  //
  //    If we know the avg normal (Navg) and the maximum angular deviation (dThetaMax)
  //     from the average normal for all triangles in a node, then:
  //     
  //     given: 0 < dThetaMax < 180 (deg)
  //     N'*Navg = cos(dThetaAvg)
  //      set ThetaC = dThetaAvg - dThetaMax    (assume avg normal tilted towards n by max deviation)
  //      if dThetaC < 0 then dThetaC = 0
  //     set N'*Nc = cos(dThetaC)    (this is the max possible norm product (lowest possible orienation error) 
  //                                  since N'*Nc <= cos(dThetaC) by reason of max deviation)
  //     
  //     =>  cost = k*(1-cos(dThetaC)) + B*dist^2
  //         maxSearchDist = sqrt([ErrorBound - k*(1-cos(dThetaC))]/B)
  //
  //     =>  search a node if the node boundary enlarged by
  //         maxSearchDist contains the point v
  //

  // Simple Bound
  //double searchDist2 = ErrorBound/posWeight;

  // Improved Bound
  double dThetaAvg = acos(n.DotProduct(node->Navg));
  double dThetaC = dThetaAvg - node->dThetaMax;
  dThetaC = dThetaC > 0.0 ? dThetaC : 0.0;   // enforce ThetaC >= 0
  double searchDist2 = (ErrorBound - k*(1 - cos(dThetaC))) / B;
  if (searchDist2 < 0.0)
  { // orientation error alone dismisses this node
    return 0;
  }

  // Rather than comparing only the x-axis value, check all coordinate directions
  //  of the node bounding box to further refine whether this node may be within 
  //  the search range of this point. Using the node coordinate frame is still 
  //  useful in this context, because it ensures a node bounding box of minimum size.
  // Conceptually, this check places another bounding box of size search distance
  //  around the point and then checks if this bounding box intersects the bounding
  //  box of this node.
  return node->Bounds.Includes(Fv, sqrt(searchDist2));


  // Search Method for Alternate Cost Function:
  // -----------------------------------------
  // Enlarge node boundary by the weighted normal distance
  //  to ensure all possible candidates are compared for cost function:
  //
  //  cost:  -k*N'*Nclosest + ||v - closest||^2 
  //
  //  Since the PD tree sorts only by position and not by the normal
  //  vector, the farthest positional search distance is obtained by assuming
  //  a triangle normal parallel to the point normal.  => cost = -k + dist^2
  //  
  //  Assuming DistBound is the smallest "cost" found so far, we have maximum
  //  search distance:  dist^2 = DistBound + k
  //                    dist = sqrt(DistBound + k)
  //
  //  NOTE: by this formulation, "cost" and "DistBound" may be negative, however,
  //        the search distance (after adding k) is always >= 0.
  //
}

// Helper Methods:

void algDirICP_IMLOP::UpdateNoiseModel(double sumSqrDist, double sumNormProducts)
{
  // Position Parameters
  // compute Gaussian variables
  //  B = 1/(2*sigma2)
  // divide sum of square distances by number of samples times 3 
  //  to get the estimate of variance because square distance is 
  //  a summation of three square Gaussian RVs, one for each coordinate axis
  double S2 = sumSqrDist / (3.0*nSamples);
  if (dynamicParamEst)
  {
    B = 1.0 / (2.0*S2);
    B = B > threshB ? threshB : B;  // threshold in case of perfect matches
    sigma2 = 1.0 / (2.0*B);
  }

  // Orientation Parameters
  // compute Fisher variables
  //
  //  MLE estimate for k is an implicit ratio of
  //   two Bessel functions => no analytical soln.
  //
  //  MLE for k:  set Ap(k) = R
  //              
  //    Ap(k) = Ip/2(k) / Ip/2-1(k)       <--- Bessel Functions
  //    p =  spatial dimension (i.e. 3)
  //    R =  Sum_i(dot(Ny,Nx))/N
  //
  //  Closed form approximation for k:  R(p-R^2)/(1-R^2)
  //
  //   NOTE: circular standard deviation of theta
  //         circSD = sqrt(-2*ln(R))
  //            where theta is angle between matched vectors
  //            i.e. theta = acos(dot(Ny,Nx))
  //

  // angular error of normal orientations
  ComputeCircErrorStatistics(sumNormProducts, Rnorm, circSD);

  // angular error of positions (wrt ctr of rotation)
  Rpos = ComputeRpos();

  // effective angular error
  // only include positional error if it reduces the value of K
  if (Rpos < Rnorm)
  {
    R = wRpos*Rpos + (1.0 - wRpos)*Rnorm;
  }
  else
  {
    R = Rnorm;
  }

  double R2 = R*R;
  if (dynamicParamEst)
  {
    // protect from division by zero
    if (R2 >= 1.0)
    {
      k = threshK;  // set k to its max value
    }
    else
    {
      k = R*(3.0 - R2) / (1.0 - R2);  // approx for k
      k = k > threshK ? threshK : k;  // perfect match threshold
    }

    // reduce k by a factor
    k = k * k_factor;
  }

#ifdef TEST_STD_ICP
  k = 0.0;
#endif
}


double algDirICP_IMLOP::ComputeRpos()
{

#define NMLZ_METHOD 1   // Normalization method

  // NOTE: could improve efficiency by computing this value as an 
  //       optional output argument in the vMFG P2P Registration
  //       function

  // Compute angular match error of sample positions relative 
  //  to the center of rotation; this is to include angular error
  //  of position matches measured from center of rotation of
  //  Procrustes problem as added indicator of rotational uncertainty

  vctDynamicVectorRef<vct3> S(samplePtsXfmd);
  vctDynamicVectorRef<vct3> C(matchPts);
  unsigned int N = S.size();
  vct3 Smean = vctCentroid(S);
  vct3 Cmean = vctCentroid(C);
  vct3 Sp, Cp;
  double dotProducts = 0.0;
  double nmlzTerm = 0.0;
  for (unsigned int i = 0; i < N; i++)
  {
    Sp.DifferenceOf(S[i], Smean);
    Cp.DifferenceOf(C[i], Cmean);

    //  Do not scale down to unit vectors, as points farther
    //    away should have greater effect on rotation.
    //  R calculation assumes unit vector normalization => need some
    //    form of normalization in the end.
    //  TODO: should normalization be a linear or square term?
#if (NMLZ_METHOD == 1)
    // use square normalization
    //  this is the best solution, as it works with the data directly
    //  i.e. the orientations do not have to be normalized prior to
    //       taking their dot product
    dotProducts += vctDotProduct(Sp, Cp);
    nmlzTerm += Sp.Norm()*Cp.Norm();
#elif (NMLZ_METHOD == 2)
    // use linear normalization
    double avgLen = (Sp.Norm() + Cp.Norm()) / 2.0;
    dotProducts += avgLen*vctDotProduct(Sp.Normalized(), Cp.Normalized());
    nmlzTerm += avgLen;
#elif (NMLZ_METHOD == 3)
    // use unit vectors
    dotProducts += vctDotProduct(Sp.Normalized(), Cp.Normalized());
#else
    std::cout << "Error: normalization method unrecognized" << std::endl;
#endif
  }
#if (NMLZ_METHOD == 3)
  nmlzTerm = N;
#endif

  double Rpos = dotProducts / nmlzTerm;
  // 0 <= Rpos <= 1
  Rpos = Rpos < 0.0 ? 0.0 : Rpos;
  Rpos = Rpos > 1.0 ? 1.0 : Rpos;
  //double posCircSD = sqrt(-2*log(Rpos));
  return Rpos;
}


// Below are various methods that were tried in an attempt to balance the k & B
//  weights s.t. k does not blow up too large

//switch (opt.errFunc)
//{
//case OptionsNorm::STDICP:
//  {
//    // orientation is not used in standard ICP => k=0
//    double sigma2 = this->S2;
//    k = 0.0;

//    // perfect match threshold
//    //  need to threshold both sigma2 & B, since sigma2 used for calculating
//    //  normalization constant.
//    //  Note:  if B > n, then sigma2 < 1/(2*n)
//    //B = B > 1e5 ? 1e5 : B;
//    sigma2 = sigma2 < 1.0/(2.0e5) ? 1.0/(2.0e5) : sigma2;
//    B = 1.0/(2.0*sigma2);

//    Compute_vMFGConstant( k,sigma2,logC );
//    err = ComputeErrorFunctionValue( B,k,logC,SumSquareDistances,SumNormProducts );
//    break;    
//  }
//case OptionsNorm::vMFG:
//  {
//    // update error function parameters
//    double sigma2;
//    Compute_vMFGParams( k,sigma2 );

//    // perfect match thresholds
//    //  need to threshold both sigma2 & B, since sigma2 used for calculating
//    //  normalization constant.
//    //  Note:  if B > n, then sigma2 < 1/(2*n)
//    //B = B > threshB ? threshB : B;
//    k = k > threshK ? threshK : k;
//    sigma2 = sigma2 < threshS2 ? threshS2 : sigma2;
//    B = 1.0/(2.0*sigma2);

//    // compute normalization constant for probability distribution
//    Compute_vMFGConstant( k,sigma2,logC );

//    // compute error function
//    err = ComputeErrorFunctionValue( B,k,logC,SumSquareDistances,SumNormProducts );
//    break;
//  }
//case OptionsNorm::kROTREG:
//  {
//    double sigma2 = this->S2;

//    // regularize k by dR
//    vctRodRot3 dR;
//    dR.From(iterData->dF.Rotation());   // convert to Rodrigues notation
//    double dAng = dR.Norm();
//    double cSD = this->circSD + dAng;
//    //  SD = sqrt(-2*log(R));
//    //   =>  R = exp(-(SD^2)/2)
//    double Rreg = exp(-(cSD*cSD)/2);
//    double Rreg2 = Rreg*Rreg;
//    k = Rreg*(3-Rreg2)/(1-Rreg2);

//    // perfect match thresholds
//    k = k > threshK ? threshK : k;
//    sigma2 = sigma2 < threshS2 ? threshS2 : sigma2;
//    B = 1.0/(2.0*sigma2);

//    Compute_vMFGConstant( k,sigma2,logC );
//    err = ComputeErrorFunctionValue( B,k,logC,SumSquareDistances,SumNormProducts );
//    break;
//  }
//case OptionsNorm::kREG:
//  {
//    double prevK = k; 
//    double sigma2;
//    Compute_vMFGParams( k,sigma2 );

//    // regularization on k
//    k = prevK + (k-prevK)*opt.kRegConst;

//    // perfect match thresholds
//    k = k > threshK ? threshK : k;
//    sigma2 = sigma2 < threshS2 ? threshS2 : sigma2;
//    B = 1.0/(2.0*sigma2);

//    Compute_vMFGConstant( k,sigma2,logC );
//    err = ComputeErrorFunctionValue( B,k,logC,SumSquareDistances,SumNormProducts );
//    break;
//  }
//case OptionsNorm::kbREG:
//  {
//    double prevK = k;
//    double prevB = B;

//    double sigma2;
//    Compute_vMFGParams( k,sigma2 );

//    // perfect match thresholds
//    k = k > threshK ? threshK : k;
//    sigma2 = sigma2 < threshS2 ? threshS2 : sigma2;
//    B = 1.0/(2.0*sigma2);

//    Compute_vMFGConstant( k,sigma2,logC );
//    err = ComputeErrorFunctionValue( B,k,logC,SumSquareDistances,SumNormProducts );
//    break;
//  }
//case OptionsNorm::kTHRESH:
//  {
//    double sigma2;
//    Compute_vMFGParams( k,sigma2 );

//    // Threshold K
//    //  Note: k = 600 equates to about 3.25 degrees circular standard deviation
//    k = k > opt.kThresh ? opt.kThresh : k;

//    // perfect match thresholds
//    sigma2 = sigma2 < threshS2 ? threshS2 : sigma2;
//    B = 1.0/(2.0*sigma2);

//    Compute_vMFGConstant( k,sigma2,logC );
//    err = ComputeErrorFunctionValue( B,k,logC,SumSquareDistances,SumNormProducts );
//    break;
//  }
//case OptionsNorm::kSCALE:
//  {
//    double sigma2;
//    Compute_vMFGParams( k,sigma2 );

//    // Scale K
//    k /= opt.kScale;

//    // perfect match thresholds
//    k = k > threshK ? threshK : k;
//    sigma2 = sigma2 < threshS2 ? threshS2 : sigma2;
//    B = 1.0/(2.0*sigma2);

//    Compute_vMFGConstant( k,sigma2,logC );
//    err = ComputeErrorFunctionValue( B,k,logC,SumSquareDistances,SumNormProducts );
//    break;
//  }
//case OptionsNorm::kFIXED:
//  {
//    k = opt.k;
//    double sigma2 = this->S2;

//    // perfect match thresholds
//    k = k > threshK ? threshK : k;
//    sigma2 = sigma2 < threshS2 ? threshS2 : sigma2;
//    B = 1.0/(2.0*sigma2);

//    Compute_vMFGConstant( k,sigma2,logC );      
//    err = ComputeErrorFunctionValue( B,k,logC,SumSquareDistances,SumNormProducts );
//    break;
//  }
//case OptionsNorm::kbFIXED:
//  {
//    k = opt.k;
//    B = opt.B;
//    double sigma2 = 1.0/(2.0*B);
//    Compute_vMFGConstant( k,sigma2,logC );
//    err = ComputeErrorFunctionValue( B,k,logC,SumSquareDistances,SumNormProducts );
//    break;
//  }
//case OptionsNorm::WRAPNORMAL:
//  {
//    // assign k according to the variance of a
//    //  wrapped normal distribution
//    double sigma2 = this->S2;
//    k = 1.0/(2.0*(this->circSD)*(this->circSD));  // ERROR: should be k = 1/circSD^2

//    // perfect match thresholds
//    k = k > threshK ? threshK : k;
//    sigma2 = sigma2 < threshS2 ? threshS2 : sigma2;
//    B = 1.0/(2.0*sigma2);

//    Compute_vMFGConstant( k,sigma2,logC );
//    err = ComputeErrorFunctionValue( B,k,logC,SumSquareDistances,SumNormProducts );
//    break;
//  }
//default:
//  {
//    std::cout << "No valid error function chosen" << std::endl;
//    assert(0);
//  }
//}