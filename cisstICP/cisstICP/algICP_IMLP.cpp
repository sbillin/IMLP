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

#include "algICP_IMLP.h"
#include "cisstICP.h"
#include "PDTreeNode.h"
#include "RegisterP2P.h"
#include "utilities.h"

#include <limits.h>

#define COMPUTE_ERROR_FUNCTION
//#define DEBUG_IMLP
//#define REMOVE_OUTLIERS


algICP_IMLP::algICP_IMLP(
  PDTreeBase *pTree,
  vctDynamicVector<vct3> &samplePts, 
  vctDynamicVector<vct3x3> &sampleCov,
  vctDynamicVector<vct3x3> &sampleMsmtCov,
  double outlierChiSquareThreshold,
  double sigma2Max)
  : algICP(pTree, samplePts),
  algPDTree(pTree)
{
  SetSamples(samplePts, sampleCov, sampleMsmtCov);
  SetChiSquareThreshold(outlierChiSquareThreshold);
  SetSigma2Max(sigma2Max);
}


void algICP_IMLP::ComputeMatchStatistics(double &Avg, double &StdDev)
{
  // return the average mahalanobis distance of the matches
  //  based on point noise models only (measurement and surface model covariances)
  //  i.e. do not include sigma2
  double sumSqrMahalDist = 0.0;
  double sumMahalDist = 0.0;
  double sqrMahalDist;
  int nGoodSamples = 0;
  vct3x3 M, Minv;
  vct3 residual;
  for (unsigned int i = 0; i < nSamples; i++)
  {
    if (outlierFlags[i]) continue;  // skip outliers

    residual = matchPts[i] - Freg * samplePts[i];
    M = Freg.Rotation() * Mxi[i] * Freg.Rotation().Transpose() + *Myi[i];
    ComputeCovInverse_NonIter(M, Minv);
    sqrMahalDist = residual*Minv*residual;

    sumSqrMahalDist += sqrMahalDist;
    sumMahalDist += sqrt(sqrMahalDist);
    nGoodSamples++;
  }

  Avg = sumMahalDist / nGoodSamples;
  StdDev = (sumSqrMahalDist / nGoodSamples) + Avg*Avg;

  //// return the average match distance of the inliers
  //double matchDist = 0.0;
  //int nGoodSamples = 0;
  //for (unsigned int i = 0; i < nSamples; i++)
  //{
  //  if (outlierFlags[i]) continue;  // skip outliers

  //  matchDist += (matchPts[i] - Freg * samplePts[i]).Norm();
  //  nGoodSamples++;
  //}
  //return matchDist / nGoodSamples;
}

// TODO: change this so that covariances for measurement noise
//       and surface model are specified independently
//       rather than specifying measurement noise model
//       and the combination of measurement noise and surface models
// NOTE: Mxi includes both MsmtMxi (measurement noise model)
//       and planar surface approximation model
void algICP_IMLP::SetSamples(
  vctDynamicVector<vct3> &argSamplePts,
  vctDynamicVector<vct3x3> &argMxi,
  vctDynamicVector<vct3x3> &argMsmtMxi)
{
  if (argMxi.size() != nSamples || argMsmtMxi.size() != nSamples)
  {
    std::cout << "ERROR: number of covariances matrices does not match number of samples" << std::endl;
  }

  // base class
  algICP::SetSamples(argSamplePts);

  Mxi = argMxi;
  MsmtMxi = argMsmtMxi;

  eigMxi.SetSize(nSamples);
  for (unsigned int s = 0; s < nSamples; s++)
  {
    ComputeCovEigenValues_SVD(argMxi[s], eigMxi[s]);
  }

  R_Mxi_Rt.SetSize(nSamples);
  R_MsmtMxi_Rt.SetSize(nSamples);
  Myi_sigma2.SetSize(nSamples);
  Myi.SetSize(nSamples);
  
  outlierFlags.SetSize(nSamples);

  residuals_PostMatch.SetSize(nSamples);
  sqrDist_PostMatch.SetSize(nSamples);
}

//// TODO: change this so that covariances for measurement noise
////       and surface model are specified independently
////       rather than specifying measurement noise model
////       and the combination of measurement noise and surface models
//// NOTE: Mi includes both MsmtMi (measurement noise model)
////       and planar surface approximation model
//void algICP_IMLP::SetSampleCovariances(
//  vctDynamicVector<vct3x3> &Mi,
//  vctDynamicVector<vct3x3> &MsmtMi)
//{
//  if (Mi.size() != nSamples || MsmtMi.size() != nSamples)
//  {
//    std::cout << "ERROR: number of covariances matrices does not match number of samples" << std::endl;
//  }
//
//  Mxi.SetSize(nSamples);
//  eigMxi.SetSize(nSamples);
//  for (unsigned int s = 0; s<nSamples; s++)
//  {
//    Mxi[s] = Mi[s];
//    ComputeCovEigenValues_SVD(Mi[s], eigMxi[s]);
//  }
//
//  this->MsmtMxi = MsmtMi;
//}

void algICP_IMLP::ICP_InitializeParameters(vctFrm3 &FGuess)
{
  // initialize base class
  algICP::ICP_InitializeParameters(FGuess);
  this->FGuess = FGuess;

  bFirstIter_Matches = true;
  nOutliers = 0;

  bTerminateAlgorithm = false;
  costFuncIncBits = 0;
  costFuncValue = std::numeric_limits<double>::max();

  sigma2 = 0.0;  

  // begin with isotropic noise model for first match, 
  //  since we don't yet have an approximation for sigma2
  
  R_Mxi_Rt.SetAll(vct3x3::Eye());
  R_MsmtMxi_Rt.SetAll(vct3x3(0.0));
  
  Myi_sigma2.SetAll(vct3x3::Eye());  
  Myi.SetAll(NULL);

  outlierFlags.SetAll(0);

  //Mi.SetSize(nSamples);
  //inv_Mi.SetSize(nSamples);
  //det_Mi.SetSize(nSamples);
  //SqrMahalDist.SetSize(nSamples);

  if (nSamples != Mxi.size() || nSamples != eigMxi.size())
  {
    std::cout << " ======> ERROR: noise model array for sample points does not match number of samples" << std::endl
      << "  Did you forget to call algICP_IMLP::SetSampleCovariances() before starting ICP?" << std::endl;
    assert(0);
  }
}

void algICP_IMLP::ICP_UpdateParameters_PostMatch()
{
  // base class
  algICP::ICP_UpdateParameters_PostMatch();

  //// compute post match errors
  //ComputeErrors_PostMatch();

  // compute sum of square distances of inliers
  sumSqrDist_Inliers = 0.0;
  //double sqrDist;
  for (unsigned int s = 0; s < nSamples; s++)
  {
    residuals_PostMatch.Element(s) = samplePtsXfmd.Element(s) - matchPts.Element(s);
    sqrDist_PostMatch.Element(s) = residuals_PostMatch.Element(s).NormSquare();

    if (!outlierFlags[s])
    {
      sumSqrDist_Inliers += sqrDist_PostMatch.Element(s);
    }
  }

  // update the match uncertainty factor
  // Should this be divide by N or divide by 3*N?
  //  Ans: it should be divide by N, because we actually want the entire
  //       square distance to closest point as the possible variance along
  //       each axis in this case. (See Estepar, et al, "Robust Generalized Total Least Squares...")
  sigma2 = sumSqrDist_Inliers / (nSamples - nOutliers);
  //sigma2 = sumSqrDist_PostMatch / nSamples;
  
  // apply max threshold
  if (sigma2 > sigma2Max)
  {
    sigma2 = sigma2Max;
  }

  // update noise models of the matches
  for (unsigned int s = 0; s < nSamples; s++)
  {
    // update target covariances
    Myi[s] = algICP::pTree->DatumCovPtr(matchDatums.Element(s));   // use pointer here for efficiency
    // target covariance with match uncertainty
    Myi_sigma2[s] = *Myi[s];
    Myi_sigma2[s].Element(0, 0) += sigma2;
    Myi_sigma2[s].Element(1, 1) += sigma2;
    Myi_sigma2[s].Element(2, 2) += sigma2;

    //// match covariance
    //Mi.Element(s) = R_Mxi_Rt.Element(s) + Myi.Element(s);
    //// match covariance decomposition
    //ComputeCovDecomposition(Mi.Element(s), inv_Mi.Element(s), det_Mi.Element(s));
    //// match square Mahalanobis distance
    //SqrMahalDist.Element(s) = Residuals.Element(s)*inv_Mi.Element(s)*Residuals.Element(s);
  }

  if (bFirstIter_Matches)
  {
    // update R*Mx*Rt to use the sample measurement noise defined for Mx
    //   (must do this here since Mx was initialized to identiy for first
    //    match step, but for first registration step, we want to use the 
    //    actual measurement noise model)
    UpdateNoiseModel_SamplesXfmd(FGuess);
  }

  // Outlier Noise Model:
  //  update measurment noise models of the transformed sample points
  //  not including the surface model covariance
  vctRot3 R(FGuess.Rotation());
  for (unsigned int s = 0; s < nSamples; s++)
  {
    R_MsmtMxi_Rt[s] = R*MsmtMxi[s] * R.Transpose();
  }

#ifdef DEBUG_IMLP
  std::cout << "My0:" << std::endl << *Myi[0] << std::endl;
  std::cout << "My1:" << std::endl << *Myi[1] << std::endl;
#endif

  bFirstIter_Matches = false;
}

void algICP_IMLP::ICP_UpdateParameters_PostRegister(vctFrm3 &Freg)
{
  // base class
  algICP::ICP_UpdateParameters_PostRegister(Freg);

  UpdateNoiseModel_SamplesXfmd(Freg);
}

void algICP_IMLP::UpdateNoiseModel_SamplesXfmd(vctFrm3 &Freg)
{
  // update noise models of the transformed sample points
  static vctRot3 R;
  R = Freg.Rotation();
  for (unsigned int s = 0; s < nSamples; s++)
  {
    R_Mxi_Rt[s] = R*Mxi[s] * R.Transpose();
  }
#ifdef DEBUG_IMLP
  std::cout << "ComputeParameters_PostReg():" << std::endl
    << "Mx0: " << std::endl << Mxi[0] << std::endl
    << "Mx1: " << std::endl << Mxi[1] << std::endl
    << "R_Mx0_Rt: " << std::endl << R_Mxi_Rt[0] << std::endl
    << "R_Mx1_Rt: " << std::endl << R_Mxi_Rt[1] << std::endl;
#endif
}


double algICP_IMLP::ICP_EvaluateErrorFunction()
{

#ifdef COMPUTE_ERROR_FUNCTION

  // TODO: do something different if "REMOVE_OUTLIERS" enabled?

  //
  // Compute negative log-likelihood of match probabilities for error function
  //
  //   likelihood function:   Product_i[ (1/sqrt((2*pi)^3*|Mi|) * exp( (-1/2)*(Yi-R*Xi-t)'*inv(Mi)*(Yi-R*Xi-t) ) ]
  //       where Mi  =  covariance of error in Yi-R*Xi-t  (Mi = R*Mxi*R' + Myi)
  //             Mxi =  covariance of measurement error for source point xi
  //             Myi =  covariance of measurement error for target point yi
  //
  //  If C's are same:
  //  -loglik:  -Sum_i[ log(Ci) ] + (1/2)*Sum_i[ di'*inv(M)*di ]
  //       where Ci  =  1/(sqrt((2*pi)^3*|Mi|)
  //
  //  -loglik = (3/2)*N*log(2*pi) + (1/2)*Sum_i[log(Mi)] + (1/2)*Sum_i[ di'*inv(M)*di ]
  //
  //   Note:  Removing constant terms and constant factors, the error functions simplifies
  //          to the error below. If using an error tolerance as the termination criteria,
  //          the constant term could create a lot of bias at low magnitudes of the variable
  //          error component => it is best to not include the constant in the error.
  //          On the other hand, it may be better to include the constant term in the error
  //          if using error tolerance as a termination condition if the algorithm has
  //          a hard time terminating, as this would ensure termination happens below
  //          some threshold (otherwise the error could keep reducing by some fixed
  //          percentage even though its basically zero, and never terminate).
  //   Note:  It is best to avoid this problem entirely, and base termination on the change
  //          in the registration parameter estimates.
  //
  //   Simplified Error = Sum_i[log(Mi) + di'*inv(M)*di]
  //
  vctDynamicVector<vct3x3>  Mi(nSamples);           // noise covariances of match (R*Mxi*Rt + Myi)
  vctDynamicVector<vct3x3>  inv_Mi(nSamples);       // inverse noise covariances of match (R*Mxi*Rt + Myi)^-1
  vctDynamicVector<double>  det_Mi(nSamples);       // determinant of noise covariances of match |R*Mxi*Rt + Myi|
  vctDynamicVector<double>  SqrMahalDist(nSamples); // square Mahalanobis distance of matches = (yi-Rxi-t)'*inv(Mi)*(yi-Rxi-t)

  // compute mahalanobis distances of the matches
  vct3 residual;
  for (unsigned int s = 0; s < nSamples; s++)
  {
    residual = samplePtsXfmd.Element(s) - matchPts.Element(s);

    // match covariance
    Mi.Element(s) = R_Mxi_Rt.Element(s) + Myi_sigma2.Element(s);
    // match covariance decomposition
    ComputeCovDecomposition_NonIter(Mi.Element(s), inv_Mi.Element(s), det_Mi.Element(s));
    // match square Mahalanobis distance
    SqrMahalDist.Element(s) = residual*inv_Mi.Element(s)*residual;
  }

  //// This is not general enough since we want to computer a correct error
  ////  both post-match and post-registration steps; this method assumes
  ////  error evaluation always happens post-match
  //// compute mahalanobis distances of the matches
  //for (unsigned int s = 0; s < nSamples; s++)
  //{
  //  // match covariance
  //  Mi.Element(s) = R_Mxi_Rt.Element(s) + Myi_sigma2.Element(s);
  //  // match covariance decomposition
  //  ComputeCovDecomposition(Mi.Element(s), inv_Mi.Element(s), det_Mi.Element(s));
  //  // match square Mahalanobis distance
  //  SqrMahalDist.Element(s) = residuals_PostMatch.Element(s)*inv_Mi.Element(s)*residuals_PostMatch.Element(s);
  //}

  //-- Here We Compute the Full Negative Log-Likelihood --//

  static double nklog2PI = nSamples*3.0*log(2.0*cmnPI);
  double logCost = 0.0;
  double expCost = 0.0;
  for (unsigned int i = 0; i<nSamples; i++)
  {
    // Compute error contribution for this sample
    //  error: log(detM) + di'*inv(Mi)*di
    expCost += SqrMahalDist.Element(i);
    logCost += log(det_Mi.Element(i));
  }

  prevCostFuncValue = costFuncValue;
  costFuncValue = (nklog2PI + logCost + expCost) / 2.0;
  //double costFunctionValue = (logCost + expCost);


  //-- Test for algorithm-specific termination --//

  // remove last iteration from monitoring variable
  //  by bit shifting one bit to the right
  costFuncIncBits >>= 1;
  if (costFuncValue > prevCostFuncValue)
  {
    // set 4th bit in monitoring variable
    //  (since we want to monitor up to 4 iterations)
    costFuncIncBits |= 0x08;

    // signal termination if cost function increased another time within 
    //  the past 3 trials and if the value has not decreased since that time
    //  TODO: better to test if the value is not the same value as before?
    if (costFuncIncBits > 0x08 && abs(prevIncCostFuncValue - costFuncValue) < 1e-10)
    {
      bTerminateAlgorithm = true;
    }
    prevIncCostFuncValue = costFuncValue;
  }
  else 
  { // record last time the cost function was non-increasing
    Fdec = Freg;
  }

  return costFuncValue;

#else
  return 0.0;
#endif
}

bool algICP_IMLP::ICP_Terminate( vctFrm3 &F )
{
  if (bTerminateAlgorithm)
  {
    // return registration for last iteration of non-increasing cost
    F = Fdec;
    return true;
  }
  else
  {
    return false;
  }
}

// TODO: change &F to Freg
//  even better: return registration rather than pass-by-reference
vctFrm3 algICP_IMLP::ICP_RegisterMatches()
{
#ifndef REMOVE_OUTLIERS

  vctFrm3 dF;
  RegisterP2P_TLS(samplePtsXfmd, matchPts,
    R_Mxi_Rt, Myi_sigma2, dF);
  Freg = dF*Freg;

  return Freg;

#ifdef DEBUG_IMLP
  std::cout << "ComputeRegistration()" << std::endl;
  std::cout << "Freg: " << Freg << std::endl;
  //std::cout << "Var: " << pICP->meanSqrDist_PostMatch << std::endl;
  //std::cout << "My0:" << std::endl << Myi_sigma2[0] << std::endl;
  //std::cout << "My12:" << std::endl << Myi_sigma2[12] << std::endl;
  //std::cout << "R_Mx0_Rt:" << std::endl << R_Mxi_Rt[0] << std::endl;
  //std::cout << "R_Mx12_Rt:" << std::endl << R_Mxi_Rt[12] << std::endl;
#endif

#else

  // remove outliers from consideration completely
  unsigned int nGoodSamples = nSamples - nOutliers;
  vctDynamicVector<vct3>    goodSamplePtsXfmd(nGoodSamples);
  vctDynamicVector<vct3>    goodMatchPts(nGoodSamples);
  vctDynamicVector<vct3x3>  goodR_Mxi_Rt(nGoodSamples);
  vctDynamicVector<vct3x3>  goodMyi_sigma2(nGoodSamples);
  unsigned int goodIndex = 0;
  for (unsigned int i = 0; i < nSamples; i++)
  {
    if (outlierFlags[i]) continue;

    goodSamplePtsXfmd.at(goodIndex) = samplePtsXfmd[i];
    goodMatchPts.at(goodIndex) = matchPts[i];
    goodR_Mxi_Rt.at(goodIndex) = R_Mxi_Rt[i];
    goodMyi_sigma2.at(goodIndex) = Myi_sigma2[i];
    goodIndex++;
  }

  vctFrm3 dF;
  RegisterP2P_TLS(goodSamplePtsXfmd, goodMatchPts,
    goodR_Mxi_Rt, goodMyi_sigma2, dF);
  Freg = dF*Freg;

  return Freg;
#endif
}

unsigned int algICP_IMLP::ICP_FilterMatches()
{	
  //
  // Filer Matches for Outliers
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
  //  When an outlier is identified, increase the variance of its noise model
  //  such that residual for that match is considered to be w/in 1 standard 
  //  deviation of its mean. This will reduce the impact of this match error
  //  on the registration result.
  //

  double StdDevExpansionFactor = 3.0;    // std dev expansion factor
  double varExpansionFactor = StdDevExpansionFactor * StdDevExpansionFactor;

  nOutliers = 0;
  vct3x3 Mo, inv_Mo;
  double sqrMahalDist = 0.0;

  for (unsigned int s = 0; s < nSamples; s++)
  {
    // compute outlier noise model based on mearurment noise and sigma2 only
    //  and not including the surface model covariance
    //
    // Note: the code below assumes that the covariance model of the target
    //       shape is comprised of only a surface model covariance with zero
    //       measurement noise; if this is not true, then the target measurement
    //       noise should be added to the outlier covariance test below as well
    //   
    Mo = R_Mxi_Rt.Element(s);
    Mo.Element(0, 0) += sigma2;
    Mo.Element(1, 1) += sigma2;
    Mo.Element(2, 2) += sigma2;
    //Mo = R_Mxi_Rt.Element(s) + Myi_sigma2.Element(s);
    //Mo.Element(0, 0) += outlier_alpha;
    //Mo.Element(1, 1) += outlier_alpha;
    //Mo.Element(2, 2) += outlier_alpha;

    // compute Mahalanobis distance
    ComputeCovInverse_NonIter(Mo, inv_Mo);
    sqrMahalDist = residuals_PostMatch.Element(s)*inv_Mo*residuals_PostMatch.Element(s);

    // check if outlier
    if (sqrMahalDist > ChiSquareThresh)
    { // an outlier
      nOutliers++;
      outlierFlags[s] = 1;
      // add isotropic outlier term to noise model for this match
      //  with magnitude of half the square match distance
      double outlierScale = 0.5 * sqrDist_PostMatch.Element(s) * varExpansionFactor;
      Myi_sigma2[s].Element(0, 0) += outlierScale;
      Myi_sigma2[s].Element(1, 1) += outlierScale;
      Myi_sigma2[s].Element(2, 2) += outlierScale;
      R_Mxi_Rt[s].Element(0, 0) += outlierScale;
      R_Mxi_Rt[s].Element(1, 1) += outlierScale;
      R_Mxi_Rt[s].Element(2, 2) += outlierScale;

      // This shouldn't be done here, because alpha term is not used
      //  in error function
      //// For Error Function Evaluation:
      //// update match covariance
      //Mi.Element(s) = R_Mxi_Rt.Element(s) + Myi.Element(s);
      //// match covariance decomposition
      //ComputeCovDecomposition(Mi.Element(s), inv_Mi.Element(s), det_Mi.Element(s));
      //// match Mahalanobis distance
      //SqrMahalDist.Element(s) = Residuals.Element(s)*inv_Mi.Element(s)*Residuals.Element(s);
    }
    else
    {
      outlierFlags[s] = 0;
    }
  }

  return nOutliers;
}



// PD Tree Methods

void algICP_IMLP::SamplePreMatch( unsigned int sampleIndex )
{
  // measurement noise
  sample_RMxRt_sigma2 = R_Mxi_Rt[sampleIndex];

  // add term for match uncertainty here to the RMxiRt
  //  term rather than to each Myi term for improved efficiency
  sample_RMxRt_sigma2.Element(0,0) += sigma2;
  sample_RMxRt_sigma2.Element(1,1) += sigma2;
  sample_RMxRt_sigma2.Element(2,2) += sigma2;

  // compute eigenvalues
  //  Note: R does not change eigenvalues of Mxi => these may be precomputed
  //        also, adding isotropic terms adds same amount to each eigenvalue
  vct3 &sampEig = eigMxi[sampleIndex];
  sample_RMxRt_sigma2_Eig[0] = sampEig[0] + sigma2;
  sample_RMxRt_sigma2_Eig[1] = sampEig[1] + sigma2;
  sample_RMxRt_sigma2_Eig[2] = sampEig[2] + sigma2;

  #ifdef DEBUG_IMLP
    if (sampleIndex == 0)
    {
      std::cout << "SamplePreMatch():" << std::endl;
      std::cout << "sigma2: " << sigma2 << std::endl;
      std::cout << "R_Mx0_Rt: " << std::endl << R_Mxi_Rt[0] << std::endl;
      std::cout << "sample_RMxRt_sigma2: " << std::endl << sample_RMxRt_sigma2 << std::endl;
    }
  #endif
}

// fast check if a datum might have smaller match error than error bound
int algICP_IMLP::DatumMightBeCloser( const vct3 &v,
                                                int datum,
                                                double ErrorBound)
{
  return true;
}

// fast check if a node might contain a datum having smaller match error
//  than the error bound
int algICP_IMLP::NodeMightBeCloser( const vct3 &v,
                                               PDTreeNode *node,
                                               double ErrorBound )
{

// uncomment the desired node bounds check method
//  NODE_SIMPLE_ELLIPSOID_BOUNDS seems much faster
//#define NODE_SPHERE_BOUNDS
#define NODE_SIMPLE_ELLIPSOID_BOUNDS

#ifdef NODE_SPHERE_BOUNDS

  //// use local variables to enable different set of values on
  ////  first match
  //double MinLogM;
  double MEigMax;
  
  if (bFirstIter_Matches)
  { // isotropic noise model for first iteration
    MinLogM = 2.0794; // log(|I+I|) = log(2*2*2) = 2.0794
    MEigMax = 2;      // max eigenvalue of M = Mx + My = 2*I
  }
  else
  { // noise model defined by node for later iterations

    // update log bound for this node if it does not share a
    //  common log bound with its parent
    if (!node->bUseParentEigRankMinBounds)
    {
      // Compute min bound on log term for error
      //  a min bound on the determinant |RMxR'+My| is found by multiplying the 
      //  sum of each eigenvalue rank pair for RMxR' and the min eigenvalues of
      //  that rank within the node
      double r0,r1,r2;
      r0 = node->EigRankMin[0] + sample_RMxRt_sigma2_Eig[0];
      r1 = node->EigRankMin[1] + sample_RMxRt_sigma2_Eig[1];
      r2 = node->EigRankMin[2] + sample_RMxRt_sigma2_Eig[2];
      MinLogM = log(r0*r1*r2);
    }
    MEigMax = node->EigMax + sample_RMxRt_sigma2_Eig[0];
  }

  // subtract min bound on the log term component of the match error
  //  from the match error bound to get an upper bound on the Mahalanobis 
  //  part of the match error for finding a better match.
  double NodeErrorBound = ErrorBound - MinLogM;

  // Get the radius of a minimal bounding sphere around the source point
  //  that fully contains the Mahalanobis bound
  double maxSqrSearchDistance = MEigMax * NodeErrorBound;
  double maxSearchDistance = sqrt(maxSqrSearchDistance);

  // transform sample point into local coordinate system of node
  vct3 Fv = node->F*v;

  // Check if point lies w/in search range of the bounding box for this node
  //  Rather than comparing only the x-axis value, check all coordinate directions
  //   of the node bounding box to further refine whether this node may be within 
  //   the search range of this point. Using the node coordinate frame is still 
  //   useful in this context, because it ensures a node bounding box of minimum size.
  //  Conceptually, this check places another bounding box of size search distance
  //   around the point and then checks if this bounding box intersects the bounding
  //   box of this node.
  return node->Bounds.Includes(Fv, maxSearchDistance);
#endif

#ifdef NODE_SIMPLE_ELLIPSOID_BOUNDS

  //double eps = 0.001;
  //bool foundPoint = false;
  //if ( abs(v[0]-26.1953) < eps && abs(v[1]-21.1742) < eps && abs(v[2]-25.8686) < eps )
  //{
  //  foundPoint = true;
  //}
  //
  //bool foundNode = false;
  //if (foundPoint == true && NodeContainsDatum( 9237 ))
  //{
  //  double foundNode = true;
  //}

  // if node only contains few datums then check the datum directly
  //  since the time to check each datum is similar to the time to
  //  compute node intersection
  if (node->NData <= 3)
  { return 1; }

  // precompute static structure for efficiency
  static const vct3x3 I_7071(vct3x3::Eye()*0.7071); // sqrt(1/2)*I
  
  if (bFirstIter_Matches)
  { // isotropic noise model for first iteration
    // M = Mx + NodeEigMax*I = 2*I = V*S*V'  =>  S = 2*I, V = I
    // Minv = N'*N = I*(1/2)*I  =>  N = sqrt(1/2)*I = 0.7071*I
    N = I_7071;
    // Dmin:
    //  M = Mx + NodeEigMax*I = 2*I = V*S*V'  =>  S = 2*I, V = I
    //  Dinv[0] = sqrt(S[0]) = sqrt(2)
    //  Dmin = 1/Dinv[0] = 1/sqrt(2) = sqrt(1/2) = 0.7071
    Dmin = 0.7071;
    MinLogM = 2.0794; // log(|I+I|) = log(2*2*2) = 2.0794
  }
  else
  { // noise model defined by node

    // update Mahalanobis bound of this node if it does not
    //  share a common covariance bound with its parent
    if (!node->bUseParentEigMaxBound)
    {
      ComputeNodeMatchCov(node);
    }
    // update log bound for this node if it does not share a
    //  common log bound with its parent
    if (!node->bUseParentEigRankMinBounds)
    {
      // Compute min bound on log component of match error
      //  a min bound on the determinant |RMxR'+My| is found by multiplying the
      //  eigenvalues of RMxR' with the min node eigenvalues of each magnitude 
      //  rank in rank order
      double r0,r1,r2;
      r0 = node->EigRankMin[0] + sample_RMxRt_sigma2_Eig[0];
      r1 = node->EigRankMin[1] + sample_RMxRt_sigma2_Eig[1];
      r2 = node->EigRankMin[2] + sample_RMxRt_sigma2_Eig[2];
      MinLogM = log(r0*r1*r2);
    }
  }
          
  //static vct3 vtemp;
  //vct3 minCorner(-1.73518,     -1.28636,    -0.305777);
  //vct3 maxCorner(2.67119,      1.61454,     0.240069);
  //vct3 d1 = minCorner - node->Bounds.MinCorner;
  //vct3 d2 = maxCorner - node->Bounds.MaxCorner;
  //double eps = 0.001;
  //if ( abs(v[0]-26.1953) < eps && abs(v[1]-21.1742) < eps && abs(v[2]-25.8686) < eps 
  //  && d1.Norm() < eps && d2.Norm() < eps)
  //{
  //  std::cout << "NodeSearch:" << std::endl;
  //  std::cout << " NumDatums = " << node->NumData() << std::endl;
  //  std::cout << " Depth = " << node->myDepth << std::endl;
  //  std::cout << " NodeParent: " << node->Parent << std::endl;
  //  std::cout << " LEq Child: " << node->LEq << std::endl;
  //  std::cout << " More Child: " << node->More << std::endl;
  //  std::cout << " ErrorBound = " << ErrorBound << std::endl;
  //  std::cout << " N = [" << std::endl << N << std::endl << "]" << std::endl;
  //  std::cout << " M = [" << std::endl << M << std::endl << "]" << std::endl;
  //  std::cout << " Dmin = " << Dmin << std::endl;
  //  std::cout << " NodeEigMax = " << *node->pEigMax << std::endl;
  //  std::cout << " NodeEigRankMin = " << *node->pEigRankMin << std::endl;
  //  std::cout << " RMxRt_sigma2 = [" << std::endl << sample_RMxRt_sigma2 << std::endl << "]" << std::endl;
  //  std::cout << " RMxRt_sigma2_Eig = "<< sample_RMxRt_sigma2_Eig << std::endl;
  //  std::cout << " Size(Mxi) = " << Mxi.size() << std::endl;
  //  std::cout << " Mxi = [" << std::endl << Mxi.at(96) << std::endl << "]" << std::endl;
  //  std::cout << " RMxiRt = [" << std::endl << R_Mxi_Rt.at(96) << std::endl << "]" << std::endl;
  //  std::cout << " eigMxi = " << eigMxi.at(96) << std::endl;
  //  //std::cout << " Node Search Mxi(96): = [" << std::endl << Mxi.at(96) << std::endl << "]" << std::endl;
  //}

#ifdef DEBUG_IMLP
  if (node == node->pMyTree->Top 
    && v.X() == samplePtsXfmd[0].X()
    && v.Y() == samplePtsXfmd[0].Y()
    && v.Z() == samplePtsXfmd[0].Z())
  {
      std::cout << "NodeMightBeCloser()" << std::endl;
      std::cout << "EigMax: " << node->EigMax << std::endl;
      std::cout << "R_Mx0_Rt: " << std::endl << R_Mxi_Rt[0] << std::endl;
      std::cout << "R_Mx_Rt_sigma2: " << std::endl << sample_RMxRt_sigma2 << std::endl;
      std::cout << "M: " << std::endl << M << std::endl;
      std::cout << "N: " << std::endl << N << std::endl;
  }
#endif

  // subtract min bound of the log term component of node error
  //  from the current best match error to get effective Mahalanobis bound
  //  for creating the boundary ellipse
  double NodeErrorBound = ErrorBound - MinLogM;

  // Test intersection between the ellipsoid and the oriented bounding
  //  box of the node
  return IntersectionSolver.Test_Ellipsoid_OBB_Intersection( v, node->Bounds, node->F, 
                                                             NodeErrorBound, N, Dmin );

#endif // NODE_SIMPLE_ELLIPSOID_BOUNDS
}


// Helper Methods

// Note: this function depends on the SamplePreMatch() function
//       to set the noise model of the current transformed sample
//       point before this function is called
void algICP_IMLP::ComputeNodeMatchCov( PDTreeNode *node )
{
  // Note:  This function is called when searching a node that is using 
  //        its own noise model rather than that of its parent node

  // Compute the effective noise model for this node, assuming the noise 
  //  model of the transformed sample point has already been computed

  // noise model of transformed sample
  M = sample_RMxRt_sigma2;
  // add the effective My for this node
  M.Element(0,0) += node->EigMax;
  M.Element(1,1) += node->EigMax;
  M.Element(2,2) += node->EigMax;

  // TODO: can this be done using only the eigen decomposition
  //       of RMxRt
  // Compute Decomposition of M
  //   M = V*S*V'
  vct3    eigenValues;
  vct3x3  eigenVectors;  
  ComputeCovEigenDecomposition_NonIter(M, eigenValues, eigenVectors);

  // Compute Decomposition of Minv = N'*N
  //   Minv = R*D^2*R' = N'*N     M = R*Dinv^2*R' => R' = V', Dinv = sqrt(S)
  //   N = D*R'      Ninv = R*inv(D)
  vct3 Dinv(
    sqrt(eigenValues[0]),
    sqrt(eigenValues[1]),
    sqrt(eigenValues[2])
    );
  N.Row(0) = eigenVectors.Column(0) / Dinv[0];
  N.Row(1) = eigenVectors.Column(1) / Dinv[1];
  N.Row(2) = eigenVectors.Column(2) / Dinv[2];
  //N.Row(0) = eigenVectors.TransposeRef().Row(0) / Dinv[0];
  //N.Row(1) = eigenVectors.TransposeRef().Row(1) / Dinv[1];
  //N.Row(2) = eigenVectors.TransposeRef().Row(2) / Dinv[2];
  Dmin = 1.0/Dinv[0]; // eigenvalues are arranged in order of decreasing magnitude
}


void algICP_IMLP::ComputeCovDecomposition_NonIter(const vct3x3 &M, vct3x3 &Minv, double &det_M)
{
  // Compute eigen decomposition of M
  //   M = V*diag(S)*V'
  vct3    eigenValues;
  vct3x3  eigenVectors;
  ComputeCovEigenDecomposition_NonIter(M, eigenValues, eigenVectors);

  // Compute Minv
  //   Minv = V*diag(1/S)*V'
  static vctFixedSizeMatrix<double, 3, 3, VCT_COL_MAJOR> V_Sinv;
  static vct3 Sinv;
  Sinv[0] = 1.0 / eigenValues[0];
  Sinv[1] = 1.0 / eigenValues[1];
  Sinv[2] = 1.0 / eigenValues[2];
  V_Sinv.Column(0) = eigenVectors.Column(0)*Sinv[0];
  V_Sinv.Column(1) = eigenVectors.Column(1)*Sinv[1];
  V_Sinv.Column(2) = eigenVectors.Column(2)*Sinv[2];
  Minv.Assign(V_Sinv * eigenVectors.TransposeRef());

  // compute determinant of M
  det_M = eigenValues.ProductOfElements();
}

void algICP_IMLP::ComputeCovDecomposition_NonIter(const vct3x3 &M, vct3x3 &Minv, vct3x3 &N, vct3x3 &Ninv, double &det_M)
{
  // Compute eigen decomposition of M
  //   M = V*diag(S)*V'
  vct3    eigenValues;
  vct3x3  eigenVectors;
  ComputeCovEigenDecomposition_NonIter(M, eigenValues, eigenVectors);

  // Compute Minv
  //   Minv = V*diag(1/S)*V'
  static vctFixedSizeMatrix<double, 3, 3, VCT_COL_MAJOR> V_Sinv;
  static vct3 Sinv;
  Sinv[0] = 1.0 / eigenValues[0];
  Sinv[1] = 1.0 / eigenValues[1];
  Sinv[2] = 1.0 / eigenValues[2];
  V_Sinv.Column(0) = eigenVectors.Column(0)*Sinv[0];
  V_Sinv.Column(1) = eigenVectors.Column(1)*Sinv[1];
  V_Sinv.Column(2) = eigenVectors.Column(2)*Sinv[2];
  Minv.Assign(V_Sinv * eigenVectors.TransposeRef());

  // Compute Decomposition of Minv = N'*N
  //   Minv = R*D^2*R' = N'*N     M = R*Dinv^2*R' => R' = V', Dinv = sqrt(S)
  //   N = D*R'      Ninv = R*inv(D)
  vct3 Dinv(
    sqrt(eigenValues[0]),
    sqrt(eigenValues[1]),
    sqrt(eigenValues[2])
    );
  N.Row(0) = eigenVectors.Column(0) / Dinv[0];
  N.Row(1) = eigenVectors.Column(1) / Dinv[1];
  N.Row(2) = eigenVectors.Column(2) / Dinv[2];
  Ninv.Column(0) = eigenVectors.Column(0)*Dinv[0];
  Ninv.Column(1) = eigenVectors.Column(1)*Dinv[1];
  Ninv.Column(2) = eigenVectors.Column(2)*Dinv[2];
  
  // Compute determinant of M
  det_M = eigenValues.ProductOfElements();
}


void algICP_IMLP::ComputeCovDecomposition_SVD( const vct3x3 &M, vct3x3 &Minv, double &det_M )
{
  // Compute SVD of M
  static vctFixedSizeMatrix<double,3,3,VCT_COL_MAJOR> A;
  static vctFixedSizeMatrix<double,3,3,VCT_COL_MAJOR> U;
  static vctFixedSizeMatrix<double,3,3,VCT_COL_MAJOR> Vt;
  static vct3 S;
  static nmrSVDFixedSizeData<3,3,VCT_COL_MAJOR>::VectorTypeWorkspace workspace;
  try 
  {
    A.Assign(M);
    nmrSVD(A, U, S, Vt, workspace);
  }
  catch(...) 
  { assert(0); }

  // Compute Minv
  //   M = U*diag(S)*V'   where U = V
  //   Minv = V*diag(1/S)*U' = U*diag(1/S)*V'
  static vctFixedSizeMatrix<double,3,3,VCT_COL_MAJOR> Sinv_Ut;
  static vct3 Sinv;
  Sinv[0] = 1/S[0];
  Sinv[1] = 1/S[1];
  Sinv[2] = 1/S[2];
  Sinv_Ut.Row(0) = Sinv[0]*Vt.Row(0);
  Sinv_Ut.Row(1) = Sinv[1]*Vt.Row(1);
  Sinv_Ut.Row(2) = Sinv[2]*Vt.Row(2);
  Minv.Assign(U*Sinv_Ut);

  // Compute determinant of M
  det_M = S[0]*S[1]*S[2];
}

void algICP_IMLP::ComputeCovDecomposition_SVD( const vct3x3 &M, vct3x3 &Minv, 
                                                      vct3x3 &N, vct3x3 &Ninv, double &det_M )
{
  // Compute SVD of M
  static vctFixedSizeMatrix<double,3,3,VCT_COL_MAJOR> A;
  static vctFixedSizeMatrix<double,3,3,VCT_COL_MAJOR> U;
  static vctFixedSizeMatrix<double,3,3,VCT_COL_MAJOR> Vt;
  static vct3 S;
  static nmrSVDFixedSizeData<3,3,VCT_COL_MAJOR>::VectorTypeWorkspace workspace;
  try 
  {
    A.Assign(M);
    nmrSVD(A, U, S, Vt, workspace);
  }
  catch(...) 
  { assert(0); }

  // Compute Minv
  //   M = U*diag(S)*V'   where U = V
  //   Minv = V*diag(1/S)*U' = U*diag(1/S)*V'
  static vctFixedSizeMatrix<double,3,3,VCT_COL_MAJOR> Sinv_Ut;
  static vct3 Sinv;
  Sinv[0] = 1/S[0];
  Sinv[1] = 1/S[1];
  Sinv[2] = 1/S[2];
  Sinv_Ut.Row(0) = Sinv[0]*Vt.Row(0);
  Sinv_Ut.Row(1) = Sinv[1]*Vt.Row(1);
  Sinv_Ut.Row(2) = Sinv[2]*Vt.Row(2);
  Minv.Assign(U*Sinv_Ut);

  // Compute Decomposition of Minv = N'*N
  //   Minv = R*D^2*R' = N'*N
  //   N = D*R'
  //   Ninv = R*inv(D)
  static vct3 Dinv; //,D;
  Dinv[0] = sqrt(S[0]);
  Dinv[1] = sqrt(S[1]);
  Dinv[2] = sqrt(S[2]);
  //D[0] = 1/Dinv[0];
  //D[1] = 1/Dinv[1];
  //D[2] = 1/Dinv[2];
  //N.Row(0) = D[0]*Vt.Row(0);
  //N.Row(1) = D[1]*Vt.Row(1);
  //N.Row(2) = D[2]*Vt.Row(2);
  N.Row(0) = Vt.Row(0)/Dinv[0];
  N.Row(1) = Vt.Row(1)/Dinv[1];
  N.Row(2) = Vt.Row(2)/Dinv[2];
  Ninv.Column(0) = U.Column(0)*Dinv[0];
  Ninv.Column(1) = U.Column(1)*Dinv[1];
  Ninv.Column(2) = U.Column(2)*Dinv[2];
  
  // Compute determinant of M
  det_M = S[0]*S[1]*S[2];
}

