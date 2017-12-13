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
#ifndef _algICP_IMLP_h
#define _algICP_IMLP_h

#include "algICP.h"
#include "Ellipsoid_OBB_Intersection_Solver.h"
#include "utilities.h"
#include <limits.h>

class algICP_IMLP : public algICP, public algPDTree
{
  //
  // This class implements algorithms for Iterative Most Likely Point
  //  registration, which is a modification of ICP to compute optimal registration
  //  in the presence of anisotropic noise.
  // Key differences in this algorithm compared to ICP are:
  //  1) total least squares framework for P2P registration (find registration that maximizes the match likelihood)
  //  2) assigned correspondences based on most likely match rather than closest match
  //
  //   likelihood function:   Product_i[ (1/sqrt((2*pi)^3*|Mi|) * exp( (-1/2)*(Yi-R*Xi-t)'*inv(Mi)*(Yi-R*Xi-t) ) ]
  //       where Mi  =  covariance of error in Yi-R*Xi-t  (Mi = R*Mxi*R' + Myi)
  //             Mxi =  covariance of measurement error for source point xi
  //             Myi =  covariance of measurement error for target point yi
  //
  //   Match Error function:  E = log|Mi| + (Yi-R*Xi-t)'*inv(Mi)*(Yi-R*Xi-t)
  //   (derived from negative log probability that Yi corresponds to R*Xi+t given the noise model Mi)
  //
  //   NOTE: negative values of the match error function 
  //         are possible (in case of very small |Mi|)
  //


  //-- Algorithm Parameters --//

protected:

  Ellipsoid_OBB_Intersection_Solver IntersectionSolver;

  vctFrm3 FGuess; // intiial guess for registration

  bool bFirstIter_Matches;  // flag that the first iteration is being run

  // dynamic noise model
  double sigma2;      // match uncertainty (added to My covariances as sigma2*I)
  double sigma2Max;   // max threshold on match uncertainty
  vctDynamicVector<vct3>  residuals_PostMatch;    // Pclosest - Psample
  vctDoubleVec            sqrDist_PostMatch;      // ||Pclosest - Psample||^2

  // noise model
  // Note: Myi values are obtained from the mesh object
  vctDynamicVector<vct3x3>  Mxi;        // noise covariances of sample points
  vctDynamicVector<vct3>    eigMxi;     // eigenvalues of sample covariances
  vctDynamicVector<vct3x3>  R_Mxi_Rt;   // noise covariances of transformed sample points
  vctDynamicVector<vct3x3*> Myi;        // noise covariances of target correspondence points
  vctDynamicVector<vct3x3>  Myi_sigma2; // noise covariances of target correspondence points with match uncertainty added
  //vctDynamicVector<vct3x3>  Mi;       // noise covariances of match (R*Mxi*Rt + Myi)
  //vctDynamicVector<vct3x3>  inv_Mi;   // inverse noise covariances of match (R*Mxi*Rt + Myi)^-1
  //vctDynamicVector<double>  det_Mi;   // determinant of noise covariances of match |R*Mxi*Rt + Myi|
  //vctDynamicVector<double>  SqrMahalDist;  // square Mahalanobis distance of matches = (yi-Rxi-t)'*inv(Mi)*(yi-Rxi-t)

  // outlier handling
  unsigned int nOutliers;
  double ChiSquareThresh; // Chi Square threshold for the outlier test  
  double sumSqrDist_Inliers;
  vctDynamicVector<int>   outlierFlags;
  // measurement noise component of the sample noise model
  //  Note: if applying measurement noise to the target shape, then
  //        MsmtMyi should be used here as well
  vctDynamicVector<vct3x3> MsmtMxi;
  // covariance model for outlier tests
  //  (does not include planar noise model)
  vctDynamicVector<vct3x3> R_MsmtMxi_Rt;

  // algorithm-specific termination
  bool bTerminateAlgorithm;
  unsigned char costFuncIncBits;
  double costFuncValue, prevCostFuncValue, prevIncCostFuncValue;
  vctFrm3 Fdec;


  // Temporary Match Variables
  //
  // These variables are recomputed each time a new sample point 
  //  is used in a PD tree search
  // Variables for current sample point undergoing a match search
  //  these are set once for each sample point by the pre-match function
  //  and used whenever a tree search routine is called thereafter.
  //  These do not change relative to the node being searched => precomputing
  //  these avoids multiple recomputations.
  vct3x3 sample_RMxRt_sigma2;      // covariance of transformed sample point plus registration uncertainty term
  vct3   sample_RMxRt_sigma2_Eig;  // eigenvalues of sample_RMxRt_sigma2 in order of decreasing magnitude

  vct3x3 M;       // effective measurement error covariance for a node & sample pair
  vct3x3 N;       // decomposition of inv(M) = N'N
  double Dmin;    // inverse sqrt of largest eigenvalue of M
                  //  (or the sqrt of the smallest eigenvalue of inv(M))
  double MinLogM; // lower bound on the log component of error for this node


  //-- Algorithm Methods --//

public:

  // constructor
  algICP_IMLP(
    PDTreeBase *pTree,
    vctDynamicVector<vct3> &samplePts, 
    vctDynamicVector<vct3x3> &sampleCov,      // full noise model (measurement noise + surface model)
    vctDynamicVector<vct3x3> &sampleMsmtCov,  // partial noise model (measurement noise only)
    double outlierChiSquareThreshold = 7.81,
    double sigma2Max = std::numeric_limits<double>::max());

  // destructor
  virtual ~algICP_IMLP() {}

  virtual void  ComputeMatchStatistics(double &Avg, double &StdDev);

  // TODO: specify surface model independently and compute Mxi
  //       by adding surface and msmt covariances
  //  Mxi ~ full noise model (measurement noise + surface model)
  //  MsmtMxi ~ partial noise model (measurement noise only)
  virtual void  SetSamples(
    vctDynamicVector<vct3> &argSamplePts, 
    vctDynamicVector<vct3x3> &argMxi, 
    vctDynamicVector<vct3x3> &argMsmtMxi);

  //void SetSampleCovariances(vctDynamicVector<vct3x3> &Mi, vctDynamicVector<vct3x3> &MsmtMi);

  // Sets Chi Square threshold for the outlier test:
  // Note:    ChiSquare(0.6) = 2.95      
  //          ChiSquare(0.7) = 3.66      
  //          ChiSquare(0.8) = 4.64      
  //          ChiSquare(0.85) = 5.32     
  //          ChiSquare(0.9) = 6.25      
  //          ChiSquare(0.925) = 6.90    
  // default  ChiSquare(0.95) = 7.81     (1.96 Std Dev)
  //          ChiSquare(0.975) = 9.35    (2.24 Std Dev)
  //          ChiSquare(0.99) = 11.34    (2.56 Std Dev)
  //          ChiSquare(0.9973) = 14.16  (3.0 Std Dev)     MATLAB: chi2inv(0.9973,3)
  void SetChiSquareThreshold(double ChiSquareValue) { ChiSquareThresh = ChiSquareValue; }
  void SetSigma2Max(double sigma2MaxValue) { sigma2Max = sigma2MaxValue; }

protected:

  void UpdateNoiseModel_SamplesXfmd(vctFrm3 &Freg);

  void ComputeNodeMatchCov(PDTreeNode *node);

  void ComputeCovDecomposition_NonIter(const vct3x3 &M, vct3x3 &Minv, double &det_M);
  void ComputeCovDecomposition_NonIter(const vct3x3 &M, vct3x3 &Minv, vct3x3 &N, vct3x3 &Ninv, double &det_M);

  void ComputeCovDecomposition_SVD(const vct3x3 &M, vct3x3 &Minv, double &det_M);
  void ComputeCovDecomposition_SVD(const vct3x3 &M, vct3x3 &Minv, vct3x3 &N, vct3x3 &Ninv, double &det_M);

  bool IntersectionSphereFace(const vct3 &n,
    const vct3 &v0, const vct3 &v1,
    const vct3 &v2, const vct3 &v3,
    double radius, double sqrRadius);

  int FindVisibleEdges(double q0, double q1,
    double v00, double v01,
    double v10, double v11,
    double v20, double v21,
    double v30, double v31,
    int *vsblEdges,
    bool ccwSequence);

  inline bool EdgeIsVisible(double q0, double q1,
    double v0, double v1,
    double n0, double n1);

  inline double SquareDistanceToEdge(const vct3 &p, const vct3 &r);

  // virtual standard routines for matching
  void SamplePreMatch(unsigned int sampleIndex);


  //--- ICP Interface Methods ---//

public:

  virtual void    ICP_InitializeParameters(vctFrm3 &FGuess);
  virtual void    ICP_UpdateParameters_PostMatch();
  virtual void    ICP_UpdateParameters_PostRegister(vctFrm3 &Freg);

  virtual vctFrm3 ICP_RegisterMatches();
  virtual unsigned int ICP_FilterMatches();  

  virtual double  ICP_EvaluateErrorFunction();
  virtual bool    ICP_Terminate(vctFrm3 &Freg);

  //virtual void  ICP_ComputeMatches();
  //virtual std::vector<cisstICP::Callback> ICP_GetIterationCallbacks();


  //--- PD Tree Interface Methods ---//

  int  NodeMightBeCloser(
    const vct3 &v,
    PDTreeNode *node,
    double ErrorBound);

  virtual double FindClosestPointOnDatum(
    const vct3 &v,
    vct3 &closest,
    int datum) = 0;

  virtual int  DatumMightBeCloser(
    const vct3 &v,
    int datum,
    double ErrorBound) = 0;

};
#endif
