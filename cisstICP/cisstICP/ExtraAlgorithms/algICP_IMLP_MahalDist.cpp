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
#include "algICP_IMLP_MahalDist.h"
#include "cisstICP.h"
#include "PDTreeNode.h"
#include "RegisterP2P.h"
#include "Ellipsoid_OBB_Intersection_Solver.h"
#include "utilities.h"

#define COMPUTE_ERROR_FUNCTION
//#define DEBUG_IMLP
//#define REMOVE_OUTLIERS



double algICP_IMLP_MahalDist::ICP_EvaluateErrorFunction()
{

#ifdef COMPUTE_ERROR_FUNCTION

  // TODO: do something different if "REMOVE_OUTLIERS" enabled?

  //
  // Compute sum of square Mahalanobis Distances of Matches
  //
  vctDynamicVector<vct3x3>  Mi(nSamples);           // noise covariances of match (R*Mxi*Rt + Myi)
  vctDynamicVector<vct3x3>  inv_Mi(nSamples);       // inverse noise covariances of match (R*Mxi*Rt + Myi)^-1
  vctDynamicVector<double>  SqrMahalDist(nSamples); // square Mahalanobis distance of matches = (yi-Rxi-t)'*inv(Mi)*(yi-Rxi-t)

  // compute mahalanobis distances of the matches
  vct3 residual;
  double expCost = 0.0;
  for (unsigned int s = 0; s < nSamples; s++)
  {
    residual = samplePtsXfmd.Element(s) - matchPts.Element(s);

    // match covariance
    Mi.Element(s) = R_Mxi_Rt.Element(s) + Myi_sigma2.Element(s);
    // match covariance decomposition
    ComputeCovInverse_NonIter(Mi.Element(s), inv_Mi.Element(s));
    // match square Mahalanobis distance
    SqrMahalDist.Element(s) = residual*inv_Mi.Element(s)*residual;

    // sum error contribution for this sample
    //  error: di'*inv(Mi)*di
    expCost += SqrMahalDist.Element(s);
  }

  prevCostFuncValue = costFuncValue;
  costFuncValue = expCost;


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



// PD Tree Methods


// fast check if a node might contain a datum having smaller match error
//  than the error bound
int algICP_IMLP_MahalDist::NodeMightBeCloser( const vct3 &v,
                                               PDTreeNode *node,
                                               double ErrorBound )
{
#define NODE_SIMPLE_ELLIPSOID_BOUNDS
#ifdef NODE_SIMPLE_ELLIPSOID_BOUNDS

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
    //MinLogM = 2.0794; // log(|I+I|) = log(2*2*2) = 2.0794
  }
  else
  { // noise model defined by node

    // update Mahalanobis bound of this node if it does not
    //  share a common covariance bound with its parent
    //  (because the bound may reference the parent's value, we must store the 
    //   effective match covariances in the nodes rather than in the algorithm class)
    if (!node->bUseParentEigMaxBound)
    {
      ComputeNodeMatchCov(node);
    }
    //// update log bound for this node if it does not share a
    ////  common log bound with its parent
    //if (!node->bUseParentEigRankMinBounds)
    //{
    //  // Compute min bound on log component of match error
    //  //  a min bound on the determinant |RMxR'+My| is found by multiplying the
    //  //  eigenvalues of RMxR' with the min node eigenvalues of each magnitude 
    //  //  rank in rank order
    //  double r0,r1,r2;
    //  r0 = node->EigRankMin[0] + sample_RMxRt_sigma2_Eig[0];
    //  r1 = node->EigRankMin[1] + sample_RMxRt_sigma2_Eig[1];
    //  r2 = node->EigRankMin[2] + sample_RMxRt_sigma2_Eig[2];
    //  MinLogM = log(r0*r1*r2);
    //}
  }

  double NodeErrorBound = ErrorBound;
  //// subtract min bound of the log term component of node error
  ////  from the current best match error to get effective Mahalanobis bound
  ////  for creating the boundary ellipse
  //double NodeErrorBound = ErrorBound - MinLogM;
  

  // Test intersection between the ellipsoid and the oriented bounding
  //  box of the node
  return IntersectionSolver.Test_Ellipsoid_OBB_Intersection( v, node->Bounds, node->F, 
                                                             NodeErrorBound, N, Dmin );

#endif // NODE_SIMPLE_ELLIPSOID_BOUNDS
}
