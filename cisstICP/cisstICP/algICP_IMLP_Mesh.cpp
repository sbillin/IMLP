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
#include "algICP_IMLP_Mesh.h"


// finds the point on this datum with lowest match error
//  and returns the match error and closest point
double algICP_IMLP_Mesh::FindClosestPointOnDatum( 
  const vct3 &point,
  vct3 &closest,
  int datum )
{
  // first iteration variables that don't change
  static const vct3x3 I_7071(vct3x3::Eye()*0.7071); // sqrt(1/2)*I
  static const vct3x3 I1_4142(vct3x3::Eye()*1.4142); // sqrt(2)*I
  static const vct3x3 I2(vct3x3::Eye()*2.0); // 2*I
  static const vct3x3 I_5(vct3x3::Eye()*0.5); // 0.5*I

  static vct3 d;
  static vct3x3 M,Minv,N,Ninv;
  double det_M;
  
  if (bFirstIter_Matches)
  { // isotropic noise model for first iteration

    // M = Mx + NodeEigMax*I = 2*I = V*S*V'  =>  S = 2*I, V = I
    // Minv = N'*N = I*(1/2)*I  =>  N = sqrt(1/2)*I = 0.7071*I
    M = I2;
    Minv = I_5;
    N = I_7071;
    Ninv = I1_4142;
    det_M = 8;
  }
  else
  { // Use noise model of node

    // compute noise model for this datum
    //  M = R*Mxi*R' + Myi
    M = sample_RMxRt_sigma2 + pMesh->TriangleCov[datum];
    ComputeCovDecomposition_NonIter(M,Minv,N,Ninv,det_M);
  }


#if 1
  // Find the closest point on this triangle in a Mahalanobis distance sense
  TCPS.FindMostLikelyPointOnTriangle(point, datum, N, Ninv, closest);
#else
  // Find the closest point on this triangle in a Euclidean distance sense
  //   Euclidean Distance:   ||x-v||
  TCPS.FindClosestPointOnTriangle( v, datum, closest );
#endif

  // return match error 
  //  log term plus square Mahalanobis distance
  d = (point-closest);
  return log(det_M) + vctDotProduct(d,Minv*d);
}


// fast check if a datum might have smaller match error than error bound
int algICP_IMLP_Mesh::DatumMightBeCloser( const vct3 &point,
                                                     int datum,
                                                     double ErrorBound)
{ 
  return true;
}
