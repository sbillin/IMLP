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

#include "algDirICP_GIMLOP_Mesh.h"


// PD Tree Methods

double algDirICP_GIMLOP_Mesh::FindClosestPointOnDatum(
  const vct3 &Xp, const vct3 &Xn,
  vct3 &closest, vct3 &closestNorm,
  int datum)
{
  // NOTE: noise model parameters are those set by the pre-match routine within the base class

#ifdef KENT_POS_ISOTROPIC
  // This finds closest point to this triangle in a Euclidean distance sense
  //   Euclidean Distance:   ||x-v||
  TCPS.FindClosestPointOnTriangle( Xp,
    pDirTree->TriangleVertexCoord(datum,0),
    pDirTree->TriangleVertexCoord(datum,1),
    pDirTree->TriangleVertexCoord(datum,2),
    -1,closest);
#else
  // Find the closest point to this triangle in a Mahalanobis distance sense
  //
  //   Mahalanobis Distance:  sqrt((x-v)'*Minv*(x-v))
  //
  TCPS.FindMostLikelyPointOnTriangle(Xp, datum, sampleN_Rt, sample_inv_N_Rt, closest);
#endif

  // norm has same value everywhere on this datum
  closestNorm = pDirTree->mesh.faceNormals[datum];

  // noise parameters set by the pre-match routine of the base class
  double major = vctDotProduct(sampleR_L.Column(0), closestNorm);
  double minor = vctDotProduct(sampleR_L.Column(1), closestNorm);

  //debugStream << "RL1*Xn: " << sampleR_L.Column(0) * Xn << std::endl
  //  << "RL2*Xn: " << sampleR_L.Column(1) * Xn << std::endl;


  // return the match error
  //  Note: add "extra" k so that match error is always >= 0
#ifdef KENT_POS_ISOTROPIC
  return sampleK * (1 - vctDotProduct(Xn, closestNorm)) - sampleB*(major*major - minor*minor) 
    + (Xp-closest).NormSquare()*sampleEmin/2.0;
#else
  return sampleK * (1 - vctDotProduct(Xn, closestNorm)) - sampleB*(major*major - minor*minor)
    + ((Xp - closest)*sampleR_InvM_Rt*(Xp - closest)) / 2.0;
#endif

}


int algDirICP_GIMLOP_Mesh::DatumMightBeCloser(
  const vct3 &Xp, const vct3 &Xn,
  int datum,
  double ErrorBound)
{
  // TODO: compute only orientation component here to check if orientation
  //       error alone excludes the point; store this error in the base
  //       class for re-use in the base class

  // doing a decent proximity check is complicated enough that it is
  //  better to just compute the full error directly
  return 1;
}