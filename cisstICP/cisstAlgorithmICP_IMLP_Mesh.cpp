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
#include "cisstAlgorithmICP_IMLP_Mesh.h"


// finds the point on this datum with lowest match error
//  and returns the match error and closest point
double cisstAlgorithmICP_IMLP_Mesh::FindClosestPointOnDatum( const vct3 &v,
                                                             vct3 &closest,
                                                             int datum )
{
  // first iteration variables that don't change
  static const vct3x3 I_7071(vct3x3::Eye()*0.7071); // sqrt(1/2)*I
  static const vct3x3 I1_4142(vct3x3::Eye()*1.4142); // sqrt(2)*I
  static const vct3x3 I2(vct3x3::Eye()*2.0); // 2*I
  static const vct3x3 I_5(vct3x3::Eye()*0.5); // 0.5*I

  static vct3 d,v0,v1,v2;
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
  //
  //   Mahalanobis Distance:  sqrt((x-v)'*Minv*(x-v))
  //
  //  Method
  //   1) Translation + Affine transform to spherical space converting
  //      Mahalanobis ellipsoid -> sphere and triangle -> triangle2
  //   2) Find closest point (c') on triangle' to the origin
  //      (using standard Euclidean means)
  //   3) Affine transform c' back to normal coordinates
  static vct3 p0,p1,p2,c;

  // 1: transform triangle to spherical coords
  pTree->TriangleVerticesCoords(datum, v0,v1,v2);
  p0 = N*(v0 - v);
  p1 = N*(v1 - v);
  p2 = N*(v2 - v);

  // 2: Find closest point on triangle to origin in spherical coords
  TCPS.FindClosestPointOnTriangle( vct3(0.0), p0,p1,p2, -1,c);

  // 3: transform closest point back to normal coords
  closest = Ninv*c + v;
#else
  // Find the closest point on this triangle in a Euclidean distance sense
  //   Euclidean Distance:   ||x-v||
  TCPS.FindClosestPointOnTriangle( v,
 			                             pTree->TriangleVertexCoord(datum,0),
 			                             pTree->TriangleVertexCoord(datum,1),
  	  	                           pTree->TriangleVertexCoord(datum,2),
	    	                           -1,closest);
#endif

  //vct3 d0,d1,d2;
  //d0 = closest - pTree->TriangleVertexCoord(datum,0);
  //d1 = closest - pTree->TriangleVertexCoord(datum,1);
  //d2 = closest - pTree->TriangleVertexCoord(datum,2);
  //if ((d0.NormSquare() > 1.0) && (d1.NormSquare() > 1.0) && (d2.NormSquare() > 1.0))
  //{
  //  std::cout << "Minv = [" << std::endl << Minv << std::endl << "]" << std::endl;
  //  std::cout << "v = [" << v << "]" << std::endl;
  //  std::cout << "p(1,:) = [" << pTree->TriangleVertexCoord(datum,0) << "]" << std::endl;
  //  std::cout << "p(2,:) = [" << pTree->TriangleVertexCoord(datum,1) << "]" << std::endl;
  //  std::cout << "p(3,:) = [" << pTree->TriangleVertexCoord(datum,2) << "]" << std::endl;
  //  std::cout << "c = [" << closest << "]" << std::endl;
  //}

  // return match error 
  //  log term plus square Mahalanobis distance
  d = (v-closest);
  return log(det_M) + vctDotProduct(d,Minv*d);
}


// fast check if a datum might have smaller match error than error bound
int cisstAlgorithmICP_IMLP_Mesh::DatumMightBeCloser( const vct3 &v,
                                                     int datum,
                                                     double ErrorBound)
{ 
  return true;
}
