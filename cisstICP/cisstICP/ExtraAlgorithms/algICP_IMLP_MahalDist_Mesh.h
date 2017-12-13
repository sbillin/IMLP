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
#ifndef _algICP_IMLP_MahalDist_Mesh_h
#define _algICP_IMLP_MahalDist_Mesh_h

#include "algICP_IMLP_MahalDist.h"
#include "PDTree_Mesh.h"
#include "TriangleClosestPointSolver.h"

class algICP_IMLP_MahalDist_Mesh : public algICP_IMLP_MahalDist
{ 
  //
  // This class implements the IMLP algorithm of the base class for 
  //  a mesh representation of the the target shape.
  //

  //-- Algorithm Parameters --//

protected:

  TriangleClosestPointSolver TCPS;
  PDTree_Mesh *pTree;
  cisstMesh *pMesh;


  //-- Algorithm Methods --//

public:

  // constructor
  algICP_IMLP_MahalDist_Mesh(
    PDTree_Mesh *pTree, 
    vctDynamicVector<vct3> &samplePts, 
    vctDynamicVector<vct3x3> &sampleCov,
    vctDynamicVector<vct3x3> &sampleMsmtCov,
    double outlierChiSquareThreshold = 7.81,
    double sigma2Max = std::numeric_limits<double>::max())
    : algICP_IMLP_MahalDist(pTree, samplePts, sampleCov, sampleMsmtCov, outlierChiSquareThreshold, sigma2Max),
    pTree(pTree),
    pMesh(pTree->MeshP),
	TCPS(*(pTree->MeshP))
  {};


  //-- PD Tree Methods --//

	double  FindClosestPointOnDatum( const vct3 &v, vct3 &closest, int datum);
	int     DatumMightBeCloser( const vct3 &v, int datum, double ErrorBound);
};
#endif
